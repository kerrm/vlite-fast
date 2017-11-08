#!/usr/bin/python
"""
Gather candidates from heimdall, coincidence them, and issue trigger
for events above threshold.
"""

host = 'vlite-nrl'
HEIMDALL_PORT = 27555
TRIGGER_PORT = 27556

import socket
from collections import deque,defaultdict
from candidate import Candidate,coincidence

# combine events overlapping (multiple triggers) provided their total
# length doesn't exceed MAX_DUMP s
MAX_DUMP = 20

# set up a listening socket for heimdall server

def make_server (nmax=15):
    s = socket.socket (socket.AF_INET, socket.SOCK_STREAM)
    s.setsockopt (socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    s.bind ( ('vlite-nrl', HEIMDALL_PORT) )
    # at most nmax queued -- set to the no. of antennas
    s.listen (nmax)
    return s

def trigger(all_cands,snthresh=8,minbeam=2,wmax=0.01,dmmin=100):
    """ Go through coincidenced beams and determine if there is an event
    satisfying the trigger criteria.
    """
    triggers = []
    coincident_count = 0
    good_count = 0
    for cand in all_cands:
        nbeam = cand.beam_mask.sum() 
        coincident_count += nbeam > 1
        c1 = nbeam >= minbeam
        c2 = cand.width < wmax
        c3 = cand.dm > dmmin
        c4 = cand.sn > snthresh
        good_count += c2 and c3 and c4
        if (c1 and c2 and c3 and c4):
            triggers.append(cand)
    print 'len(all_cands)=%d'%(len(all_cands))
    print 'coincident_count = %d'%(coincident_count)
    print 'good_count = %d'%(good_count)
    print 'len(triggers)=%d'%(len(triggers))

    return triggers

if __name__ == '__main__':

    server_socket = make_server()
    utc_groups = dict()
    utc_sent_triggers = defaultdict(set)
    while (True):
        clientsocket, address = server_socket.accept ()
        print 'Received a connection from ', address
        payload = deque()
        while (True):
            msg = clientsocket.recv (4096)
            if len(msg) == 0:
                break
            payload.append(msg)
        lines = filter(lambda l: len(l) > 0,
            map(str.strip,''.join(payload).split('\n')))
        print lines[0]
        for line in lines[1:]:
            print Candidate(None,line)

        # do I want empty entries?
        if len(lines) == 1:
            continue

        # this is file start
        toks = lines[0].split()
        # NB at least right now this appears to be local time
        utc = toks[0]
        beam = int(toks[3]) - 1 # heimdall wants 0-base

        # check UTC for start of new observation
        if (utc not in utc_groups.keys()):
            utc_groups[utc] = defaultdict(deque)

        cgroups = utc_groups[utc]


        # add in Candidate objects to the appropriate beam
        cgroups[beam].extend((Candidate(None,l) for l in lines[1:]))
        if len(cgroups.keys()) < 2:
            #print 'Only one beam, skipping coincidence/triggering.'
            #continue
            pass

        # coincidence them
        all_cands = coincidence(cgroups.values())

        # get triggers
        sent_triggers = utc_sent_triggers[utc]
        triggers = set(trigger(all_cands,minbeam=1)).difference(sent_triggers)

        if len(triggers) == 0:
            continue

        # group up triggers within MAX_DUMP seconds
        # for now, just allow one trigger per processing block (~30s)
        # could loop over sent triggers if needed
        i0 = min((t.i0 for t in triggers))
        i1 = -1
        for trig in triggers:
            if (trig.i1 - i0)*trig.tsamp < MAX_DUMP:
                sent_triggers.add(trig)
                i1 = max(i1,trig.i1)

        # send a trigger based on active_utc, i0, i1        
        dump_offs = int(i1*trig.tsamp)
        dump_len = int( (i1-i0)*trig.tsamp ) + 1

        print 'Sending trigger for UTC %s with offset %d and length %d.'%(utc,dump_offs,dump_len)

