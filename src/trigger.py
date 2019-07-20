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
import ctypes
import struct

import time
import calendar

# a ctypes implementation for reference
#class Trigger(ctypes.BigEndianStructure):
#    _pack = 1
#    _fields_ = [("t0",ctypes.c_int64),("t1",ctypes.c_int64), ("meta", ctypes.c_char*128)]
#t0c = ctypes.c_int64(t0)
#t1c = ctypes.c_int64(t1)
#meta = ctypes.create_string_buffer(s,size=128)
#t = Trigger(t0,t1,s)
#print ctypes.sizeof(t)

# combine events overlapping (multiple triggers) provided their total
# length doesn't exceed MAX_DUMP s
MAX_DUMP = 20
DM_DELAY = 4.15e-3*(0.320**-2-0.384**-2)

# set up a listening socket for heimdall server

def make_server (nmax=18):
    s = socket.socket (socket.AF_INET, socket.SOCK_STREAM)
    s.setsockopt (socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    s.bind ( ('vlite-nrl', HEIMDALL_PORT) )
    # at most nmax queued -- set to the no. of antennas
    s.listen (nmax)
    return s

def trigger(all_cands,snthresh=8,minbeam=3,wmax=0.01,dmmin=70):
    """ Go through coincidenced beams and determine if there is an event
    satisfying the trigger criteria.
    """
    triggers = []
    coincident_count = 0
    good_count = 0
    we = [1e100,-1e100]
    for cand in all_cands:
        nbeam = (cand.beam_mask>0).sum() 
        coincident_count += nbeam > 1
        if cand.width < we[0]:
            we[0] = cand.width
        if cand.width > we[1]:
            we[1] = cand.width
        c1 = nbeam >= minbeam
        c2 = cand.width < wmax
        c3 = cand.dm > dmmin
        c4 = cand.sn > snthresh
        good_count += c2 and c3 and c4
        if (c1 and c2 and c3 and c4):
            triggers.append(cand)
    print 'min/max width: ',we[0],we[1]
    print 'len(all_cands)=%d'%(len(all_cands))
    print 'coincident_count = %d'%(coincident_count)
    print 'good_count = %d'%(good_count)
    print 'len(triggers)=%d'%(len(triggers))

    return triggers

trigger_group = ('224.3.29.71',20003)

def send_trigger(trigger_struct):

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.settimeout(0.2)
    ttl = struct.pack('b', 1)
    sock.setsockopt(socket.IPPROTO_IP, socket.IP_MULTICAST_TTL, ttl)
    sock.sendto(trigger_struct,trigger_group)
    #sock.shutdown(socket.SHUT_RDWR)
    sock.close()

if __name__ == '__main__':

    server_socket = make_server()
    utc_groups = dict()
    utc_sent_triggers = defaultdict(set)
    output = file('/home/vlite-master/mtk/trigger_log.asc','a')
    while (True):
        clientsocket, address = server_socket.accept ()
        print 'Received a connection from %s:%s.\n'%(address[0],address[1])
        output.write('Received a connection from %s:%s.\n'%(address[0],address[1]))
        payload = deque()
        while (True):
            msg = clientsocket.recv (4096)
            if len(msg) == 0:
                break
            payload.append(msg)
        lines = filter(lambda l: len(l) > 0,
            map(str.strip,''.join(payload).split('\n')))
        output.write('\n'.join(lines))
        print 'Received %d new candidates.'%(len(lines)-1)
        #print lines[0]
        # add this temporarily, to allow reconstruction of candidates
        # that are received
        #for line in lines[1:]:
        #    print line
        #for line in lines[1:]:
        #    print Candidate(None,line)

        # do I want empty entries?
        if len(lines) == 2:
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
        cgroups[beam].extend((Candidate(None,l) for l in lines[2:]))
        if len(cgroups.keys()) < 2:
            #print 'Only one beam, skipping coincidence/triggering.'
            #continue
            pass

        print 'UTC',utc

        # coincidence them
        all_cands = coincidence(cgroups.values())
        if all_cands is None:
            continue

        # get triggers
        sent_triggers = utc_sent_triggers[utc]
        current_triggers = trigger(all_cands,snthresh=7.5,minbeam=2,wmax=0.5,dmmin=20)
        new_triggers = set(current_triggers).difference(sent_triggers)
        print 'new_triggers len: ',len(new_triggers) # DEBUG

        if len(new_triggers) == 0:
            continue

        for trig in new_triggers:
            print 'TRIGGERING ON CANDIDATE:',trig
            i0,i1 = trig.i0,trig.i1

            # send a trigger based on active_utc, i0, i1        
            dm_delay = trig.dm*DM_DELAY
            dump_offs = i0*trig.tsamp
            dump_len = (i1-i0)*trig.tsamp + dm_delay

            # TODO -- it would be nice to print out the latency between the candidate
            # peak time and the time the trigger is sent; it is 40-50 s with current 
            # gulp settings
            print 'Sending trigger for UTC %s with offset %d and length %.2f.'%(utc,dump_offs,dump_len)
            s = "Trigger at UTC %s + %d"%(utc,dump_offs)
            t = time.strptime(utc,'%Y-%m-%d-%H:%M:%S')
            # add in 100ms buffer in case heimdall isn't perfectly accurate!
            t0 = calendar.timegm(t) + dump_offs - 0.1
            t1 = t0 + dump_len + 0.2
            print 't0=',t0,' t1=',t1
            t = struct.pack('dd128s',t0,t1,s)
            send_trigger(t)
            sent_triggers.add(trig)

# TODO shut down socket on interrupt
