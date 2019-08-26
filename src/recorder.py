#!/usr/bin/python
"""
Enable recording of data by issuing periodic triggers.
"""

host = 'vlite-nrl'
HEIMDALL_PORT = 27555
TRIGGER_PORT = 27556

import socket
import struct
import sys
import time

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

    # works with two modes
    # python recorder.py 10  -- record 10 seconds of data, writing out 1s/s
    # python recorder.py 10 1 -- record 10s of data *AT ONCE* (burst)

    # NB one issue here is that if we select the current time, that 1s buffer
    # will still be having data written to it, so by default record from the
    # buffer 1s in the past.  For burst mode, do the same.

    nsec = int(sys.argv[1])
    if len(sys.argv) > 2:
        burst = True
    else:
        burst = False

    if burst:
        t0 = time.time()-1 # time in Unix epoch
        t1 = t0-nsec
        t0,t1 = t1,t0
        s = 'Burst mode record of %d seconds.'%nsec
        print 't0=',t0,' t1=',t1
        t = struct.pack('dd128s',t0,t1,s)
        send_trigger(t)
        sys.exit(0)

    for i in range(nsec):
        t0 = time.time()-1 # time in Unix epoch
        t1 = t0 + 1e-6
        s = 'Recorded data segment number %02d.'%(i)
        print 't0=',t0,' t1=',t1
        t = struct.pack('dd128s',t0,t1,s)
        send_trigger(t)
        time.sleep(1)
