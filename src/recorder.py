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

    nsec = int(sys.argv[1])

    for i in range(nsec):
        t0 = time.time() # time in Unix epoch
        t1 = t0 + 0.1
        s = 'Recorded data segment number %02d.'%(i)
        print 't0=',t0,' t1=',t1
        t = struct.pack('dd128s',t0,t1,s)
        send_trigger(t)
        time.sleep(1)
