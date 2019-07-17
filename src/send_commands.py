#!/usr/bin/python
""" Allow user to send arbitrary commands over the control socket.
"""
import socket
import struct
import time

#define MULTI_OBSINFO_PORT 53001
reader_group = ('224.3.29.71',20000)
writer_group = ('224.3.29.71',20001)

def send_char(char,multicast_group):

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.settimeout(0.2)
    ttl = struct.pack('b', 1)
    sock.setsockopt(socket.IPPROTO_IP, socket.IP_MULTICAST_TTL, ttl)
    sock.sendto(char,multicast_group)
    #sock.shutdown(socket.SHUT_RDWR)
    sock.close()

if __name__ == '__main__':
    print 'Enter a command at the line and it will be sent to all writer and process_baseband processes.'
    print 'Valid commands are Q=Quit, F=Fake_obs, G=Stop_fake_obs.'
    while True:
        s = raw_input()
        if len(s) < 1:
            continue
        c = s[0]
        if c == 'Q':
            print 'Sending Quit.'
        elif c == 'F':
            print 'Sending a fake ObservationDocument. [START]'
        elif c == 'G':
            print 'Sending a fake ObservationDocument. [FINISH]'
        else:
            print 'Did not understand command.'
            continue
        send_char(c,writer_group)
        time.sleep(0.01)
        send_char(c,reader_group)

