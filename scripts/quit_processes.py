#!/usr/bin/python
""" Multicast a signal to all writer and process_baseband process to 
shutdown.
"""
import socket
import struct
import time

#define MULTI_OBSINFO_PORT 53001
reader_group = ('224.3.29.71',20000)
writer_group = ('224.3.29.71',20001)

def send_quit(multicast_group):

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.settimeout(0.2)
    ttl = struct.pack('b', 1)
    sock.setsockopt(socket.IPPROTO_IP, socket.IP_MULTICAST_TTL, ttl)
    sock.sendto('Q',multicast_group)
    #sock.shutdown(socket.SHUT_RDWR)
    sock.close()

if __name__ == '__main__':
    send_quit(writer_group)
    time.sleep(1.0)
    send_quit(reader_group)

