import socket
import struct
import sys
import xml.dom.minidom
from xml.etree import ElementTree
import time
import argparse
import os

def setup_socket():

    multicast_group = '239.192.3.2'
    server_address = ('', 53001)

    # Create the socket
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

    # Bind to the server address
    sock.bind(server_address)

    # Tell the operating system to add the socket to the multicast group
    # on all interfaces.
    group = socket.inet_aton(multicast_group)
    mreq = struct.pack('4sL', group, socket.INADDR_ANY)
    sock.setsockopt(socket.IPPROTO_IP, socket.IP_ADD_MEMBERSHIP, mreq)

    return sock

def get_projid(data):
    e = ElementTree.fromstring(data)
    return e.attrib['datasetId'].split('.')[0]

def print_obsdoc(data):
    e = ElementTree.fromstring(data)
    elems = e.getchildren()
    tags = [x.tag for x in elems]
    s = '\n'.join((
    (" Observation:"),
    ("    datasetId = %s"%(e.attrib['datasetId'])),
    ("    configId = %s"%(e.attrib['configId'])),
    ("    startTime = %10.8f"%(float(e.attrib['startTime']))),
    ("    name = %s"%(elems[tags.index('name')].text)),
    ("    ra = %10.8f"%(float(elems[tags.index('ra')].text))),
    ("    dec = %10.8f"%(float(elems[tags.index('dec')].text))),
    #print("    dra = %10.8f\n", od->dra);
    #print("    ddec = %10.8f\n", od->ddec);
    #print("    azoffs = %10.8f\n", od->azoffs);
    #print("    eloffs = %10.8f\n", od->eloffs);
    ("    startLST = %10.8f"%(float(elems[tags.index('startLST')].text))),
    ("    scanNo = %s"%(elems[tags.index('scanNo')].text)),
    ("    subscanNo = %s"%(elems[tags.index('subscanNo')].text)),
    #print("    primaryBand = %s\n", od->primaryBand);
    #print("    usesPband = %d\n", od->usesPband);
    ))
    print s

def send_email(projid,finish=False):
    if finish:
        cmd = "echo '' | mail -s 'Project ID %s block now finished.' matthew.kerr@gmail.com"%projid
    else:
        cmd = "echo '' | mail -s 'Project ID %s block now beginning.' matthew.kerr@gmail.com"%projid
    os.system(cmd)

def main_loop(sock,projids):

    monitor_all_projids = projids[0] == 'all'
    if monitor_all_projids:
        projids = []

    sent_email = dict()
    for projid in projids:
        sent_email[projid] = False

    while True:
        print >>sys.stderr, '\nwaiting to receive message'
        data, address = sock.recvfrom(8192)
                
        print >>sys.stderr, 'received %s bytes from %s' % (len(data), address)
        print_obsdoc(data)

        if monitor_all_projids:
            projid = get_projid(data)
            print 'Adding %s to list of monitored projids.'%projid
            if projid not in projids:
                sent_email[projid] = False
                projids += [projid]

        for projid in projids:
            if projid in data:
                # signal end of schedule block
                if 'FINISH' in data:
                    send_email(projid,finish=True)
                    sent_email[projid] = False
                elif sent_email[projid]:
                    continue
                else:
                    send_email(projid,finish=False)
                    sent_email[projid] = True
                    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="monitor XML multicast for particular program ids and send an email alert.")
    parser.add_argument("projids",help="Project IDs to monitor.  Can include other keywords, too.  [Source name most likely.]  Enter 'all' to receive an email for all project ids.",type=str,nargs='+')
    #parser.add_argument('--cooldown',type=float,default=3600,help="Interval (s) to wait before sending an email about the same project id.")
    args = parser.parse_args()

    sock = setup_socket()
    main_loop(sock,args.projids)

