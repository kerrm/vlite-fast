import socket
import struct
import sys
import xml.dom.minidom
from xml.etree import ElementTree

#define MULTI_OBSINFO_PORT 53001
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

def print_obsdoc(data):
    e = ElementTree.fromstring(data)
    elems = e.getchildren()
    tags = [x.tag for x in elems]
    s = '\n'.join((
    ("    datasetId = %s\n"%(e.attrib['datasetId'])),
    ("    configId = %s\n"%(e.attrib['configId'])),
    ("    startTime = %10.8f\n"%(float(e.attrib['startTime']))),
    ("    name = %s\n"%(elems[tags.index('name')].text)),
    ("    ra = %10.8f\n"%(float(elems[tags.index('ra')].text))),
    ("    dec = %10.8f\n"%(float(elems[tags.index('dec')].text))),
    #print("    dra = %10.8f\n", od->dra);
    #print("    ddec = %10.8f\n", od->ddec);
    #print("    azoffs = %10.8f\n", od->azoffs);
    #print("    eloffs = %10.8f\n", od->eloffs);
    ("    startLST = %10.8f\n"%(float(elems[tags.index('startLST')].text))),
    ("    scanNo = %s\n"%(elems[tags.index('scanNo')].text)),
    ("    subscanNo = %s\n"%(elems[tags.index('subscanNo')].text)),
    #print("    primaryBand = %s\n", od->primaryBand);
    #print("    usesPband = %d\n", od->usesPband);
    ))
    print s

# Receive/respond loop
while True:
    print >>sys.stderr, '\nwaiting to receive message'
    data, address = sock.recvfrom(8192)
            
    print >>sys.stderr, 'received %s bytes from %s' % (len(data), address)
    #dom = xml.dom.minidom.parseString(data)
    #print dom.toprettyxml()
    print_obsdoc(data)

