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

# Tell the operating system to add the socket to the multicast group
# on all interfaces.
group = socket.inet_aton(multicast_group)
sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
mreq = struct.pack('4sL', group, socket.INADDR_ANY)
sock.setsockopt(socket.IPPROTO_IP, socket.IP_ADD_MEMBERSHIP, mreq)

# Bind to the server address
sock.bind(server_address)


def stringify_obsdoc(data):
    e = ElementTree.fromstring(data)
    elems = e.getchildren()
    tags = [x.tag for x in elems]
    # get time in clock format UT for easy of reading
    time = e.attrib['startTime']
    # convert to hhmmss
    fracday = float(time)-int(float(time))
    fracday *= 24
    hh = int(fracday)
    fracday -= hh
    fracday *= 60
    mm = int(fracday)
    fracday -= mm
    fracday *= 60
    ss = int(round(fracday))
    clock_string = '%02d:%02d:%02d'%(hh,mm,ss)
    s = '\n'.join((
    (" Observation:"),
    ("    datasetId = %s"%(e.attrib['datasetId'])),
    ("    configId = %s"%(e.attrib['configId'])),
    ("    startTime = %10.8f [%s]"%(float(e.attrib['startTime']),clock_string)),
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
    return s

# Receive/respond loop
while True:
    print >>sys.stderr, '\nwaiting to receive message'
    data, address = sock.recvfrom(8192)
            
    #print >>sys.stderr, 'received %s bytes from %s' % (len(data), address)
    #dom = xml.dom.minidom.parseString(data)
    #print dom.toprettyxml()
    try:
        print stringify_obsdoc(data)
    except Exception as e:
        print e

