import socket
import struct
import sys
import time
import xml.dom.minidom
from xml.etree import ElementTree

#define MULTI_OBSINFO_PORT 53001
multicast_group = '239.192.3.1'
server_address = ('', 53000)

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


def stringify_antdoc(data):
    e = ElementTree.fromstring(data)
    elems = e.getchildren()
    tags = [x.tag for x in elems]
    print tags
    return
    # get time in clock format UT for easy of reading
    eop = e.attrib['startTime']
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
    sec = int(time.time())
    file('/users/mkerr/mtk/logs/antprop/antrop_%d.xml'%(sec),'w').write(data)
            
    #print >>sys.stderr, 'received %s bytes from %s' % (len(data), address)
    dom = xml.dom.minidom.parseString(data)
    print dom.toprettyxml()
    try:
        print stringify_obsdoc(data)
    except Exception as e:
        print e

"""
<AntennaPropertyTable configuration="A" creation="58725.65837644676" datasetId="19A-107.sb36897542.eb37158321.58725.65837579861">                                                                 
<EopSet>                                                                                 
<eopday>                                                                         
<epoch>58723.0</epoch>                                                   
<tai_utc>37.0</tai_utc>                                                  
<ut1_utc>-0.153017</ut1_utc>                                             
<x_pole>0.21601</x_pole>                                                 
<y_pole>0.35796</y_pole>                                                 
</eopday>                                                                        
<eopday>                                                                         
<epoch>58724.0</epoch>                                                   
<tai_utc>37.0</tai_utc>                                                  
<ut1_utc>-0.152843</ut1_utc>                                             
<x_pole>0.2157</x_pole>                                                  
<y_pole>0.35616</y_pole>                                                 
</eopday>                                                                        
<eopday>                                                                         
<epoch>58725.0</epoch>                                                   
<tai_utc>37.0</tai_utc>                                                  
<ut1_utc>-0.152935</ut1_utc>                                             
<x_pole>0.21545</x_pole>                                                 
<y_pole>0.3544</y_pole>                                                  
</eopday>                                                                        
<eopday>                                                                         
<epoch>58726.0</epoch>                                                   
<tai_utc>37.0</tai_utc>                                                  
<ut1_utc>-0.153287</ut1_utc>                                             
<x_pole>0.21518</x_pole>                                                 
<y_pole>0.35272</y_pole>                                                 
</eopday>                                                                        
<eopday>                                                                         
<epoch>58727.0</epoch>                                                   
<tai_utc>37.0</tai_utc>                                                  
<ut1_utc>-0.153834</ut1_utc>                                             
<x_pole>0.21497</x_pole>                                                 
<y_pole>0.35106</y_pole>                                                 
</eopday>                                                                        
</EopSet>                                                                                
<AntennaProperties name="ea14">                                                          
<widarID>14</widarID>                                                            
<pad>E72</pad>                                                                   
<X>16724.5138</X>                                                                
<Y>-10408.1207</Y>                                                               
<Z>-7275.8817</Z>                                                                
<offset>-1.9E-12</offset>                                                        
</AntennaProperties>                                                                     
<AntennaProperties name="ea13">                                                          
<widarID>13</widarID>                                                            
<pad>W56</pad>                                                                   
<X>-12070.0011</X>                                                               
<Y>-635.6011</Y>                                                                 
<Z>-6329.9819</Z>                                                                
<offset>6.3E-12</offset>                                                         
</AntennaProperties>                                                                     
<AntennaProperties name="ea12">                                                          
<widarID>12</widarID>                                                            
<pad>E16</pad>                                                                   
<X>1259.2935</X>                                                                 
<Y>-795.4632</Y>                                                                 
<Z>-556.0902</Z>                                                                 
<offset>-3.4E-12</offset>                                                        
</AntennaProperties>                                                                     
<AntennaProperties name="ea11">                                                          
<widarID>11</widarID>                                                            
<pad>N32</pad>                                                                   
<X>404.3592</X>                                                                  
<Y>2630.084</Y>                                                                  
<Z>3885.6242</Z>                                                                 
<offset>1.18E-11</offset>                                                        
</AntennaProperties>                                                                     
<AntennaProperties name="ea10">                                                          
<widarID>10</widarID>                                                            
<pad>E08</pad>                                                                   
<X>383.4825</X>                                                                  
<Y>-241.8706</Y>                                                                 
<Z>-169.4508</Z>                                                                 
<offset>-2.6E-12</offset>                                                        
</AntennaProperties>                                                                     
<AntennaProperties name="ea28">                                                          
<widarID>28</widarID>                                                            
<pad>W08</pad>                                                                   
<X>-428.6825</X>                                                                 
<Y>-24.1537</Y>                                                                  
<Z>-223.392</Z>                                                                  
<offset>-7.0E-12</offset>                                                        
</AntennaProperties>                                                                     
<AntennaProperties name="ea27">                                                          
<widarID>27</widarID>                                                            
<pad>E32</pad>                                                                   
<X>4132.27</X>                                                                   
<Y>-2627.1754</Y>                                                                
<Z>-1816.9196</Z>                                                                
<offset>-2.0E-13</offset>                                                        
</AntennaProperties>                                                                     
<AntennaProperties name="ea26">                                                          
<widarID>26</widarID>                                                            
<pad>W72</pad>                                                                   
<X>-18571.9216</X>                                                               
<Y>-960.1587</Y>                                                                 
<Z>-9755.5078</Z>                                                                
<offset>-5.1E-12</offset>                                                        
</AntennaProperties>                                                                     
<AntennaProperties name="ea25">                                                          
<widarID>25</widarID>                                                            
<pad>N08</pad>                                                                   
<X>37.4541</X>                                                                   
<Y>243.6678</Y>                                                                  
<Z>360.0444</Z>                                                                  
<offset>-1.19E-11</offset>                                                       
</AntennaProperties>                                                                     
<AntennaProperties name="ea24">                                                          
<widarID>24</widarID>                                                            
<pad>W32</pad>                                                                   
<X>-4623.2423</X>                                                                
<Y>-252.588</Y>                                                                  
<Z>-2416.6965</Z>                                                                
<offset>5.0E-12</offset>                                                         
</AntennaProperties>                                                                     
<AntennaProperties name="ea23">                                                          
<widarID>23</widarID>                                                            
<pad>E24</pad>                                                                   
<X>2522.3102</X>                                                                 
<Y>-1603.8748</Y>                                                                
<Z>-1108.8848</Z>                                                                
<offset>-2.5E-12</offset>                                                        
</AntennaProperties>                                                                     
<AntennaProperties name="ea22">                                                          
<widarID>22</widarID>                                                            
<pad>W48</pad>                                                                   
<X>-9266.5161</X>                                                                
<Y>-493.6189</Y>                                                                 
<Z>-4854.8491</Z>                                                                
<offset>7.1E-12</offset>                                                         
</AntennaProperties>                                                                     
<AntennaProperties name="ea20">                                                          
<widarID>20</widarID>                                                            
<pad>N24</pad>                                                                   
<X>255.3233</X>                                                                  
<Y>1661.0954</Y>                                                                 
<Z>2454.5031</Z>                                                                 
<offset>1.37E-11</offset>                                                        
</AntennaProperties>                                                                     
<AntennaProperties name="ea09">                                                          
<widarID>9</widarID>                                                             
<pad>N16</pad>                                                                   
<X>123.4366</X>                                                                  
<Y>801.6189</Y>                                                                  
<Z>1182.1394</Z>                                                                 
<offset>-6.38E-11</offset>                                                       
</AntennaProperties>                                                                     
<AntennaProperties name="ea08">                                                          
<widarID>8</widarID>                                                             
<pad>W16</pad>                                                                   
<X>-1407.4663</X>                                                                
<Y>-77.486</Y>                                                                   
<Z>-735.192</Z>                                                                  
<offset>2.88E-11</offset>                                                        
</AntennaProperties>                                                                     
<AntennaProperties name="ea07">                                                          
<widarID>7</widarID>                                                             
<pad>E48</pad>                                                                   
<X>8291.3217</X>                                                                 
<Y>-5251.6189</Y>                                                                
<Z>-3654.6842</Z>                                                                
<offset>-8.8E-12</offset>                                                        
</AntennaProperties>                                                                     
<AntennaProperties name="ea06">                                                          
<widarID>6</widarID>                                                             
<pad>N72</pad>                                                                   
<X>1627.4704</X>                                                                 
<Y>10581.1453</Y>                                                                
<Z>15618.8606</Z>                                                                
<offset>2.4E-12</offset>                                                         
</AntennaProperties>                                                                     
<AntennaProperties name="ea05">                                                          
<widarID>5</widarID>                                                             
<pad>E56</pad>                                                                   
<X>10804.792</X>                                                                 
<Y>-6832.7449</Y>                                                                
<Z>-4767.4485</Z>                                                                
<offset>1.9E-12</offset>                                                         
</AntennaProperties>                                                                     
<AntennaProperties name="ea04">                                                          
<widarID>4</widarID>                                                             
<pad>N48</pad>                                                                   
<X>810.5123</X>                                                                  
<Y>5273.2899</Y>                                                                 
<Z>7791.9922</Z>                                                                 
<offset>2.96E-11</offset>                                                        
</AntennaProperties>                                                                     
<AntennaProperties name="ea03">                                                          
<widarID>3</widarID>                                                             
<pad>W64</pad>                                                                   
<X>-15176.1943</X>                                                               
<Y>-793.0198</Y>                                                                 
<Z>-7964.4554</Z>                                                                
<offset>-4.8E-12</offset>                                                        
</AntennaProperties>                                                                     
<AntennaProperties name="ea02">                                                          
<widarID>2</widarID>                                                             
<pad>W40</pad>                                                                   
<X>-6777.0613</X>                                                                
<Y>-360.7018</Y>                                                                 
<Z>-3550.9465</Z>                                                                
<offset>0.0</offset>                                                             
</AntennaProperties>                                                                     
<AntennaProperties name="ea01">                                                          
<widarID>1</widarID>                                                             
<pad>W24</pad>                                                                   
<X>-2823.3431</X>                                                                
<Y>-158.3194</Y>
<Z>-1472.1916</Z>
<offset>1.2E-12</offset>
</AntennaProperties>
<AntennaProperties name="ea19">
<widarID>19</widarID>
<pad>E40</pad>
<X>6060.4661</X>
<Y>-3851.9618</Y>
<Z>-2665.2164</Z>
<offset>-7.5E-12</offset>
</AntennaProperties>
<AntennaProperties name="ea18">
<widarID>18</widarID>
<pad>N64</pad>
<X>1329.7112</X>
<Y>8645.1221</Y>
<Z>12760.7159</Z>
<offset>-6.2E-12</offset>
</AntennaProperties>
<AntennaProperties name="ea17">
<widarID>17</widarID>
<pad>N40</pad>
<X>592.6474</X>
<Y>3856.1703</Y>
<Z>5698.9555</Z>
<offset>-6.4E-12</offset>
</AntennaProperties>
<AntennaProperties name="ea16">
<widarID>16</widarID>
<pad>E64</pad>
<X>13585.1885</X>
<Y>-8598.3687</Y>
<Z>-5990.5096</Z>
<offset>-1.61E-11</offset>
</AntennaProperties>
<AntennaProperties name="ea15">
<widarID>15</widarID>
<pad>N56</pad>
<X>1056.9944</X>
<Y>6873.369</Y>
<Z>10148.7743</Z>
<offset>4.1E-12</offset>
</AntennaProperties>
</AntennaPropertyTable>
"""
