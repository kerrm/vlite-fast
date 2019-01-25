""" Search a log for a position with a given tolerance, and print
relevant information about the observation and recorded data."""

import sys
from collections import deque
from math import cos

if __name__ == '__main__':
    usage = """
Usage: search_log_pos logfile ra dec [tol]
    ra, dec, and tolerance all in radians.
    tol -- default = 0.5 degrees
"""

    try:
        logfile = sys.argv[1]
        ra0 = float(sys.argv[2])
        de0 = float(sys.argv[3])
    except:
        print usage
        sys.exit(0)
    try:
        tol = float(sys.argv[4])
    except:
        tol = 3.14/180*0.5
    tolsq = tol**2

    lines = filter(lambda l: len(l) > 0,
            map(str.strip,file(logfile).readlines()))
    for iline,line in enumerate(lines):

        if line.split()[0] == 'SIGPROC_FILE':
            ra = float(lines[iline-13].split()[1])
            de = float(lines[iline-12].split()[1])
            distsq = (cos(de)*(ra-ra0))**2 + (de-de0)**2
            if distsq < tolsq:
                print 'Found a matching entry:'
                print lines[iline-11] # NAME
                print lines[iline-13] # RA 
                print lines[iline-12] # DEC
                print lines[iline-11] # MJD
                print lines[iline-3] # UTC
                print lines[iline+0] # file info
                # find out where it was written to
                for iteriline in range(10):
                    if 'Wrote' in lines[iline+iteriline].split():
                        print lines[iline+iteriline] # recording info

"""
[2018-12-15-09:21:50] Received header with size 4096.
[2018-12-15-09:21:50] psrdada header:
STATIONID    23                     
NCHAN        1                      
BANDWIDTH    -64.000000             
CFREQ        352.000000             
NPOL         2                      
NBIT         8                      
TSAMP        0.007812               
RA           3.840734               
DEC          -0.003056              
NAME         GRB_080310             
SCANSTART    58467.685983           
UTC_START    2018-12-14-23:21:50    
[2018-12-15-09:21:50] Beginning new observation.

[2018-12-15-09:21:51] Starting trim at frame 8391 and second 58910.
[2018-12-15-09:21:52] Source GRB_080310 not on target list, disabling filterbank data.
[2018-12-15-09:21:52] Filterbank output disabled.  Would have written to /home/vlite-master/mtk/data/20181215_162152_muos_ea23.fil.
[2018-12-15-09:21:52] Writing no-RFI-excision filterbanks to /dev/null.
[2018-12-15-09:21:52] Writing RFI-excision filterbanks to /dev/null.
[2018-12-15-09:21:52] STATIONID    23                     
RA           3.840734               
DEC          -0.003056              
NAME         GRB_080310             
SCANSTART    58467.685983           
NCHAN        4096                   
BANDWIDTH    -41.936330             
CFREQ        340.978403             
NPOL         1                      
NBIT         2                      
TSAMP        781.250000             
UTC_START    2018-12-15-16:21:52    
VDIF_MJD     58467                  
VDIF_SEC     5051607712             
SIGPROC_FILE /home/vlite-master/mtk/data/20181215_162152_muos_ea23_kur.fil   
[2018-12-15-09:21:52] Starting sec=58912, thread=0
[2018-12-15-09:22:04] /home/vlite-master/mtk/bin/heimdall -nsamps_gulp 45000 -gpu_id 0 -dm 2 1000 -boxcar_max 64 -group_output -zap_chans 0 190 -zap_chans 3900 4096 -beam 23 -k 42 -coincidencer vlite-nrl:27555
[2018-12-15-09:33:33] Wrote 918.81 MB (701.00 s) to /dev/null
[2018-12-15-09:33:33] Proc Time...702.933
[2018-12-15-09:33:33] Waiting for DADA header.
"""
