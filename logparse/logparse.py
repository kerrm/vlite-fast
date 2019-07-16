#!/usr/bin/python
import re
import glob
import datetime as dt
import os
import time

ROOT_Hashirama='/home/shining/study/MS/vLITE/logs/logs/'
ROOT_Nrl='/home/vlite-master/mtk/logs/'
ROOT=ROOT_Hashirama
NE='/tmp/'
PROMFILE='one.prom'
NODENAME='vlite-difx12'
TIMESLEEP=5
GLOBEXP=ROOT+'*'+NODENAME+'_writer*.log'
TIMESTAMP_REGEX = r"^(\[\d{4}-\d{2}-\d{2}-\d{1,2}:\d{2}:\d{2}\])"
TIMESTAMP_SIG = "[%Y-%m-%d-%H:%M:%S]"

stime = lambda x : dt.datetime.strptime(x, TIMESTAMP_SIG).strftime('%s')

entry = re.compile(
        TIMESTAMP_REGEX +
        r"(.*\n)"
        )

def tailf(fline):
    fline.seek(0,os.SEEK_END)
    while True:
        line = fline.readline()
        if not line:
            time.sleep(TIMESLEEP)
            continue
        yield line

def parser(lines, match):
    log = []
    for xlin in tailf(lines):
        sre = match(xlin)
        if sre:
            yield sre.group(1),log
            log = []
            log.append(sre.group(2).strip())
        else:
            log.append(xlin.strip())
    yield sre.group(1),log


def getter(gEXP):
    fil = max(glob.iglob(gEXP), key=os.path.getctime)
    print fil
    fl = open(fil,'r')
    for ts,x in parser(fl, entry.match):
        # print x
        if "ObservationDocument:" in x:
            srcname = x[4].split(' = ')[1]
            ret = "#HELP obspar Observation parameters\n"
            ret = ret + "#TYPE obspar gauge\n"
            ret = ret + '{0} obspar {{"sourcename"="{1}";"nodename"="{2}"}} 1\n'.format(stime(ts),srcname,NODENAME)
            ofile.write(ret)
            ofile.flush()
        xfil = max(glob.iglob(gEXP), key=os.path.getctime)
        if xfil == fil:
            continue
        else:
            fl.close()
            fl = open(xfil,'r')
            fil = xfil

if __name__ == '__main__':
    ofile = open(NE+PROMFILE,'wa+')
    getter(GLOBEXP)
    ofile.close()
