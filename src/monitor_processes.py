""" An adhoc tool for monitoring the status of various processes.

The primary purpose is to restart the telescope in the case one of the
processes hangs.  Similar functionality could also be used for determining
when observations have finished and then launch necessary post-processing.
"""

from __future__ import print_function
import argparse
import time
import glob
import os

import numpy as np

class Configuration(object):
    """ Configuration for a single antenna processing."""

    def __init__(self,confline):
        toks = confline.split()
        self.hostname = toks[0]
        self.iface = toks[1]
        self.gpu = toks[2]
        self.reader_port = toks[3]
        self.writer_port = toks[4]
        self.info_port = toks[5]
        self.baseband_dada_key = toks[6]
        self.fb_dada_key = toks[7]
        self.write_fb = toks[8]
        self.nbit = toks[9]

    def match_process_baseband(self,line):
        """ Match a line from a process_baseband log to determine if 
            they correspond.
        """
        pass


def parse_config(fname):
  
    """ Read a configuration file and return Configuration objects.
  
        Configuration file entries are of the form:

# hostname iface gpu reader_port writer_port info_port baseband_dada_key fb_dada_key write_fb nbit
vlite-difx1.evla.nrao.edu eth0 0 30001 30101 30201 40 42 2 2

"""

    with file(fname,'r') as f:
        lines = filter(lambda x: len(x) > 0 and x[0] != '#',
                    map(str.strip,f.readlines()))
    return map(Configuration,lines)

def find_logs(t0,after=True,close=True):
    """ Find log files matching the invocation time of this process. 
    
    By default, it is assumed that these logs will be generated after this
    process is started, and will be created relatively close in time.
    """

    logs = sorted(glob.glob('/home/vlite-master/mtk/logs/*process*'))

    # calculate timestamps in UTC for each file
    timestr = [os.path.split(log)[1][:15] for log in logs]
    timestamps = [time.strptime(ts,'%Y%m%d_%H%M%S') for ts in timestr]
    epoch_times = np.asarray([time.mktime(ts) for ts in timestamps])
    # mktime always makes a local time, convert to UTC
    epoch_times -= time.timezone

    # find maximum time offset and then select all logs that are within 60s
    amax = np.argmax(epoch_times-t0)
    if epoch_times[amax]-t0 < 0:
        # no logs occur before invocation
        return
    good_indices = np.argwhere((epoch_times[amax]-epoch_times) <  60)[:,0]
    return [logs[idx] for idx in good_indices]

if __name__ == '__main__':

    # save invocation time; NB file strings are in UTC
    t0 = time.time()

    parser = argparse.ArgumentParser(description="Monitor signal-processing for VLITE-Fast.")
    parser.add_argument("config",help="Configuration file.")
    args = parser.parse_args()

    config_fname = args.config.strip()
    print('Found configuration file {0}'.format(config_fname))
    configs = parse_config(config_fname)

    logs = find_logs(t0)

