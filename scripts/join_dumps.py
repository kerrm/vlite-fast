""" Rename files files contiguous in time for easier managment."""
import glob
import os
import datetime

fnames = sorted(glob.glob('/mnt/ssd/dumps/*.vdif'))
base_names = [os.path.split(fname)[1] for fname in fnames]
# split files by antenna
antennas = sorted(set([name.split('_')[2] for name in base_names]))
for antenna in antennas:
    ant_names = [fname for fname in base_names if antenna in fname]
    time_stamps = [fname[:15] for fname in base_names]
    dts = [datetime.datetime.strptime(time_stamp,'%Y%m%d_%H%M%S') for time_stamp in time_stamps]
    previous_dt = dts[0]
    current_ts = time_stamps[0]
    new_names = []
    counter = 0
    for ts,dt,old_name in zip(time_stamps,dts,base_names):
        delta = (dt-previous_dt).seconds
        if delta > 1:
            current_ts = ts
            counter = 0
        else:
            counter += 1
        name = '%s_%s_buff%02d.vdif'%(current_ts,antenna,counter)
        new_names.append(name)
        previous_dt = dt
        print 'renaming %s to %s'%(old_name,name)
