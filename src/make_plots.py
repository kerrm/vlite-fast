""" Monitor files as they are written and launch plots.
"""

import time
import subprocess
import os

import matplotlib
matplotlib.use('Agg')
from sigpyproc.Readers import FilReader
import pylab as pl
import numpy as np

sleep_cycle = 20

data_dir = '/mnt/ssd/fildata'

def fb_avg(fname,tbingoal=1000,fgoal=128,fscr=None,last_minute=False):
    """ Downsample a filterbank file to appropriate resolution for plotting.
    """
    fr = FilReader(fname)
    if fscr is None:
        fscr = int(fr.header['nchans']/fgoal)
    nsamp = fr.header['nsamples']
    if (nsamp*fr.header['tsamp']) < 60:
        last_minute = False

    if last_minute:
        nsamp = int(60./fr.header['tsamp'])+1
    try:
        tscr = 2**int(np.log2(nsamp/tbingoal))
    except OverflowError:
        tscr = 1
    gulp = tscr

    nchan = fr.header['nchans']
    nblock = nsamp/tscr
    nsub = nchan/fscr
    output = np.empty( (nsub, nblock) )
    start = 0
    if last_minute:
        start = max(0,int(fr.header['nsamples']-nsamp)-1)
    for nsamps,i,data in fr.readPlan(gulp,start=start,verbose=False):
        if nsamps != gulp:
            #print 'breaking at %d'%i
            break
        t = data.reshape((gulp,nchan)).transpose().astype(np.float32).mean(axis=1)

        output[:,i] = np.average(t.reshape(nchan/fscr,fscr),axis=1)

    freqs = np.arange(fr.header['nchans'])*fr.header['foff']+fr.header['fch1']
    freqs = freqs.reshape(len(freqs)/fscr,fscr).mean(axis=1)
    times = (start + np.arange(nblock)*gulp)*fr.header['tsamp']

    return output,freqs,times

def plot_coarse_data(fname,outpath=None,fignum=2):
    """ Produce the waterfall plot for the beams at very coarse resolution.
    """

    # this is the "observation group"
    path,basename = os.path.split(fname)
    antenna = 'ea'+basename.split('ea')[1].split('_')[0][:2]
    outpath = outpath or path
    plot_name = os.path.join(outpath,
            basename.replace('.fil','.waterfall.png'))

    data,freqs,times = fb_avg(fname)

    quantiles = np.percentile(np.ravel(data),[2,5,10,50,90,95,98])
    print quantiles

    figx = 14
    figy = 4.0
    pl.figure(fignum,(figx,figy)); pl.clf()
    pl.subplots_adjust(left=0.06,right=0.98,top=0.98,bottom=0.13)
    ax = pl.gca()

    # NB formal min/max is 0.826,1.386 (?)
    ax.imshow(data,extent=[times[0],times[-1],freqs[0],freqs[-1]],vmin=quantiles[0],vmax=quantiles[-1],cmap='Greys',aspect='auto',origin='lower',interpolation='nearest')
    ax.set_ylabel('Freq (MHz)',size='large')
    ax.text(times[int(len(times)*0.95)],328,antenna,size='large',color='white')
    ax.set_xlabel('Elapsed Time (s)',size='large')
    ax.set_yticks([330,340,350])

    pl.savefig(plot_name)
    os.system('chmod 00644 %s'%plot_name)

    # do the same, but only plot most recent 60s
    data,freqs,times = fb_avg(fname,last_minute=True)
    quantiles = np.percentile(np.ravel(data),[2,5,10,50,90,95,98])
    print quantiles
    pl.clf()
    pl.subplots_adjust(left=0.06,right=0.98,top=0.98,bottom=0.13)
    ax = pl.gca()

    ax.imshow(data,extent=[times[0],times[-1],freqs[0],freqs[-1]],vmin=quantiles[0],vmax=quantiles[-1],cmap='Greys',aspect='auto',origin='lower',interpolation='nearest')
    ax.set_ylabel('Freq (MHz)',size='large')
    ax.text(times[int(len(times)*0.95)],328,antenna,size='large',color='white')
    ax.set_xlabel('Elapsed Time (s)',size='large')
    ax.set_yticks([330,340,350])

    plot_name = plot_name.replace('waterfall','waterfall_60s')
    pl.savefig(plot_name)
    os.system('chmod 00644 %s'%plot_name)

if __name__=='__main__':

    # continue to cycle and look for new observations

    # use a bit of a hack to do this.  To avoid any timestamping issues,
    # use a dummy file to establish a timestamp and then use this with find
    os.system('touch %s/timestamp_file'%data_dir)
    cmd = 'find %s -type f -cnewer %s/timestamp_file'%(data_dir,data_dir)

    while(True):

        time.sleep(sleep_cycle)

        # find all new files
        output = subprocess.check_output(cmd,shell=True)
        if len(output) == 0:
            continue
        print output

        files = [x for x in output.split('\n') if len(x) > 0]

        # process out observations -- NB this isn't quite right yet since
        # it treats everything as the same rather than distinguishing "kur"
        # observations.  Ultimately, it would be better to put kur up with
        # the data and muos so it distinguishes its own observation.
        observations = set([os.path.split(f)[-1].split('_ea')[0] for f in files])
        for observation in observations:

            # make sure there's a space for the outbound observation
            year = observation[:4]
            month = observation[4:6]
            outbound_dir = '/home/vlite-master/mtk/outbound/%s/%s/%s'%(year,month,observation)
            os.system('mkdir -p %s && chmod 00755 %s'%(outbound_dir,outbound_dir))

            print 'working on %s'%observation
            obs_files = [x for x in files if observation in x]

            fil_files = [f for f in obs_files if f.endswith('.fil')]
            for fil in fil_files:
                print 'Working on %s'%fil
                try:
                    plot_coarse_data(fil,outpath=outbound_dir)
                except Exception as e:
                    print e
                    print 'Failed on %s'%fil


            cand_files = sorted([f for f in obs_files if 'cand' in f])
            if len(cand_files) > 0:
                out_cand_file = os.path.split(cand_files[0])[-1][:-4]
                cand_files = ' '.join(cand_files)
                full_cand = '%s/%s'%(outbound_dir,out_cand_file)
                cat_cmd = 'cat %s > %s && chmod 00644 %s'%(cand_files,full_cand,full_cand)
                os.system(cat_cmd)

        os.system('touch %s/timestamp_file'%data_dir)

    os.system('rm %s/timestamp_file'%data_dir)
