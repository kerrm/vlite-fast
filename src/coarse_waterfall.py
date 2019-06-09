""" Stripped down version of plotting to make a coarse waterfall plot.
"""

import matplotlib
matplotlib.use('Agg')
from sigpyproc.Readers import FilReader
import pylab as pl
import numpy as np
import argparse
import os

def fb_avg(fname,tgoal=1000,fgoal=128,fscr=None):
    """ Downsample a filterbank file to appropriate resolution for plotting.
    """
    fr = FilReader(fname)
    if fscr is None:
        fscr = int(fr.header['nchans']/fgoal)
    tscr = 2**int(np.log2(fr.header['nsamples']/tgoal))
    gulp = tscr

    nchan = fr.header['nchans']
    nblock = fr.header['nsamples']/tscr
    nsub = nchan/fscr
    output = np.empty( (nsub, nblock) )
    for nsamps,i,data in fr.readPlan(gulp,verbose=False):
        if nsamps != gulp:
            #print 'breaking at %d'%i
            break
        t = data.reshape((gulp,nchan)).transpose().astype(np.float32).mean(axis=1)

        output[:,i] = np.average(t.reshape(nchan/fscr,fscr),axis=1)

    freqs = np.arange(fr.header['nchans'])*fr.header['foff']+fr.header['fch1']
    freqs = freqs.reshape(len(freqs)/fscr,fscr).mean(axis=1)
    times = np.arange(nblock)*fr.header['tsamp']*gulp

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
    tmin = np.inf
    tmax = times[-1]+(times[1]-times[0])
    tmin = times[0]

    figx = 14
    figy = 4.0
    pl.figure(fignum,(figx,figy))
    pl.subplots_adjust(left=0.06,right=0.98,top=0.98,bottom=0.13)
    ax = pl.gca()

    ax.imshow(data,extent=[times[0],times[-1],freqs[0],freqs[-1]],vmin=0.9,vmax=1.1,cmap='Greys',aspect='auto',origin='lower')
    ax.set_ylabel('Freq (MHz)',size='large')
    ax.text(times[int(len(times)*0.95)],328,antenna,size='large',color='white')
    ax.set_xlabel('Elapsed Time (s)',size='large')
    ax.set_yticks([330,340,350])

    pl.savefig(plot_name)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Make a coarse waterfall plot of indicated file.")
    parser.add_argument('fname',help="Absolute path and name of file to plot.")
    parser.add_argument('--outpath',help="Save plots to this directory.  (Default is same directory as data.",default=None)

    args = parser.parse_args()

    plot_coarse_data(args.fname,outpath=args.outpath)
