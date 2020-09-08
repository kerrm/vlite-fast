""" Do preliminary processing:
        Load and analyze voltage dataset
        Produce filterbanked total power time series
        Using metadata from candidate file (DM, pulse time, pulse width),
            re-optimize pulse location and DM
        Output information in a pickle for Step 1
"""
from __future__ import print_function,division
import glob
import pickle
import time
import os

import numpy as np
import pylab as pl

import baseband
from baseband import VDIFHeader
import beamforming
import utils
from utils import inplace_roll

"""
From Surya, RE timestamping:
RE i0 and i1
I don't think I get what you asked. i0 and i1 is set by the trigger dispatch. Current convention is:
    i0 = t - 0.2
    i1 = i0 + dm_delay + 0.15 + (30-dm-delay)
    t is the arrival time in the reference frequency (highest frequency). 0.2 determines the peak-time in following fbson/dbson.
    dm_delay is the trigger dm delay. 0.15 is the width consideration (our max width is 0.1, extra 0.05 added for safety).
    The 30-dm-delay is dm_delay for a DM=30pc/cc added for bowtie plane consideration (we only go over trigger dm by 25pc/cc, 5pc/cc for safety)

That is for TRIGGER.  But note that for "dbson" times (i1-i0=0.2) the peak
should be at ~0.1.
"""


# dedisperse our data.  From FFTI, the highest channel is 320 MHz.

def get_vlite_chan_freqs(nchan):
    """ Return VLITE channel frequencies in MHZ."""
    return ((np.arange(nchan))*64./nchan + 320)[::-1]

def dedisperse(array,dm,tsamp,ref_freq=320):
    """ Incoherently dedisperse filterbank data.

    Parameters
    ----------
    array
        Filterbank data array: antenna x channel x time
    dm
        Dispersion measure to which to dedisperse
    tsamp
        Filterbank sample time (s)
    ref_freq
        Frequency to at which dispersion delay is 0 (MHz)
    
    NB there is some discrepancy between the reference frequency conventions
    we use.  Surya typically dedisperses to 320 MHz, whereas I have been
    doing ~361 MHz, the effective top of the band before MUOS RFI kicks in.
    """
    nchan = array.shape[1]
    chan_freqs = get_vlite_chan_freqs(nchan)
    sample_delays = -np.round(dm*4.15e-3*((chan_freqs*1e-3)**-2-(ref_freq*1e-3)**-2)/tsamp).astype(int)
    for ichan in range(nchan):
        inplace_roll(array[:,ichan],sample_delays[ichan])

def build_filterbanks(data,nsamp=12500,navg=10,antennas=None,pol=-1):
    """ Quick and dirty method for filterbanking and square-law detecting
    the voltage data.
    """

    nb = beamforming.NewBaseband(data)
    if pol == -1:
        bbi_p0 = nb.get_iterator(nsamp,thread=0,antennas=antennas)
        bbi_p1 = nb.get_iterator(nsamp,thread=1,antennas=antennas)
        ffti_p0 = beamforming.FFTIterator(bbi_p0)
        ffti_p1 = beamforming.FFTIterator(bbi_p1)
    else:
        bbi_p0 = nb.get_iterator(nsamp,thread=pol,antennas=antennas)
        ffti_p0 = beamforming.FFTIterator(bbi_p0)
        bbi_p1 = None
        ffti_p1 = None

    nout = len(ffti_p0)//navg
    tmp = ffti_p0[0]
    results = np.zeros([tmp.shape[0],tmp.shape[1],nout],dtype=np.float32)
    prof_t1 = time.time()
    counter = 0
    tmp = np.empty(tmp.shape,dtype=np.float32)
    for i in range(nout):
        for j in range(navg):
            px = ffti_p0[counter]
            np.abs(px,out=tmp)
            tmp *= tmp
            results[...,i] += tmp
            if pol >= 0:
                counter += 1
                continue
            py = ffti_p1[counter]
            np.abs(py,out=tmp)
            tmp *= tmp
            results[...,i] += tmp
            counter += 1

    prof_t2 = time.time()
    tsamp = (nsamp*navg)/128e6
    print('Filterbanking took %.2fs.'%(prof_t2-prof_t1))
    return results,tsamp

def chan_mask():
    # explicitly zero out some RFI for the nchan=6251 case.
    mask = np.zeros(6251)
    mask[2350:6200] = 1
    mask[3123:3128] = 0
    mask[4297] = 0
    mask[4988:4992] = 0
    return mask

def optimize_pulse(ts,i0,i1,wmax=32):
    """ Determine optimal extraction parameters in terms of tophat width
    and pulse location.

    Parameters
    ----------
    ts
        frequency scrunched time series
    i0
        sample no. at start of dbson-like interval (~0.1s before pulse rise)
    i1
        sample no. at end of dbson-like interval (~0.1s after pulse rise)
    """
    widths = np.arange(1,wmax+1,2)
    nsamp = i1-i0
    # take the first and last 25% of the samples
    s = np.append(ts[i0:i0+int(nsamp*0.25)],ts[i0+int(nsamp*0.75):i1])
    # calculate Qn as robust estimate of variance
    mean = np.median(s)
    std = utils.qn(s)
    sns = np.zeros(len(widths))
    locs = np.zeros(len(widths))
    for iw,w in enumerate(widths):
        sts = utils.tophat_smooth(ts[i0:i1],w)[w:-w]
        a = np.argmax(sts) 
        locs[iw] = a+w
        sns[iw] = (sts[a]-mean)/std*w**0.5
    return widths,sns,locs

def refine_dm(filterbanks,tsamp,dm0,i0,i1,delta_dm=0.001,dm_range=0.1,
        width=None):
    """ Quick and dirty brute force search through DM space to refine the
    DM determined by the on-line candidate search.
    """
    # go ahead and average antennas
    F = filterbanks.mean(axis=0)[None,...]
    nstep = int(round(2*dm_range/delta_dm))
    dms = dm0 + np.linspace(-dm_range,dm_range,nstep)
    mask = chan_mask()
    sns = np.empty_like(dms)
    for idm,dm in enumerate(dms):
        f = F.copy()
        dedisperse(f,dm,tsamp,ref_freq=fref)
        q = (f[0]*mask[:,None]).sum(axis=0)
        mean = np.median(q[i0:i1])
        std = utils.qn(q[i0:i1])
        if width is not None:
            q = utils.tophat_smooth(q,width)
            std *= 1./width**0.5
        sns[idm] = (q[i0:i1].max()-mean)/std
    return dms,sns


# comment in dedisperse method -- this reference frequency for dedispersion
# is at the top of the band
fref = 361.9414488190918

#fnames = sorted(glob.glob('/data/kerrm/voltages/B0329+54/20200316*.vdif'))
fnames = sorted(glob.glob('/data/kerrm/voltages/crab/20200213*vdif'))
# these values copied from Surya's candidate metadata, obviously can
# automate this...
#t0, dm0, width = 1584408723.203906,26.613115310668945,0.012500000186264515
t0, dm0, width = 1581646177.3398435,56.771188, 0.00625
# NB from Surya -- requested dump time is 1581646177.24

# make output directory -- these are catalogued by their UNIX timestamp
# to millisecond precision.  Currently, I'm just keying in this timestamp
# manually in subsequent steps, but obviously if these are chained as a
# pipeline that would be passed along.
outdir = 'localization_output/%.3f'%(t0)
os.system('mkdir %s'%outdir)

# process VDIF files to determine metadata, alignment, etc.
data = beamforming.load_dataset(fnames)
antennas = data['antennas']

# preliminary processing
# (1) build filterbanks
# (2) using "dbsons" style, find optimal width + location for pulse extraction
# (3) using optimal width, find optimal DM
# (4) proceed with extraction and use location + position for optimal extraction


filterbanks,tsamp = build_filterbanks(data,navg=8,pol=-1)
process_filterbanks = filterbanks.mean(axis=0)[None,...] # Antenna averaged
orig_filterbanks = process_filterbanks.copy()
mask = chan_mask()
dedisperse(process_filterbanks,dm0,tsamp,ref_freq=fref)
q1 = (process_filterbanks[0]*mask[:,None]).sum(axis=0)

i0 = int(round((t0-data[data['antennas'][0]]['header'].get_unix_timestamp())/tsamp))
i1 = int(round((t0-data[data['antennas'][0]]['header'].get_unix_timestamp()+0.2)/tsamp))

# iterate with DM refinement
widths,sns,locs = optimize_pulse(q1,i0,i1)
width = widths[np.argmax(sns)]
dms,sns = refine_dm(orig_filterbanks,tsamp,dm0,i0,i1,delta_dm=0.02,dm_range=0.2,width=width)
pl.figure(1); pl.clf()
pl.plot(dms,sns)
dm = dms[np.argmax(sns)]
dms,sns = refine_dm(orig_filterbanks,tsamp,dm,i0,i1,delta_dm=0.001,dm_range=0.02,width=width)
pl.plot(dms,sns)
pl.xlabel('DM')
pl.ylabel('S/N')
pl.savefig(os.path.join(outdir,'step0_dmopt.png'))
dm = dms[np.argmax(sns)]
process_filterbanks[:] = orig_filterbanks
dedisperse(process_filterbanks,dm,tsamp,ref_freq=fref)
q1 = (process_filterbanks[0]*mask[:,None]).sum(axis=0)
widths,sns,locs = optimize_pulse(q1,i0,i1)
pl.clf()
pl.plot(widths,sns)
pl.xlabel('Width (Samples)')
pl.ylabel('S/N')
pl.savefig(os.path.join(outdir,'step0_widthopt.png'))
width = widths[np.argmax(sns)]
loc = int(locs[np.argmax(sns)])
peak_loc = loc + i0

# process filterbanks with final DM and width
dedisperse(filterbanks,dm,tsamp,ref_freq=fref)
q = (filterbanks*mask[None,:,None]).mean(axis=1)
qs = utils.tophat_smooth(q,width,axis=1)


# make some plots and export a package for Step 1
pl.clf()
x = q1[i0:i1]
tmp = np.append(x[:loc-2*width],x[loc+2*width:])
std = utils.qn(tmp)
mean = np.median(tmp)
x = (x-mean)/std
pl.plot(x)
x = utils.tophat_smooth(q1,width)[i0:i1]
x_coadd = (x-mean)/std*width**0.5
pl.plot(x_coadd)
pl.axvline(loc-width,color='C2')
pl.axvline(loc+width,color='C2')
pl.savefig(os.path.join(outdir,'step0_profile.png'))
pl.clf()
for iantenna,antenna in enumerate(antennas):
    x = q[iantenna][i0:i1]
    tmp = np.append(x[:loc-2*width],x[loc+2*width:])
    std = utils.qn(tmp)
    mean = np.median(tmp)
    x = (qs[iantenna][i0:i1]-mean)/std*width**0.5
    resid = x-x_coadd/np.sqrt(len(antennas))
    bad = resid[loc] < 0
    label = '%s %s'%(antenna,'b' if bad else ' ')
    pl.plot(resid,label=label)
pl.axis([loc-5*width,loc+5*width,pl.axis()[2],pl.axis()[3]])
pl.legend(loc='lower left',ncol=3)
pl.savefig(os.path.join(outdir,'step0_profile2.png'))


# peak location information
meta = [t0,tsamp,peak_loc,width,dm,fref]
data = dict()
data['t0'] = t0
data['tsamp'] = tsamp
data['peak'] = peak_loc
data['width'] = width
data['dm'] = dm
data['fref'] = fref
data['fnames'] = fnames
data['antennas'] = antennas
pickle.dump(data,open(os.path.join(outdir,'step0.pickle'),'wb'))


# some old analyses preserved for my benefit
#fnames = sorted(glob.glob('/data/kerrm/voltages/B1749-28/20200304_20075[5,6]_*_158340633[4,5].vdif'))
#dm0 = 50.27
#t0=1583406334.19
#t1=1583406335.34
#width = 0.00078125*12 # replace with real data

#fnames = sorted(glob.glob('/data/kerrm/voltages/crab/20200213*.vdif'))
#dm0=56.75
#dm0 = 56.7357
#t0=1581646177.24
#width = 0.00078125*8 # replace with real data

#fnames = sorted(glob.glob('/data/kerrm/voltages/B0329+54/20200312_10131*_*5[3].vdif'))
#1584034751.0203123  26.613115310668945  0.010937499813735485
#t0, dm0, width = 1584034753.1609373,26.88109588623047,0.010937499813735485
#1584034753.8749998  26.88109588623047   0.010937499813735485
#1584042943.7757812  177.4746856689453   0.09843750298023224
#1584034751.7265623  27.149089813232422  0.01406249962747097
