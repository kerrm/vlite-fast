""" Using information from Step 0, do coherent dedispersion and trim
    the coherently dedispersed voltages down to just the pulse and a bit
    of the surrounds for context.
"""

from __future__ import print_function,division
import glob
import os
import pickle
import time

import numpy as np
import pylab as pl
import pyfftw
from scipy.interpolate import interp1d  
from scipy.signal import fftconvolve

import baseband
from baseband import VDIFHeader
import beamforming
import utils

# TODO -- take T0 as input argument
#t0 = 1584408723.204
t0 = 1581646177.340

# Load information about dataset and peak positions
outdir = 'localization_output/%.3f'%(t0)
p = pickle.load(open(os.path.join(outdir,'step0.pickle'),'rb'))
t0 = p['t0']
dm = p['dm']
peak_loc = p['peak']
width = p['width']
fref = p['fref']
tsamp_fb = p['tsamp']
antennas = p['antennas']

data = beamforming.load_dataset(p['fnames'],antennas=antennas)

freq_hi = 384.
freq_lo = 320.
t_dm = dm/2.41e-4*(1./(freq_lo*freq_lo)-1./(freq_hi*freq_hi))
tsamp = 1./128e6
n_dm_samp = int(round(t_dm/tsamp))

# don't confuse "overlap" with "advance".  We need "advance" to be
# < n_dm_samp, which means that "overlap" may be large (approaching
# the length of the FFT)

nsamp = int(128e6)
nice_factors = np.asarray([20,10,8,5,4,2])
m = (nsamp-nsamp/nice_factors)/n_dm_samp > 1
factor = (nice_factors[m]).min()
overlap = nsamp-nsamp//factor

buffs = dict()
nb = beamforming.NewBaseband(data)
outbuff = np.empty(nb.nsamp,dtype=np.float32)
outbuff[:] = np.nan
fbuff = pyfftw.empty_aligned(nsamp//2+1,dtype=np.complex64)

nchan = nsamp//2 + 1
freqs = (np.arange(nchan)*(64./nchan))[::-1] # upper sideband
freq0 = 320.

# these phases will dedisperse to 320 MHz, which makes the overlap in the
# dedispersion much easier to handle, and later on we will correct
# the time delay for the reference frequency used in Step 1
arg = (2*np.pi*dm/2.41e-10)*freqs*freqs/(freq0*freq0*(freq0+freqs))
kernel = np.empty(len(arg),dtype=np.complex64)
kernel.real = np.cos(arg)
kernel.imag = np.sin(arg)
del(arg)

for antenna in antennas:
    for pol in [0,1]:

        print('Working on antenna=%s, pol=%d'%(antenna,pol))

        # choose a number of samples that will fit nicely in with the data,
        # 12500 is a good factor; we also don't want to have too little
        # efficiency.  For DM ~50, 1s is convenient. 
        bbi_p0 = nb.get_iterator(nsamp,thread=pol,overlap=overlap,
                antennas=[antenna])

        buff = bbi_p0.buff[0]
        assert(nchan==buff.shape[0]//2+1)
        # for now, let's do this unencapsulated
        try:
            wisdom = pickle.load(open('pyfftw.wisdom','rb'))
            pyfftw.import_wisdom(wisdom)
        except Exception as e:
            print('Could not load wisdom because:')
            print(e)
            pass
        fftw = pyfftw.FFTW(buff,fbuff,direction='FFTW_FORWARD',flags=[
            'FFTW_MEASURE','FFTW_DESTROY_INPUT'],planning_timelimit=60)
        fftwi = pyfftw.FFTW(fbuff,buff,direction='FFTW_BACKWARD',flags=[
            'FFTW_MEASURE','FFTW_DESTROY_INPUT'],planning_timelimit=60)
        wisdom = pyfftw.export_wisdom()
        pickle.dump(wisdom,open('pyfftw.wisdom','wb'))

        # forward transform
        fftw(buff,fbuff)

        # use the first sample to remove narrowband RFI and shape bandpass
        # NB do I really want to flatten the bandpass?

        # first, remove narrowband RFI and shape bandpass
        pa = utils.fave(np.abs(fbuff)[1:],100) # remember to mask out 0th channel going forward
        q = np.abs(pa[1:]-pa[:-1])
        narrow_mask = q > np.median(q)*10
        narrow_mask = np.append(True,narrow_mask) | np.append(narrow_mask,True)
        x = np.arange(len(pa))
        # replace the narrowband channels with the interpolated values
        pa[narrow_mask] = np.interp(x[narrow_mask],x[~narrow_mask],pa[~narrow_mask])
        # make a coarse version of the bandpass to estimate interpolants
        pa2 = utils.fave(pa,100)
        fa2 = utils.fave(freqs[1:],10000)
        ip = interp1d(fa2,pa2,kind='linear',bounds_error=False,fill_value=(pa2[-1],pa2[0]))
        bp = ip(freqs)
        bp_mask = np.zeros(len(bp),dtype=bool)
        m = np.ravel(np.argwhere(narrow_mask))
        for idx in m:
            bp_mask[idx*100:(idx+1)*100] = True

        # with this bandpass scaling, real/imag both are unit normal in Four dom.
        bp *= 1./2**0.5

        # now begin processing
        for ichunk in range(bbi_p0.nchunk):
            #print('processing chunk %d/%d'%(ichunk,bbi_p0.nchunk))

            if ichunk > 0:
                buff = bbi_p0[ichunk][0]
                fftw(buff,fbuff)

            # now, backtrack this to modify / scale the kernel
            fbuff.real[bp_mask] = np.random.randn(bp_mask.sum())
            fbuff.imag[bp_mask] = np.random.randn(bp_mask.sum())

            fbuff *= kernel
            fbuff /= bp
            # zero out MUOS
            fbuff[:int(nchan*24./64)] = 0

            fftwi(fbuff,buff)

            if ichunk == 0:
                outbuff[:nsamp] = buff
            else:
                # good data (with freq=320) lies in buff[n_dm_samp:]
                i0 = ichunk*(nsamp-overlap)
                outbuff[i0+n_dm_samp:i0+nsamp] = buff[n_dm_samp:]


        # this point, we have processed everything and just need to get the
        # pulse
        tstart = tsamp_fb*(peak_loc-2*width) + dm*4.15e-3*(0.320**-2-(fref*1e-3)**-2)
        tstop = tsamp_fb*(peak_loc+2*width) + dm*4.15e-3*(0.320**-2-(fref*1e-3)**-2)
        samp0 = int(round(tstart/tsamp))
        samp1 = int(round(tstop/tsamp))

        buffs['%s_pol%d'%(antenna,pol)] = outbuff[samp0:samp1].copy()
        # TODO -- store additional metadata on pulse position in the buffer
        # to abstract away DM etc. downstream

pickle.dump(buffs,open(os.path.join(outdir,'step1.pickle'),'wb'))
