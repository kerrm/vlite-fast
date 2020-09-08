""" Estimate total voltage sample offsets + sub-sample delays by analyzing
    the total cross power on each baseline, i.e. with intensity
    interferometry.
    
    I have tried to break the bands up into individual channels, but so
    far haven't made this robust (e.g. getting the right # of channels for
    a given S/N, etc.) so it remains TODO.
"""

import os
import pickle

import numpy as np
import pylab as pl
from scipy.signal import fftconvolve

import utils

def proc_buffs(p):
    rvals = dict()
    antennas = set()
    for key in p.keys():
        antennas.add(key.split('_')[0])
    antennas = sorted(list(antennas))
    rvals['antennas'] = antennas
    for antenna in antennas:
        rvals[antenna] = [p['%s_pol0'%antenna],p['%s_pol1'%antenna]]
    return rvals

#t0 = 1584408723.204
t0 = 1581646177.340
outdir = 'localization_output/%.3f'%(t0)
p0 = pickle.load(open(os.path.join(outdir,'step0.pickle'),'rb'))
p1 = pickle.load(open(os.path.join(outdir,'step1.pickle'),'rb'))
width = p0['width']
tsamp_fb = p0['tsamp']
buffs = proc_buffs(p1)

antennas = buffs['antennas']
delays = dict()

# we have saved 4 "widths", so the central two can be taken as the pulse,
# and the bounding 2 as offpulse

b0 = buffs[antennas[0]][0]
i0 = int(len(b0)*3./8)
i1 = int(len(b0)*5./8)
if (i1-i0)%2 == 1:
    i1 += 1
nchan = (i1-i0)//2+1
# number of delays per baseline
ndelay = 1
#if ndelay != 1:
#    raise ValueError('Stick with 1 delay across band for now.')

tmp = np.empty(nchan,dtype=np.complex128)
freqs = (np.fft.rfftfreq(i1-i0)*(2.0j*np.pi))
freqs_x0 = np.zeros(len(freqs))
edge_mask = int(round(len(freqs)*60./64)) # top edge of band
muos_mask = int(round(len(freqs)*24./64)) # bottom edge of band
delay_edges = np.round(len(freqs)*np.linspace(24,60,ndelay+1)*(1./64)).astype(int)


for iant in range(len(antennas)):
    for jant in range(iant+1,len(antennas)):
        for pol in [0,1]:
            all_x0 = np.empty(ndelay)
            all_err = np.empty(ndelay)
            ant1 = antennas[iant]
            ant2 = antennas[jant]
            baseline = '%s-%s-pol%d'%(ant1,ant2,pol)
            b1 = buffs[ant1][pol][i0:i1]
            b2 = buffs[ant2][pol][i0:i1]

            # do initial correlation to remove sample offset
            f1 = np.fft.rfft(b1)
            f2 = np.fft.rfft(b2)
            f1[edge_mask:] = 0
            f1[:muos_mask] = 0
            f2[edge_mask:] = 0
            f2[:muos_mask] = 0
            conv = np.fft.irfft(f1*f2.conj())
            idx0 = np.argmax(np.abs(conv))
            if idx0 > len(b1)//2:
                idx0 = idx0-len(b1)
            b2 = np.roll(buffs[ant2][pol],idx0)[i0:i1]
            f2 = np.fft.rfft(b2)
            f2[edge_mask:] = 0
            f2[:muos_mask] = 0

            # now compute visibility/cross-power over a range of delays
            viz = f1*f2.conj()
            # TODO -- restrict this range
            dom = np.linspace(-5,5,101)
            cod = np.empty((len(dom),ndelay))
            for ix,x in enumerate(dom):
                np.exp(freqs*x,out=tmp)
                tmp *= viz
                for idelay in range(ndelay):
                    chan_vals = tmp[delay_edges[idelay]:delay_edges[idelay+1]].real
                    cod[ix,idelay] = np.abs(np.sum(chan_vals))

            pl.clf()
            for idelay in range(ndelay):
                x = cod[:,idelay]
                pl.plot(dom,(x-x.min())/(x.max()-x.min()))

                a = np.argmax(x)
                # handle bad (low S/N) baselines where there's no local max
                if a == 0:
                    a = 1
                if a == len(x)-1:
                    a = len(x)-2
                p = np.polyfit(dom[a-1:a+2],x[a-1:a+2],2)
                x0 = -p[1]/p[0]*0.5
                fmax = -p[1]**2/(4*p[0])+p[2]
                err = (-0.5*fmax/p[0])**0.5
                freqs_x0[delay_edges[idelay]:delay_edges[idelay+1]] = x0
                all_x0[idelay] = x0
                all_err[idelay] = err

            pl.savefig(os.path.join(outdir,'new_delay1_%s.png'%baseline))
            
            # rotate visibility to be all real, to the extent possible
            bestviz = viz*np.exp(freqs*freqs_x0)
            total_angle = np.arctan2(bestviz.imag.sum(),bestviz.real.sum())
            bestviz *= np.exp(-1j*total_angle)
            nsmooth = int(nchan/64)
            print('nsmooth=',nsmooth)
            xi = utils.tophat_smooth(bestviz.imag,nsmooth)
            xr = utils.tophat_smooth(bestviz.real,nsmooth)
            xt = utils.tophat_smooth(np.abs(bestviz),nsmooth)
            phase = np.arctan2(xi,xr)

            # calculate final signal to noise ratio from IFFT
            acorr = np.fft.irfft(bestviz)
            acorr_std = np.std(acorr[int(len(acorr)*0.05):int(len(acorr)*0.95)])
            if (acorr.max() != acorr[0]):
                print('Autocorrelation max not at delay max!')
            acorr_sn = acorr[0]/acorr_std

            # estimate coherence from real parts vs. total intensity
            # NB this measurement is kindof tough, best to used smoothed
            # parts, I think
            #coherence = np.mean(bestviz.real)/np.mean(np.abs(bestviz))
            #coherence = np.mean(np.abs(bestviz.imag))/np.mean(np.abs(bestviz))
            #mask = np.abs(bestviz) > 0
            #coherence = np.mean(bestviz.real[mask]/np.abs(bestviz[mask]))
            mask = xt > 0
            #coherence = np.mean(xr[mask]/xt[mask])
            coherence = 1-np.mean(np.abs(xi[mask])/xr[mask])
            coherence = 1-np.mean(np.abs(xi[mask]))/np.mean(xr[mask])
            coherence2 = np.mean(xr[mask])/np.mean(xt[mask])
            print('Accor S/N = %.2f, Coherence = %.2f'%(acorr_sn,coherence))

            pl.figure(2); pl.clf()
            pl.plot(xt[::nsmooth//2]/xt.max(),color='k',label='Total Visibility')
            pl.plot(xr[::nsmooth//2]/xt.max(),color='C0',label='Real',ls='-')
            pl.plot(xi[::nsmooth//2]/xt.max(),color='C1',label='Imaginary',ls='--')
            pl.plot(phase[::nsmooth//2],color='red',label='Phase = %.2f'%total_angle,alpha=0.5)
            pl.axis([0,len(bestviz)//(nsmooth//2),-0.8,1.2])
            pl.xlabel('Channel')
            pl.figtext(0.15,0.82,'S/N = %.2f'%acorr_sn,size='large')
            pl.figtext(0.15,0.78,'Coherence = %.2f / %.2f'%(coherence,coherence2),size='large')
            pl.legend(loc='lower left')
            pl.title('Baseline %s'%(baseline))
            pl.savefig(os.path.join(outdir,'step2_qa_%s.png'%baseline))
            # save an extra version of the plot to make it easy to identify
            # bad antennas
            toks = baseline.split('-')
            switched_baseline = '%s-%s-pol%s'%(toks[1],toks[0],pol)
            pl.title('Baseline %s'%switched_baseline)
            pl.figtext(0.85,0.92,'MIRROR',size='large')
            pl.savefig(os.path.join(outdir,'new_delay2_%s.png'%switched_baseline))

            # Yes, the total power not agreeing with "real" is real!  There
            # seem to be lots of small fluctuations (possibly real structure
            # in the pulse?) that are smoothed over in the plot, but apparently
            # that don't follow the global structure of the bandpass.

            delays[baseline] = (idx0,all_x0[0],all_err[0],total_angle,
                    acorr_sn,coherence,coherence2,bestviz)

pickle.dump(delays,open(os.path.join(outdir,'step2_test.pickle'),'wb'))
