from __future__ import print_function,division
import numpy as np

def inplace_roll(a, shift):
    """ An in-place roll along last axis."""
    if shift == 0:
        return a
    mshift = a.shape[-1]-abs(shift)
    if shift < 0:
        shift = -shift
        shift,mshift = mshift,shift
    # copy out last elements
    tmp = a[...,mshift:].copy()
    # move other elements
    a[...,shift:] = a[...,:mshift]
    a[...,:shift] = tmp
    return a

def shift_multiply(a,b,shift,out=None,conjugate=False):
    """ Shift array b by shift slots and multiply.
    
    Equivalent to a*np.roll(b,shift) but with fewer copies/temps.  This is
    verified with timing, but need also to cook up something with 
    conjugation for it to be useful.
    """
    if out is None:
        out = np.empty_like(a)
    if shift == 0:
        if conjugate:
            np.conjugate(b,out=out)
            np.multiply(a,out,out=out)
        else:
            np.multiply(a,b,out)
        return out
    mshift = len(a)-abs(shift)
    if shift < 0:
        shift = -shift
        shift,mshift = mshift,shift
    if conjugate:
        np.conjugate(b[mshift:],out=out[:shift])
        np.conjugate(b[:mshift],out=out[shift:])
        np.multiply(a,out,out=out)
    else:
        np.multiply(a[:shift],b[mshift:],out=out[:shift])
        np.multiply(a[shift:],b[:mshift],out=out[shift:])
    return out

def roll_multiply(a,b,shift,out=None,conjugate=False):
    if conjugate:
        return np.multiply(a,np.roll(b.conj(),shift),out=out)
    else:
        return np.multiply(a,np.roll(b,shift),out=out)

def test_shift_multiply():

    a = np.arange(100000,dtype=np.complex64)
    b = (a - 1000).copy()
    b.imag[:] = 3
    c = np.empty_like(a)
    print(np.allclose(shift_multiply(a,b,2,c),roll_multiply(a,b,2,c)))
    print(np.allclose(shift_multiply(a,b,-2,c),roll_multiply(a,b,-2,c)))
    print(np.allclose(shift_multiply(a,b,-2,c,conjugate=True),
            roll_multiply(a,b,-2,c,conjugate=True)))
    print(np.allclose(shift_multiply(a,b,5,c,conjugate=True),
            roll_multiply(a,b,5,c,conjugate=True)))

def time_shift(a,delta):
    """
    Shift the real time series a by an amount delta (in samples)  using FFT.
    """
    kernel = np.exp(np.fft.rfftfreq(len(a))*(delta*2.0j*np.pi))
    return np.fft.irfft(kernel*np.fft.rfft(a))

def tophat_smooth(a,n,axis=0,out=None):
    """ Smooth array using a running mean of n samples.

    This can be relatively efficiently done using a series of cumulative
    sums.  Edge effects are dealt with by ramping up to the full convolution
    width, at the expense of additional noise (due to fewer samples.)

    n is the running mean width, and for symmetry and ease of index 
    manipulation will be forced to be odd.
    """
    ndim = len(a.shape)
    if(n%2==0):
        #print('Warning! n must be odd, so %d --> %d.'%(n,n+1))
        n += 1
    lena = a.shape[axis]
    if lena < n:
        return a
    c = np.cumsum(a,axis=axis)
    if out is None:
        out = np.empty_like(a)

    # multi-dimensional generalization
    if (ndim > 1):

        outslice = [slice(None),]*ndim
        outslice[axis] = slice(n//2+1,lena-n//2,None)
        cslice1 = [slice(None),]*ndim
        cslice1[axis] = slice(n,None,None)
        cslice2 = [slice(None),]*ndim
        cslice2[axis] = slice(None,lena-n,None)
        out[tuple(outslice)] = (c[tuple(cslice1)]-c[tuple(cslice2)])*(1./n)

        outslice[axis] = slice(None,n//2+1)
        cslice1[axis] = slice(None,n,2)
        normslice = [np.newaxis,]*ndim
        normslice[axis] = slice(None)
        out[tuple(outslice)] = c[tuple(cslice1)]/np.arange(1,n+1,2)[tuple(normslice)]

        outslice[axis] = slice(lena-n//2,None)
        cslice1[axis] = slice(lena-1,lena)
        cslice2[axis] = slice(lena-n+1,None,2)
        out[tuple(outslice)] = (c[tuple(cslice1)]-c[tuple(cslice2)])/(np.arange(1,n,2)[::-1])[tuple(normslice)]
        return out

    # single dimension slicing version
    out[n//2+1:len(a)-n//2] = (c[n:]-c[:len(a)-n])*(1./n)
    out[:n//2+1] = c[:n:2]/np.arange(1,n+1,2)
    out[len(a)-n//2:] = (c[-1]-c[len(a)-n+1::2])/np.arange(1,n,2)[::-1]
    return out

def unwrap_phase(phi):
    # TODO -- make this more efficient so it pre-calculates the phase at
    # each point and does a single addition, i.e. O(N) instead of O(N^2)
    dphi = phi[1:]-phi[:-1]
    change_points = np.ravel(np.argwhere(np.abs(dphi)>np.pi))
    for idx in change_points:
        if dphi[idx] > 0:
            phi[idx+1:] -= np.pi*2
        elif dphi[idx] < 0:
            phi[idx+1:] += np.pi*2
    return phi

def fave(spectrum,nbins,axis=-1,mask=None):
    """ Average nbins together."""
    if spectrum.shape[axis] % nbins != 0:
        raise ValueError('Spectrum not commensurate with averaging factor!')
    nout = spectrum.shape[axis] // nbins
    out_shape = np.copy(spectrum.shape)
    out_shape[axis] = spectrum.shape[axis] // nbins
    out = np.empty(out_shape,dtype=spectrum.dtype)
    if mask is None:
        for i in range(nout):
            sin = [slice(None),]*len(out_shape)
            sin[axis] = slice(i*nbins,(i+1)*nbins)
            sout = [slice(None),]*len(out_shape)
            sout[axis] = slice(i,i+1)
            sav = [slice(None),]*len(out_shape)
            sav[axis] = np.newaxis
            out[tuple(sout)] = np.average(spectrum[tuple(sin)],axis=axis)[tuple(sav)]
    else:
        for i in range(nout):
            sin = [slice(None),]*len(out_shape)
            sin[axis] = slice(i*nbins,(i+1)*nbins)
            sout = [slice(None),]*len(out_shape)
            sout[axis] = slice(i,i+1)
            sav = [slice(None),]*len(out_shape)
            sav[axis] = np.newaxis
            # assumes weights have shape of nchan
            try:
                out[tuple(sout)] = np.average(spectrum[tuple(sin)],axis=axis,weights=mask[sin[axis]])[tuple(sav)]
            except ZeroDivisionError:
                out[tuple(sout)] = np.nan
    return out

def dft(t,f):
    """ Brute force DFT with arbitrary frequencies.
    
    Frequencies are scaled to tsamp.  Scaling and sign convention are the
    same as numpy.fft.rfft.
    """
    f = f*(-2j*np.pi)
    rvals = np.empty(len(f),dtype=np.complex128)
    tmp1 = np.empty(len(t),dtype=np.complex128)
    tmp2 = np.empty_like(tmp1)
    tmp3 = np.empty_like(tmp1)
    dts = np.arange(len(t))
    for i in range(len(f)):
        np.multiply(f[i],dts,out=tmp1)
        np.exp(tmp1,out=tmp2)
        np.multiply(t,tmp2,out=tmp3)
        rvals[i] = np.sum(tmp3)
    return rvals

def qn(s):
    # calculate Qn as robust estimate of variance
    diffs = np.zeros((len(s)*(len(s)-1))//2)
    counter = 0
    for i in range(len(s)-1):
        newsamps = len(s)-i-1
        diffs[counter:counter+newsamps] = np.abs(s[i]-s[i+1:])
        counter += newsamps
    return 2.2219*np.percentile(diffs,25)
