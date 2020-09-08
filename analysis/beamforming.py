# Aug 19, 2019
# attempting to make a working example of beamforming with good, online
# source positioning, gains estimation, and a working example with B0329

# the rough order of operations of legacy code is
# (1) process the .vdif files into a single, aligned dataset with a common
#     length over all of the beams (baseband.re_order_baseband) --> .bb
# (2) use an OrderedBasebandIterator; delays can be input at this time to
#     give an iterator with (approximately) removed delays.
#     See baseband.solve_delays

# in fine channelization, use a channel mask, and cross correlation to
# estimate delays.  This will include clock delays *and* phase delays, but
# that's OK, since we can calculate the differential clock delays between
# 3C whatever and B0329.

# have just verified that the mask works well, replacing the silly
# channel range I had before...

# how do we actually do beamforming localization?  I think the same idea as
# with imaging: dedisperse, and then compare power on-pulse to power
# off-pulse.  Could form "typical" visibilities from a number of off-pulse
# chunks (but still phased up), then do classical beamforming to make the
# antenna patterns and subtract the off from the on.

# file wrangling: one has the thought that a better way to do this would
# be to take a set of baseband files, ingest them all at once, and for each
# file map out the sample range (e.g. by using a int64 to take
# second*frame to assign a unique sample, possibly using -ve to represent
# the other polarization, and with thise determine which files exactly
# need to be accessed to get which data

from __future__ import print_function,division

import numpy as np
import glob
import pyfftw

import baseband
from baseband import VDIFHeader



FRAMESPERSEC = 25600
FRAMELEN = 5032
FRAMESAMP = 5000
INVFRAMESAMP = 1./5000

class FileData(object):

    def __init__(self,fname,i0,i1):
        self.fname = fname
        self.i0 = i0
        self.i1 = i1
        self.thread = self.i0%2
        self.bufflen = (self.i1-self.i0+1)*FRAMESAMP
        self._samples = None
        # set up the first frames for each polarization individually
        if self.thread == 0:
            self._ti0s = [self.i0,self.i0+1]
        else:
            self._ti0s = [self.i0+1,self.i0]
        if (self.i1 & 1) == 1:
            self._ti1s = [self.i1-1,self.i1]
        else:
            self._ti1s = [self.i1,self.i1-1]

    def follows(self,other):
        return (self.i0-other.i1) == 1

    def _load_data(self):
        if self._samples is not None:
            return self._samples

        # we have at least some overlap, so load in data
        d = np.memmap(self.fname,dtype=np.uint8)
        # trim off header
        d = d.reshape((len(d)//FRAMELEN,FRAMELEN))[:,32:]
        # separate by polarization
        t0 = d[::2].reshape((d.shape[0]//2*d.shape[1]))
        t1 = d[1::2].reshape((d.shape[0]//2*d.shape[1]))
        if self.thread==1:
            t0,t1 = t1,t0
        self._samples = (t0,t1)
        return self._samples

    def contiguous(self):
        d = np.memmap(self.fname,dtype=np.uint32)
        h0 = VDIFHeader(d[:8])
        stride = h0.frame_length//4
        second = d[0::stride] & (2**30-1)
        frame = d[1::stride] & (2**24-1)
        threadid = (d[3::stride] & (2**26-2**16)) >> 16
        frame_idx = 2*(second.astype(np.int64)*FRAMESPERSEC+frame) + thread
        return np.all(np.diff(frame_idx) == 1)


    def fill_buffer(self,buff,frame_idx,frame_off,max_samp):
        """ Fill buffer with at most max_samp samples starting at the
            indicated frame_idx and frame_offs.

            Will determine the overlap of the requested data range with
            this file and return the number of samples filled.

            If the buffers are None for a particularly polarization, no
            data from that polarization is returned.

            buff == [slice pointing to thread 0 data, or None,
                     slice pointing to thread 1 data, or None]

            frame_idx == frame index for start of data for thread 0.  Will
                be incremented by 1 for thread 1.

            frame_off == sample offset from frame_idx (or frame_idx+1)

            max_samp  == maximum samples with which to fill buffer
        """

        if (frame_idx & 1) == 1:
            raise ValueError("frame_idx must be even! (thread=0)")

        frame_idx = np.asarray([frame_idx,frame_idx+1])

        # map the frame index onto internal buffers
        idx0 = (frame_idx-self._ti0s)//2*FRAMESAMP+frame_off
        idx1 = idx0 + max_samp
        no_overlap = np.logical_or(idx1 <= 0,idx0 >= self.bufflen//2)
        if np.all(no_overlap):
            return (0,0)

        idx0[idx0<0] = 0
        idx1[idx1>self.bufflen//2] = self.bufflen//2
        nsamp = idx1-idx0

        samples = self._load_data()

        for i in range(2):
            if buff[i] is None:
                nsamp[i] = 0
                continue
            #print i,nsamp[i],idx0[i],idx1[i],samples[i].shape,buff[i].shape
            buff[i][:nsamp[i]] = samples[i][idx0[i]:idx1[i]]

        return nsamp

class DataSet(dict):

    def __init__(self,datadict):
        self.data = datadict
        self.update(datadict)

    def get_fnames(self):
        fnames = []
        for ant in self['antennas']:
            for fd in self.data[ant]:
                fnames.append(fd.fname)
        return fnames

    def get_headers():
        pass

def load_dataset(fnames,check_contiguity=False,antennas=None): 
    """ Take a list of VDIF names and analyze the files.

    Parameters
    ----------
    fnames
        unsorted list of all VDIF files to consider
    check_contiguity
        exclude files with missing frames
    antennas
        an optional list of antennas to select; default is to include all
        that are implied by the fnames parameters

    Returns
    -------
    dataset
        a DataSet object, essentially a glorified dictionary
    """
        
    data = dict()
    headers = []
    for fname in fnames:
        
        d = np.memmap(fname,dtype=np.uint32)
        h0 = VDIFHeader(d[:8])
        if (antennas is not None) and (h0.antenna() not in antennas):
            continue
        h1 = VDIFHeader(d[len(d)-h0.frame_length//4:])
        stride = h0.frame_length//4

        # get start and stop indices
        i0 = h0.unique_idx()
        i1 = h1.unique_idx()

        # check that the data have length equal to the implied number of frames
        nframe = i1-i0+1
        data_ok = nframe*stride == len(d)
        if not data_ok:
            print('Found a problem with file %s.!  Discarding it.'%fname)
            continue
        
        # data for file:
        dat = FileData(fname, i0, i1)

        if check_contiguity:
            if not dat.contiguous():
                print('File %s was not contiguous!  Discarding it.'%fname)
                continue
        try:
            data[h0.antenna()]['data'].append(dat)
        except KeyError:
            data[h0.antenna()] = dict()
            data[h0.antenna()]['data'] = []
            data[h0.antenna()]['data'].append(dat)
            data[h0.antenna()]['header'] = h0

    antennas = list(data.keys())
    data['antennas'] = antennas
    for ant in antennas:
        segs = data[ant]['data']
        for s1,s2 in zip(segs[:-1],segs[1:]):
            # check thread 0
            if not s2.follows(s1):
                print('Warning!  Files %s and %s seem not to follow.'%(s1.fname,s2.fname))

    starting_frames = np.asarray([data[ant]['data'][0].i0 for ant in antennas])
    starting_frame = max(starting_frames)
    starting_frame += starting_frame%2
    stopping_frames = np.asarray([data[ant]['data'][-1].i1 for ant in antennas])
    stopping_frame = min(stopping_frames)
    stopping_frame -= (stopping_frame%2==0)

    # now, go through and find the earliest and latest complete frames 
    # shared by all antennas.  To ensure the same number of samples
    # between polarizations, always start on thread=0 and end on thread=1.
    data['starting_frame'] = starting_frame
    data['stopping_frame'] = stopping_frame

    return DataSet(data)

class NewBaseband(object):

    def __init__(self,dataset):
        self.antennas = dataset['antennas']
        self.nsamp = (dataset['stopping_frame']-dataset['starting_frame']+1)//2
        self.nsamp *= FRAMESAMP
        self._frame0 = dataset['starting_frame']
        self.dataset = dataset
        self.fsamp = 128000000

    def get_data(self,nsamp,off=0,thread=0,antennas=None,buff=None,
            antenna_offsets=None):
        """ Sample numbers are specified relative to the internal 0th
            frame, i.e. they start at 0.
        """

        if antennas is None:
            antennas = self.antennas
        antennas = sorted(set(antennas).intersection(antennas))

        if antenna_offsets is not None:
            aos = [key for key in antenna_offsets.keys() if key in antennas]
            aos = np.asarray([antenna_offsets[x] for x in aos])
        else:
            aos = np.zeros(len(antennas),dtype=int)

        if (off + max(aos) + nsamp) > len(self):
            raise ValueError("Requested data lies beyond end of data set.")

        if thread >= 0:
            buff_shape = (len(antennas),nsamp)
        else:
            buff_shape = (len(antennas),2,nsamp)
        if buff is None:
            buff = np.empty(buff_shape,dtype=np.float32)
        else:
            if not buff.shape == buff_shape:
                raise ValueError("Provided buffer has wrong shape!")

        for iant,ant in enumerate(antennas):
            aoff = off + aos[iant]
            b = buff[iant,off:]
            acc0 = 0
            acc1 = 0
            buff_t0 = buff_t1 = None
            if thread == -1:
                buff_t0 = buff[iant,0,acc0:]
                buff_t1 = buff[iant,1,acc1:]
            elif thread == 0:
                buff_t0 = buff[iant,acc0:]
            else:
                buff_t1 = buff[iant,acc1:]
            for fd in self.dataset[ant]['data']:
                n0,n1 = fd.fill_buffer([buff_t0,buff_t1],self._frame0,
                        aoff,nsamp)
                if n0 > 0:
                    buff_t0 = buff_t0[n0:]
                    acc0 += n0
                if n1 > 0:
                    buff_t1 = buff_t1[n1:]
                    acc1 += n1

        buff -= 128
        return buff

    def __len__(self):
        return self.nsamp

    def get_iterator(self,nsamp,global_offset=0,thread=0,overlap=0,
            antennas=None,antenna_offsets=None):

        return BasebandIterator(self,nsamp,global_offset=global_offset,
                thread=thread,overlap=overlap,antennas=antennas,
                antenna_offsets=antenna_offsets)


class BasebandIterator(object):
    """ Enumerate over baseband data from start to stop with given chunk
    and antenna size."""

    def __init__(self,bb,nsamp,global_offset=0,thread=0,overlap=0,
            antennas=None,antenna_offsets=None):
        self.bb = bb
        self.fsamp = bb.fsamp
        self.nsamp = nsamp
        self.off = global_offset
        self.overlap = overlap
        self.thread = thread
        self.antennas = antennas
        self.antenna_offsets = antenna_offsets
        if antenna_offsets is not None:
            eff_nsamp = len(bb)-global_offset-max(antenna_offsets.values())
        else:
            eff_nsamp = len(bb)-global_offset
        if nsamp > len(bb):
            raise ValueError('Requested step size longer than data.')
        if overlap > 0:
            self.nchunk = 1 + (eff_nsamp-nsamp)//(nsamp-overlap)
        else:
            self.nchunk =  eff_nsamp//nsamp
        # initialize buffer
        self.buff = bb.get_data(nsamp,off=global_offset,
                antenna_offsets=antenna_offsets,
                thread=self.thread,antennas=self.antennas)
        self._last_idx = None
        self.counter = 0

    def __len__(self):
        return self.nchunk

    def __iter__(self):
        self.counter = 0
        return self

    def next(self):
        self.counter += 1
        if self.counter == self.nchunk:
            raise StopIteration
        off = self.off + (self.counter-1)*(self.nsamp-self.overlap)
        return self.bb.get_data(self.nsamp,off=off,
                thread=self.thread,antennas=self.antennas,buff=self.buff,
                antenna_offsets=self.antenna_offsets)

    def __getitem__(self,idx):
        if idx == self._last_idx:
            return self.buff
        if idx < 0:
            idx = idx + self.nchunk
        if (idx >= self.nchunk) or (idx < 0):
            raise ValueError
        off = self.off + idx*(self.nsamp-self.overlap)
        return self.bb.get_data(self.nsamp,off=off,
                thread=self.thread,antennas=self.antennas,buff=self.buff,
                antenna_offsets=self.antenna_offsets)

    def get_shape(self):
        return self.buff.shape

    def samp_per_chunk(self):
        return self.nsamp

class FFTIterator(object):

    # TODO -- make this inherit baseband iterator

    def __init__(self,bbi,use_window=False):
        self.bbi = bbi
        self.nchunk = bbi.nchunk
        nbeam = bbi.buff.shape[0]
        nfft = bbi.buff.shape[1] 
        nchan = nfft//2+1
        # set up FFTs
        self.fbuff = pyfftw.empty_aligned((nbeam,nchan),dtype=np.complex64)
        #fbuff_conj = np.empty_like(fbuff)
        self.fftw = pyfftw.FFTW(bbi.buff[0],self.fbuff[0],
                direction='FFTW_FORWARD',flags=['FFTW_UNALIGNED'])
        if use_window:
            self.window = hamming(nfft,sym=False)
        else:
            self.window = None
        self._last_idx = None
        self.counter = 0

    def __len__(self):
        return self.nchunk

    def __iter__(self):
        self.counter = 0
        return self

    def _dofft(self,data):

        if self.window is not None:
            np.multiply(data,self.window[None,:],out=data)

        for ibeam in range(self.fbuff.shape[0]):
            self.fftw(data[ibeam],self.fbuff[ibeam])

        return self.fbuff

    def next(self):
        self.counter += 1
        if self.counter == self.nchunk:
            raise StopIteration
        return self._dofft(self.bbi[self.counter-1])

    def __getitem__(self,idx):
        if idx == self._last_idx:
            return self.fbuff
        if idx < 0:
            idx = idx + self.nchunk
        if (idx >= self.nchunk) or (idx < 0):
            raise ValueError
        return self._dofft(self.bbi[idx])

    def get_normalization(self):
        if self.window is None:
            return 1./(self.bbi.samp_per_chunk())
        else:
            return 1./np.sum(self.window**2)

    def get_shape(self):
        return self.fbuff.shape

def correlate(bbi,nchunk=None,normalize=True,use_window=False,alpha=0,
        use_time_window=False):
    """ Correlate the data in the baseband iterator.
    
    If alpha != 0, compute the cyclic correlation at frequency shift alpha.
    """

    ffti = FFTIterator(bbi,use_window=use_window)
    if nchunk is None:
        nchunk = len(bbi)
    nbeam = bbi.buff.shape[0]
    nfft = bbi.buff.shape[1] 
    nchan = nfft//2+1

    df = bbi.fsamp/nfft
    alpha_idx = int(round(alpha/df))
    if (alpha_idx*df != alpha):
        raise ValueError('Cyclic frequency not commensurate with FFT.')

    # set up FFTs
    fbuff_conj = np.empty_like(ffti.fbuff)

    # set up correlation matrix
    cmatrix = np.zeros((nbeam,nbeam,nchan),dtype=ffti.fbuff.dtype)
    tmp_cmatrix = np.empty_like(cmatrix)

    # evaluate window if being used
    if use_time_window:
        time_window = hamming(nchunk,sym=False)
    else:
        time_window = np.ones(nchunk)

    # calculate for each chunk
    for ichunk in range(nchunk):
        fbuff = ffti[ichunk]

        np.conjugate(fbuff,out=fbuff_conj)
        if alpha_idx != 0:
            inplace_roll(fbuff_conj,alpha_idx)

        np.multiply(fbuff[:,None,:],fbuff_conj,out=tmp_cmatrix)
        if use_time_window:
            tmp_cmatrix *= time_window[ichunk]
        cmatrix += tmp_cmatrix

    if normalize:
        np.multiply(cmatrix,ffti.get_normalization()/nchunk,out=cmatrix)
        if use_time_window:
            cmatrix *= 1./np.mean(time_window**2)

    return cmatrix

if __name__ == '__main__':
    # first, ingest a whole list of files gT
    #fnames = sorted(glob.glob('/data/kerrm/voltages/20190711/20190711_160008_*_156288276*.vdif'))
    fnames = sorted(glob.glob('/data/kerrm/voltages/crab/*vdif'))
    data = load_dataset(fnames)


    #fd = data['ea15'][0]
    nb = NewBaseband(data)
    bbi = nb.get_iterator(12500,thread=0)
    ffti = FFTIterator(bbi)

    tmp = ffti[0]
    results = np.empty([tmp.shape[0],tmp.shape[1],len(ffti)],dtype=np.complex64)
    import time
    t1 = time.time()
    for i in range(len(ffti)):
        results[...,i] = ffti[i]
    t2 = time.time()
    print(t2-t1)

#cmatrix = correlate(bbi,nchunk=100)

# thoughts on gain cal: probably need to do a truncated Fourier series to
# model the scalar gain.  Test this with the beamformed power?  Does it
# do better than the Nth degree polynomial?
