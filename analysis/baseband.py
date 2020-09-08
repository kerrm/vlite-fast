from __future__ import print_function,division
from calendar import timegm
from collections import deque
import glob
import time

import numpy as np  
import pyfftw
import pylab as pl
from scipy.signal import fftconvolve,hamming,hanning
from scipy.optimize import leastsq

from utils import inplace_roll

FRAMESPERSEC = 25600

class VDIFHeader(object):

    def __init__(self,hdr):
        d = hdr.copy()
        self.second = d[0] & (2**30-1)
        self.epoch = (d[1] & (2**30-2**24)) >> 24
        self.frame = d[1] & (2**24-1)
        self.frame_length = frame_length = (d[2] & (2**24-1))*8
        self.frame_nsamp = frame_length-32
        self.station = d[3] & (2**16-1)
        self.threadid = (d[3] & (2**26-2**16)) >> 16
        self.thread = int(self.threadid != 0)

    def __str__(self):
        s = ['Reference epoch = %d'%self.epoch,
            'Second = %d'%self.second,
            'Threadid = %d'%self.threadid,
            'Thread = %d'%self.thread,
            'Frame = %d'%self.frame,
            'Frame length (samples) = %d'%self.frame_nsamp,
            'Station = %d'%(self.station),
            'UTC = %s'%(self.get_utc_str()),]
            #'File length (samples) = %d'%(self.bufflen)]
        return '\n'.join(s)

    def unique_idx(self):
        return 2*(self.second*FRAMESPERSEC+self.frame) + self.thread

    def antenna(self):
        return 'ea%02d'%(self.station)

    def __cmp__(self,other):
        if other.second < self.second:
            return 1
        if other.second > self.second:
            return -1
        if other.frame < self.frame:
            return 1
        if other.frame > self.frame:
            return -1
        if other.thread < self.thread:
            return 1
        if other.thread > self.thread:
            return -1
        return 0

    def follows(self,other):

        if other.threadid == self.threadid:
            return False

        if other.second == self.second:
            if self.thread==1:
                return self.frame-other.frame == 0
            else:
                return self.frame-other.frame == 1

        # can only follow if we are wrapping to a new second
        if ((other.second == self.second-1) and
            (other.frame == (FRAMESPERSEC -1)) and 
            (self.frame == 0)):
            return True
        return False

    def get_utc_str(self,astropy=False):
        """ Easiest to simply use UNIX epoch and then take advantage of
            time system.  Don't worry about frame offset.

            We know the UNIX epoch in MJD, and until 2032 the VDIF epoch
            will be Jan 1, 2000, so simply add the required number of
            seconds to the VDIF seconds field to convert this to a UNIX
            time_t, then get the GMT string.
        """
        # first compute offset using a time struct
        yy = 2000 + self.epoch//2
        mm = (self.epoch%2)*6 + 1
        dd = 1
        epoch_sec = timegm([yy,mm,dd,0,0,0,0,0,0])
        dt_sec_unix = epoch_sec + self.second
        gmt = time.gmtime(dt_sec_unix)
        if astropy:
            return time.strftime('%Y-%m-%d-%H:%M:%S',gmt)
        else:
            return time.strftime('%Y-%m-%d %H:%M:%S',gmt)

    def get_unix_timestamp(self):
        """ Easiest to simply use UNIX epoch and then take advantage of
            time system.  Don't worry about frame offset.

            We know the UNIX epoch in MJD, and until 2032 the VDIF epoch
            will be Jan 1, 2000, so simply add the required number of
            seconds to the VDIF seconds field to convert this to a UNIX
            time_t, then get the GMT string.
        """
        # first compute offset using a time struct
        yy = 2000 + self.epoch//2
        mm = (self.epoch%2)*6 + 1
        dd = 1
        epoch_sec = timegm([yy,mm,dd,0,0,0,0,0,0])
        return epoch_sec + self.second


class BasebandFragment(object):
    """ Encapsulate a single baseband file."""

    def __init__(self,fname):
        self.fname = fname
        self._load_header()
        # this is the *internal* 
        #self.thread = self.header.thread

    def __len__(self):
        """ Return length in samples."""
        return self.bufflen

    def get_nframe(self,thread=None):
        if thread is None:
            return self.nframe
        else:
            return self.thread_nframe[thread]

    def get_nsamp(self,thread=0):
        return self.thread_nframe[thread]*self.frame_nsamp

    #def find_frame(self,frame):
        #""" Return the sample offset
    
    def check_contiguity(self):
        hdrs = self.get_headers()
        for h2,h1 in zip(hdrs[1:],hdrs[:-1]):
            if not h2.follows(h1):
                return False
        return True

    def _load_header(self):
        d = np.memmap(self.fname,dtype=np.uint32)
        h = self.header = VDIFHeader(d)
        fl = h.frame_length
        self.__dict__.update(h.__dict__)
        if len(d) % (fl//4) != 0:
            raise ValueError('Non-standard dump length!')
        self.nframe = len(d)/(fl//4)
        self.thread_nframe = [self.nframe//2]*2
        if self.nframe % 2 != 0:
            self.thread_nframe[self.header.thread()] += 1
        self.bufflen = (len(d)*4)//fl*self.frame_nsamp

        h0 = h
        h1 = VDIFHeader(d[fl:fl+8])
        i0 = len(d)-fl//4
        self.last_header = VDIFHeader(d[i0:i0+8])

    def get_headers(self,max_headers=None):
        d = np.memmap(self.fname,dtype=np.uint32)
        nframe = self.nframe
        if max_headers is None:
            max_headers = nframe
        d = np.reshape(d,(nframe,self.frame_length//4))[:max_headers,:4]
        hdrs = map(VDIFHeader,d)
        return hdrs

    def antenna(self):
        return 'ea%02d'%self.station

    def __str__(self):
        s = [self.fname,
            'Reference epoch = %d'%self.epoch,
            'Second = %d'%self.second,
            'Threadid = %d'%self.threadid,
            'Frame = %d'%self.frame,
            'Frame length (samples) = %d'%self.frame_length,
            'Station = %d'%(self.station),
            'File length (samples) = %d'%(self.bufflen)]
        return '\n'.join(s)

    def __cmp__(self,other):
        if other.station < self.station:
            return 1
        if other.station > self.station:
            return -1
        if other.second < self.second:
            return 1
        if other.second > self.second:
            return -1
        if other.frame < self.frame:
            return 1
        if other.frame > self.frame:
            return -1
        return 0

    def follows(self,other):
        """ Return True iff the first frame of this fragment directly
            follows the last frame of "other", including threads.
        """

        return self.header.follows(other.last_header)

    def get_buff(self,nsamp=None,thread=0):
        if nsamp is None:
            nsamp = self.bufflen//2
        if thread < 0:
            return np.empty((2,nsamp),dtype=np.float32)
        return np.empty((1,nsamp),dtype=np.float32)
        
    def get_data(self,thread=0,nsamp=None,offs=None,buff=None):
        """ Get the unpacked data. 
        
        thread==0,1 single pol, thread==-1, both.
        offset: sample offset from beginning of data, must be provided
                for BOTH polarizations if requesting both

        NB if thread==-1, an offset for *both* polarizations must be
        specified.
        """

        d1 = np.memmap(self.fname,dtype=np.uint8)
        if nsamp is None:
            if thread != -1:
                nsamp = self.thread_nframe[thread]*self.frame_nsamp-offs
            else:
                nsamp = min(self.thread_nframe)*self.frame_nsamp-offs
        if offs is None:
            offs = [0,0]
        else:
            if thread == -1:
                try:
                    if not len(offs)==2:
                        raise ValueError('Must provide an offset for both polarizations!')
                except TypeError:
                    raise ValueError('Must provide an offset for both polarizations!')
            else:
                offs = [offs]*2
        offs = np.asarray(offs)


        # check that requested sample range lies within this data
        if thread == -1:
            if (nsamp+offs[0]) > self.get_nsamp(thread=0):
                raise ValueError('Requested sample range too large!')
            if (nsamp+offs[1]) > self.get_nsamp(thread=1):
                raise ValueError('Requested sample range too large!')
        else:
            if (nsamp+offs[thread]) > self.get_nsamp(thread=thread):
                raise ValueError('Requested sample range too large!')

        # determine which frames to load
        frame0 = min(offs)//self.frame_nsamp
        frame1 = (nsamp + max(offs))//self.frame_nsamp
        if (frame1*self.frame_nsamp) < nsamp+max(offs):
            frame1 += 1
        offs -= frame0*self.frame_nsamp

        # slice the correct range out of the memory map, then reshape it
        # to trim off the headers
        by_frames = d1[frame0*2*self.frame_length:(frame1+1)*2*self.frame_length]
        nframe = len(by_frames)//self.frame_length
        by_frames = by_frames.reshape(nframe,self.frame_length)[:,32:]

        # by_frames is now an array of 2N x 5000 chunks of data.  Divide
        # further by thread, then reshape back to an array of samples
        t0 = by_frames[::2].reshape((nframe//2*self.frame_nsamp))
        t1 = by_frames[1::2].reshape((nframe//2*self.frame_nsamp))
        if self.thread == 1:
            t0,t1 = t1,t0

        if buff is None:
            buff = self.get_buff(nsamp=nsamp,thread=thread)
        else:
            npol = 2 if thread==-1 else 1
            if buff.shape != (npol,nsamps):
                raise ValueError(
                        'Memory shape incompatible! (%s != %s)'%(
                            buff.shape,(npol,nsamps)))
        if thread == -1:
            buff[0,:] = t0[offs[0]:offs[0]+nsamp]
            buff[1,:] = t1[offs[1]:offs[1]+nsamp]
        elif thread == 0:
            buff[0,:] = t0[offs[0]:offs[0]+nsamp]
        else:
            buff[0,:] = t1[offs[1]:offs[1]+nsamp]
        buff -= 127.5
        return buff

class BasebandFragments(object):
    """ Encapsulate dumps for a single antenna."""

    def __init__(self,bbfs):
        """ bbfs -- a time-sorted list of BasebandFragment objects.
        """
        self.bbfs = sorted(bbfs)
        self.antenna = self.bbfs[0].antenna()
        nsamps_t0 = np.cumsum([bbf.get_nsamp(thread=0) for bbf in self.bbfs])
        nsamps_t1 = np.cumsum([bbf.get_nsamp(thread=1) for bbf in self.bbfs])
        self.nsamps_threads = [nsamps_t0,nsamps_t1]

    def __len__(self):
        return sum((len(bbf) for bbf in self.bbfs))

    def check_contiguity(self,fast=True):
        for b2,b1 in zip(self.bbfs[1:],self.bbfs[:-1]):
            if not b2.follows(b1):
                return False
        if not fast:
            for bbf in self.bbfs:
                if not bbf.check_contiguity():
                    return False
        return True

    def get_nframe(self):
        return sum((bbf.get_nframe() for bbf in self.bbfs))

    def __getitem__(self,index):
        return self.bbfs[index]

    def get_data(self,nsamp,offs=0,thread=0,buff=None):
        """ Return nsamp samples, offset by specified number of samples."""
        i0 = offs+nsamp
        i1 = gg
        bbf = self.bbfs[0]
        chunk_nsamp = bbf.bufflen//2
        if (samp_offset + nsamp) > chunk_nsamp*len(self.bbfs):
            raise ValueError('Requested chunk exceeds data range!')
        start_chunk = samp_offset // chunk_nsamp
        samp_offset -= start_chunk*chunk_nsamp
        nchunk = int(np.ceil(float(nsamp+samp_offset)/chunk_nsamp))
        if nchunk==1:
            skip_right = chunk_nsamp-nsamp-samp_offset
            return self.bbfs[start_chunk].get_data(
                    pol=pol,buff=buff,skip_left=samp_offset,skip_right=skip_right)
        else:
            if buff is None:
                buff = bbf.get_buff(nsamp)
            # do special case of first chunk
            self.bbfs[start_chunk].get_data(skip_left=samp_offset,pol=pol,buff=buff[:,:(chunk_nsamp-samp_offset)])
            for i in range(1,nchunk-1):
                self.bbfs[start_chunk+i].get_data(skip_left=0,pol=pol,buff=buff[:,i*chunk_nsamp-samp_offset:((i+1)*chunk_nsamp-samp_offset)])
            # do special case of last chunk
            remaining_samps = nsamp-(nchunk-1)*chunk_nsamp+samp_offset
            self.bbfs[start_chunk+nchunk-1].get_data(skip_left=0,skip_right=chunk_nsamp-remaining_samps,pol=pol,buff=buff[:,(nsamp-remaining_samps):])
        return buff


class Baseband(object):
    """ Encapsulate a set of baseband dumps.
    
    Facilitate online I/O for correlation, filtering, etc.
    """

    def __init__(self,fnames,fsamp=128000000,check_contiguity=True,
            samp_offsets=None):
        self.fsamp = fsamp
        self.fragments = sorted(map(BasebandFragment,fnames))
        self.by_antenna = dict()
        self.antennas = []
        for frag in self.fragments:
            ant = frag.antenna()
            if ant not in self.by_antenna.keys():
                self.by_antenna[ant] = []
                self.antennas.append(ant)
            self.by_antenna[ant].append(frag)
        self.bbfs = [BasebandFragments(self.by_antenna[ant]) for ant in self.antennas]
        frames = np.asarray([bbf[0].frame for bbf in self.bbfs])
        self.frame_offsets = frames.max()-frames
        #if enforce_order:
            #self._threadid_check()
        if check_contiguity:
            for bbf in self.bbfs:
                if not bbf.check_contiguity():
                    print('Warning! Antenna %s reports non-contiguous data.'%(bbf.antenna()))

    def _threadid_check(self):
        # check for threadid consistency
        threadids = np.asarray([bbf[0].threadid for bbf in self.bbfs])
        if not np.all(threadids==threadids[0]):
            raise ValueError('Different polarizations in starting frames!')

    def get_data(self,nsamp,samp_start=0,pol=0,antennas=None,buff=None,
            samp_offsets=None):
        """ NB pol=-1 gets both polarizations
        """
        if antennas is not None:
            bbfs,frame_offsets = zip(*[(bbf,fo) for bbf,fo in zip(self.bbfs,self.frame_offsets) if bbf.antenna in antennas])
        else:
            bbfs,frame_offsets = self.bbfs,self.frame_offsets
        if samp_offsets is not None:
            if len(samp_offsets) != len(bbfs):
                raise ValueError("Must have a sample offset for each antenna!")
        #samp_start = int(tstart*self.fsamp)
        samp_stop = samp_start + nsamp
        #samp_stop = int((tstart+tspan)*self.fsamp)
        #nsamp = samp_stop-samp_start
        npol = 1 + int(pol<0)
        nbeam = len(bbfs)*npol
        if buff is None:
            buff = np.empty((nbeam,nsamp),dtype=np.float32)
        else:
            if buff.shape != (nbeam,nsamp):
                raise ValueError("Memory shapes are incompatible!")
        for ibbf,(bbf,foff) in enumerate(zip(bbfs,frame_offsets)):
            samp_offset = foff*bbf[0].frame_nsamp
            bbf.get_data(nsamp=nsamp,samp_offset=samp_offset+samp_start,pol=pol,buff=buff[ibbf*npol:(ibbf+1)*npol,:])
        return buff

    def getlen(self):
        """ Return the data length, subject to alignment."""
        bufflen = self.bbfs[0].bbfs[0].bufflen
        frame_nsamp = self.bbfs[0].bbfs[0].frame_nsamp
        nchunks = len(self.bbfs[0].bbfs)
        return nchunks*bufflen//2-frame_nsamp*np.max(self.frame_offsets)

    def get_iterator(self,nsamp,samp_start=0,pol=0,antennas=None,
            samp_offsets=None):
        """ Set up a buffer/iterator object to iterate over the data."""
        return BasebandIterator(self,nsamp,samp_start=samp_start,pol=pol,antennas=antennas,samp_offsets=samp_offsets)

class BasebandIterator(object):
    """ Enumerate over baseband data from start to stop with given chunk
    and antenna size."""

    def __init__(self,bb,nsamp,samp_start=0,pol=0,antennas=None,
            samp_offsets=None):
        self.bb = bb
        self.fsamp = bb.fsamp
        self.nsamp = nsamp
        self.samp_start = samp_start
        self.pol = pol
        self.antennas = antennas
        self.nchunk = int(bb.getlen()//nsamp)
        self.samp_offsets = samp_offsets
        # initialize buffer
        self.buff = bb.get_data(nsamp,samp_start=samp_start,
                samp_offsets=samp_offsets,
                pol=self.pol,antennas=self.antennas)
        self._last_idx = None

    def __len__(self):
        return self.nchunk

    def __iter__(self):
        self.counter = 0
        return self

    def next(self):
        self.counter += 1
        if self.counter == self.nchunk:
            raise StopIteration
        samp_start = self.samp_start + (self.counter-1)*self.nsamp
        return self.bb.get_data(self.nsamp,samp_start=samp_start,
                pol=self.pol,antennas=self.antennas,buff=self.buff)

    def __getitem__(self,idx):
        if idx == self._last_idx:
            return self.buff
        if idx < 0:
            idx = idx + self.nchunk
        if (idx >= self.nchunk) or (idx < 0):
            raise ValueError
        samp_start = self.samp_start + idx*self.nsamp
        return self.bb.get_data(self.nsamp,samp_start=samp_start,
                pol=self.pol,antennas=self.antennas,buff=self.buff)

    def get_shape(self):
        return self.buff.shape

class OrderedBasebandIterator(object):
    """ Enumerate over baseband data from start to stop with given chunk
    and antenna size."""

    def __init__(self,stem,nsamp,path=None,samp_start=0,pol=0,
            antennas=None,sample_offsets=None,fsamp=128000000):
        if path is None:
            path = '/data/kerrm/voltages/%s'%(stem.split('_')[0])
        fnames = sorted(glob.glob('%s/%s*.bb'%(path,stem)))
        if antennas is not None:
            # check against file names
            def fname_okay(fname):
                ant = 'ea'+fname.split('ea')[1][:2]
                return ant in antennas
            fnames = filter(fname_okay,fnames)
        self.antennas = ['ea'+fname.split('ea')[1][:2] for fname in fnames]
        if sample_offsets is None:
            sample_offsets = np.zeros(len(self.antennas),dtype=int)
        else:
            if len(sample_offsets) != len(self.antennas):
                raise ValueError('Must have 1-1 mapping of offsets to antennas.')
        self.fnames = fnames
        self.stem = stem
        self.nsamp = nsamp
        self.samp_start = samp_start
        self.pol = pol
        self.sample_offsets = sample_offsets
        self.fsamp = fsamp

        self.mmaps = [np.memmap(fname,dtype=np.uint8,offset=off) for fname,off in zip(fnames,sample_offsets)]
        # allow for sample offsets
        lenpol = min((len(mmap)//2 for mmap in self.mmaps))
        if pol == -1:
            #self.pol_mmaps = zip(*[(x[0::2],x[1::2]) for x in self.mmaps])
            self.pol_mmaps = zip(*[(x[:lenpol],x[lenpol:]) for x in self.mmaps])
        elif pol == 0:
            self.pol_mmaps = [x[:lenpol] for x in self.mmaps]
        elif pol == 1:
            self.pol_mmaps = [x[lenpol:] for x in self.mmaps]
        else:
            raise ValueError('Invalid polarization specificiation ', pol)
        self.nchunk = lenpol//nsamp
        self.buff = np.empty((len(self.pol_mmaps),nsamp),dtype=np.float32)
        self._last_idx = None

    def new_iterator(self,samp_start=0,pol=0,sample_offsets=None):
        """ Return an iterator with similar characteristics but different
        polarization, starting position, offsets, etc.
        """
        d = dict(
                samp_start=samp_start,
                pol=pol,
                sample_offsets=sample_offsets,
                fsamp=self.fsamp,
                antennas=self.antennas)
        return OrderedBasebandIterator(self.stem,self.nsamp,**d)

    def __len__(self):
        return self.nchunk

    def __iter__(self):
        self.counter = 0
        return self

    def _load_data(self,samp_start,buff=None):
        if buff is not None:
            if buff.shape[0] != len(self.pol_mmaps):
                raise ValueError('Buffer not compatible with data.')
            nsamp = buff.shape[1]
        else:
            buff = self.buff
            nsamp = self.nsamp
        for immap,mmap in enumerate(self.pol_mmaps):
            buff[immap] = mmap[samp_start:samp_start+nsamp]
        np.subtract(buff,128,out=buff)
        return buff

    def get_data(self,nsamp,samp_start=0):
        """ An ad hoc method to return nsamp samples regardless of 
        iteration state.

        This should really be factored into a base class, but meh.
        """
        output = np.empty((self.buff.shape[0],nsamp),dtype=self.buff.dtype)
        return self._load_data(samp_start,buff=output)

    def next(self):
        self.counter += 1
        if self.counter == self.nchunk:
            raise StopIteration
        samp_start = self.samp_start + (self.counter-1)*self.nsamp
        return self._load_data(samp_start)

    def __getitem__(self,idx):
        if idx == self._last_idx:
            return self.buff
        if idx < 0:
            idx = idx + self.nchunk
        if (idx >= self.nchunk) or (idx < 0):
            raise ValueError
        samp_start = self.samp_start + idx*self.nsamp
        return self._load_data(samp_start)

    def nbeam(self):
        return self.buff.shape[0]

    def samp_per_chunk(self):
        return self.nsamp

    def nchan(self):
        nfft = self.buff.shape[1]
        return nfft//2+1

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

def re_order_baseband(fnames,max_frames=None):
    """ Take output dumps, which are sorted by frames and have somewhat
    arbitrary offsets, and write out the raw data with a common starting
    frame and appropriate length.  Order is [pol00...pol0N,pol10...pol1N].
    
    max_frames -- maximum frames (for a single polarization) to write
    """
    bb = Baseband(fnames)

    # use the headers from the first 100 frames to find find the frame
    # offset for each of the beams
    all_hdrs = [bbf[0].get_headers(100) for bbf in bb.bbfs]
    # find latest header at start and offsets to it
    ref = all_hdrs[np.argmax([hdrs[0] for hdrs in all_hdrs])][0]
    offsets = np.asarray([np.searchsorted(hdrs,ref) for hdrs in all_hdrs])
    if ref.threadid == 1:
        offsets += 1
    print(ref)
    print(ref.unique_idx())
    print(offsets)
    return
    maximum_frames = min((bbf.get_nframe() for bbf in bb.bbfs))-offsets.max()
    # always write out the same # of polarizations!
    maximum_frames -= maximum_frames % 2
    if max_frames is not None:
        maximum_frames = min(max_frames*2,maximum_frames)

    # for each beam write out maximum_frames frames to new array,
    # stripping headers
    for ibbfs,bbfs in enumerate(bb.bbfs):

        frames_to_write = maximum_frames
        frames_written = 0
        output = bbfs.bbfs[0].fname.split('buff')[0][:-1]+'.bb'
        outmmap = np.memmap(output,dtype=np.uint8,mode='w+',
                shape=((2,frames_to_write//2*bbfs.bbfs[0].frame_nsamp)))
        outmmap[:] = np.nan
        current_pol = 0
        out_offset = [0,0]
        in_off = offsets[ibbfs]

        for ibbf,bbf in enumerate(bbfs):
            inmmap = np.memmap(bbf.fname,dtype=np.uint8)
            #off = offsets[ibbfs] if ibbf==0 else 0
            by_frames = inmmap.reshape(bbf.nframe,bbf.frame_length)[in_off:,32:]
            # use no offset in all subsequent files
            in_off = 0
            for iframe,frame in enumerate(by_frames):
                if frames_written == frames_to_write:
                    break
                off = out_offset[current_pol]
                try:
                    outmmap[current_pol,off:off+bbf.frame_nsamp] = np.ravel(frame)
                except ValueError as e:
                    print(e)
                    print(iframe,frames_written,outmmap.shape,frames_to_write,off)
                    raise ValueError()
                out_offset[current_pol] += bbf.frame_nsamp
                frames_written += 1
                current_pol = 0 if current_pol==1 else 1

            del(inmmap)
            if frames_written == frames_to_write:
                break

        del(outmmap)

def get_delays_new(cmatrix,mask,frame_delays=None,normalize=True):
    """ Crude estimate of delays from correlation matrix."""
    nbeam = cmatrix.shape[0]
    nchan = cmatrix.shape[2]
    cmatrix = np.where(mask,cmatrix,0)
    delays = np.zeros((nbeam,nbeam))
    freqs = np.fft.fftfreq(nchan)
    for i in range(0,nbeam):
        for j in range(i+1,nbeam):
            c = cmatrix[i,j]
            if normalize:
                c /= (cmatrix[i,i].real*cmatrix[j,j].real)**0.5
            y = np.fft.fft(c)
            a = np.argmax(np.abs(y))
            delays[i,j] = freqs[a]
            delays[j,i] = -delays[i,j]
    return delays

def get_delays(cmatrix,chan_min=0,chan_max=None,chan_mask=None,
        frame_delays=None):
    """ Crude estimate of delays from correlation matrix.

        This is very susceptible to RFI so either provide an RFI-blanked
        matrix or set the channel boundaries (kwargs here) appropriately.

        See e.g. get_standard_mask.
    """

    if chan_mask is None:
        chan_mask = np.ones(cmatrix.shape[2],dtype=np.float32)
    if chan_max is not None:
        chan_mask[chan_max] = 0
    if chan_min is None:
        chan_mask[:chan_min+1] = 0
    if frame_delays is not None:
        assert(frame_delays.shape[0]==cmatrix.shape[0])
    else:
        frame_delays = np.zeros(cmatrix.shape[0])
    #cmatrix = cmatrix[...,chan_min:chan_max]
    delays = np.zeros((cmatrix.shape[0],cmatrix.shape[0]))
    freqs = np.fft.fftfreq(cmatrix.shape[2])
    for i in range(0,cmatrix.shape[0]):
        for j in range(i+1,cmatrix.shape[0]):
            c = cmatrix[i,j]/(cmatrix[i,i]*cmatrix[j,j])**0.5
            y = np.fft.fft(c*chan_mask)
            a = np.argmax(np.abs(y))
            delays[i,j] = freqs[a] + frame_delays[i] - frame_delays[j]
            delays[j,i] = -delays[i,j]
    return delays

def fit_delays(delays,guess=0):
    """ Determine the per-antenna delay that minimizes the r.m.s. of the
        predicted model delay on each baseline and that derived from the
        sine/cosine decomposition of the cross-powers.
    """
    def form_delays(delay_vector):
        delay_vector = np.append(0,delay_vector)
        d = np.zeros_like(delays)
        for i in range(len(delay_vector)):
            for j in range(len(delay_vector)):
                d[i,j] = delay_vector[i] - delay_vector[j]
        #for i in range(len(delay_vector)):
            #d[i] = delay_vector[i] - delay_vector
        return d
    idx = np.triu_indices(delays.shape[0],k=1)
    data = delays[idx]
    def chi(delay_vector):
        model = form_delays(delay_vector)
        return data-model[idx]
    guess = np.zeros(delays.shape[0]-1)+guess
    m = leastsq(chi,guess)[0]
    return np.append(0,m),form_delays(m)

def solve_delays(obi):
    """ Return a baseband iterator with approximate delay correction,
    along with the applied sample offsets.
    """
    cmatrix = correlate(obi,use_window=False,nchunk=1000)
    d0 = get_delays(cmatrix,chan_min=32000,chan_max=50000)
    m,model = fit_delays(d0)
    samp_offsets = np.round((m.max()-m)*1e-3*obi.fsamp).astype(int)
    print(samp_offsets)
    return obi.new_iterator(sample_offsets=samp_offsets)

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

def apply_filter(bbi,projector,nchunk=None,normalize=True):
    """
    Use baseband iterator and projection matrix to produce filtered and
    unfiltered power spectra.
    
    Started 12/7/2017.
    Last edited 12/7/2017.
    """

    ffti = FFTIterator(bbi,use_window=False)
    if nchunk is None:
        nchunk = len(bbi)
    nbeam = bbi.nbeam()
    nfft = bbi.buff.shape[1] 
    nchan = nfft//2+1

    if (projector.shape[0] != nchan):
        raise ValueError

    # set up correlation matrices (do both filtered and unfiltered)
    cmatrix0 = np.zeros((nbeam,nbeam,nchan),dtype=ffti.fbuff.dtype)
    cmatrix1 = np.zeros((nbeam,nbeam,nchan),dtype=ffti.fbuff.dtype)
    fbuff_conj = np.empty_like(ffti.fbuff)
    tmp_cmatrix = np.empty_like(cmatrix0)

    projector = projector.astype(ffti.fbuff.dtype)

    # calculate for each chunk
    for ichunk in range(nchunk):

        # compute unfiltered power spectrum
        fbuff = ffti[ichunk]
        np.conjugate(fbuff,out=fbuff_conj)
        np.multiply(fbuff[:,None,:],fbuff_conj,out=tmp_cmatrix)
        cmatrix0 += tmp_cmatrix

        # compute filtered power spectrum
        # proj has shape = nchan, nbeams, nbeams --> ijk
        # fbeams has shape = nbeams,nchan --> ki
        np.einsum('ijk,...ki->...ji',projector,fbuff,out=fbuff_conj)
        np.conjugate(fbuff_conj,out=fbuff)
        np.multiply(fbuff_conj[:,None,:],fbuff,out=tmp_cmatrix)
        cmatrix1 += tmp_cmatrix

    if normalize:
        norm = ffti.get_normalization()/nchunk
        np.multiply(cmatrix0,norm,out=cmatrix0)
        np.multiply(cmatrix1,norm,out=cmatrix1)

    return cmatrix0,cmatrix1


def real_to_complex(samples,shift_band=True,flip_sideband=False):
    """ Return a complex baseband (centred at 0 Hz) version of the real
    samples provided.
    
    shift_band -- if True, shift by fs/4 so that the two halves of the
    band are frequency ordered.  Otherwise, the +ve and -ve frequencies
    will be swapped relative to the physical band order.  False has the
    property of leaving the DC channel where it was originally sampled.

    flip_sideband -- conjugate the spectrum, shifting the freq. sense
    """

    # convert to analytic signal
    t = np.fft.fft(samples)
    t[1:len(t)//2] *= 2
    t[len(t)//2+1:] = 0

    # downsample
    x = np.empty(len(t)//2,dtype=np.complex64)
    x[:] = np.fft.ifft(t)[::2]

    if shift_band:
        # multiply by complex tone to shift to baseband
        # tone is exp(-i2pi * fs/4 * ts*2) = exp(-ipi) = [1, -1, 1, ...]
        x[1::2] *= -1

    if flip_sideband:
        x.imag *= -1

    return x

def fscrunch(filterbanks,n=3):
    if n == 0:
        return
    #if filterbanks.shape[-1] % 2**n != 0:
        #raise ValueError('Cannot evenly average channels.')
    single_axis = len(filterbanks.shape) == 1
    if single_axis:
        filterbanks = filterbanks[np.newaxis,:]
    for i in range(n):
        new_fb = np.empty((filterbanks.shape[0],filterbanks.shape[1]//2 + 1),
                dtype=filterbanks.dtype)
        new_fb[:,0] = filterbanks[:,0]
        new_fb[:,1:] = 0.5*(filterbanks[:,1::2] + filterbanks[:,2::2])
        filterbanks = new_fb
    if single_axis:
        filterbanks = filterbanks[0]
    return filterbanks

def filterbank(samples,nfft=12500*8,detect=True,favg=0,tavg=False):
    NFFT = nfft

    if detect:
        filterbanks = np.empty([samples.shape[0]//NFFT,NFFT//2+1],
                dtype=np.float64)
    else:
        filterbanks = np.empty([samples.shape[0]//NFFT,NFFT//2+1],
                dtype=np.complex128)

    for i in range(filterbanks.shape[0]):
        if detect:
            filterbanks[i] = np.abs(np.fft.rfft(samples[i*NFFT:(i+1)*NFFT]))**2
        else:
            filterbanks[i] = np.fft.rfft(samples[i*NFFT:(i+1)*NFFT])

    if not detect:
        return filterbanks

    for i in range(favg):
        new_fb = np.empty((filterbanks.shape[0],filterbanks.shape[1]//2 + 1),
                dtype=np.float64)
        new_fb[:,0] = filterbanks[:,0]
        new_fb[:,1:] = 0.5*(filterbanks[:,1::2] + filterbanks[:,2::2])
        filterbanks = new_fb

    if tavg:
        filterbanks = np.average(filterbanks,axis=0)

    return filterbanks

def get_cross_spectrum(s1,s2,nfft=50000,maxchunk=None):
    N = nfft
    nchunk = min(len(s1)//N,len(s2)//N)
    if maxchunk is not None:
        nchunk = min(maxchunk,nchunk)
    rvals = np.zeros(N//2,dtype=np.complex128)
    bp1 = np.zeros(N//2,dtype=np.float64)
    bp2 = np.zeros(N//2,dtype=np.float64)
    for i in range(nchunk):
        f1 = np.fft.rfft(s1[i*N:(i+1)*N])
        f2 = np.fft.rfft(s2[i*N:(i+1)*N])
        p = f1*f2.conj()
        rvals += p[1:]
        #bp1 += (f1*f1.conj()).real[1:]
        #bp2 += (f2*f2.conj()).real[1:]
        bp1 += np.abs(f1)[1:]
        bp2 += np.abs(f2)[1:]
    return rvals * (1./nchunk), bp1*(1./nchunk), bp2*(1./nchunk), nchunk

def get_lag_spectrum(s1,s2,nfft=50000):
    cc,bp1,bp2,nchunk = get_cross_spectrum(s1,s2,nfft=nfft)
    t = np.fft.fft(cc/(bp1*bp2)**0.5)
    return t/nchunk**2

def do_filter(beams,fbeams=None):

    nbeams = len(beams)
    tmat = 2 # integration time for correlation matrix

    # (0) select first 2*tmat s of data
    beams = [beam[:128000000*tmat*2] for beam in beams]

    # (1) get FFTs for each beam
    if fbeams is None:
        nfft = 128000
        fbeams = np.asarray(
                [filterbank(beam,nfft=nfft,detect=False) for beam in beams])
    else:
        nfft = 2*(fbeams[0].shape[1]-1)

    return fbeams

    nchan = nfft//2 + 1 # what to do with DC channel?

    # (2) now, build correlation matrix from first tmat of data
    nsamp = int(128e6*tmat/nfft)
    cmatrix = np.empty((nchan,len(beams),len(beams)),dtype=np.complex128)

    for ichan in range(nchan):
        if ichan % 1000 == 0:
            print(ichan)
        for iant in range(nbeams):
            di = fbeams[iant][:nsamp,ichan]
            for jant in range(iant,nbeams):
                dj = fbeams[jant][:nsamp,ichan].conj()
                cmatrix[ichan,iant,jant] = np.average(di*dj)
        for iant in range(nbeams):
            for jant in range(iant):
                cmatrix[ichan,iant,jant] = cmatrix[ichan,jant,iant].conj()

    # (3) assemble the orthogonal projector
    proj = np.empty((nchan,len(beams),len(beams)),dtype=np.complex128)
    dim = len(beams)
    for ichan in range(nchan):
        evals,evecs = np.linalg.eigh(cmatrix[ichan])
        v = np.matrix(evecs[np.argmax(evals)]).transpose()
        vt = v.transpose().conj()
        proj[ichan] = np.eye(dim)-(v/(vt*v))*vt

    # (4) project data
    ffilt = np.empty_like(fbeams)

    for ichan in range(nchan):
        ffilt[...,ichan] = np.matrix(proj[ichan])*np.matrix(fbeams[...,ichan])


    return cmatrix,proj,beams,fbeams,ffilt
        
def cyclic_spectrum(samples,nfft=12500*8,lags=500):
    # basically generate FFTs in chunks, and for each chunk, calculate a
    # lag spectrum

    # NB this version currently lacks twiddle factors to account for the
    # fact that the FFTs start at different times.  Since the frequency
    # lags are one channel, this is 1/2nfft.  The time shift is equal to
    # nfft, so each lag/block needs a twiddle factor
    # exp (-i*2pi*ilag/2nfft*iblock*nfft) =
    # exp (-i*pi*ilag*iblock).

    # since the algorithm works lag wise, this is equivalent to multiplying
    # each entry by either 1, -1, 1, -1, ... if ilag is odd, or
    #                      1,  1, 1,  1, ... if ilag is even

    rvals = np.zeros((lags,nfft//2+1),dtype=np.complex128)
    filterbanks = filterbank(samples,nfft=nfft,detect=False) 
    twiddles = np.ones((2,filterbanks.shape[0]))
    twiddles[1,1::2] = -1
    fbconj = filterbanks.conj()
    for lag in range(lags):
        rvals[lag,:] = np.average(filterbanks*np.roll(fbconj,lag,axis=1)*twiddles[lag%2,:,np.newaxis],axis=0)
    return rvals

def cyclic_spectrum_fsm(samples,nchan=64):
    """ Implement cyclic spectrum with a FFT + smoothing method.
    
    With this method, use the FFT to transform samples to freq domain.
    Then, for each channel, multiply with the conjugate to obtain an
    estimate of SC(nu,alpha) over the full range of alpha.  Do this for
    fbin channels and average the result over nu, yielding a top hat
    smoothed and decimated version of the cyclic spectrum.
    
    The overall complexity is O(N^2).
    """

    X = np.fft.rfft(samples).astype(np.complex64)
    Xc = X.conj()
    maxchan = int(5./64*len(X))
    output = np.zeros((nchan,maxchan),dtype=np.complex64)
    stride = len(X)//nchan
    print(maxchan,len(X),stride)
    for i in range(nchan):
        for j in range(stride):
            # multiply current channel by c.c. and roll such that alpha=0
            # is aligned with 0th bin
            #idx = i*stride+j
            output[i] += X[i*stride+j] * Xc[j:(j+maxchan)]
            #try:
            #    output[i] += X[idx] * Xc[idx:(idx+maxchan)]
            #except ValueError:
            #    continue
        Xc = np.roll(Xc,-stride)
        #Xc = np.roll(Xc,-i*stride)
    return output

    


def cyclic_spectrum2(samples,nfft=12500*8,lags=500):
    # attempt a time domain method

    pass
    rvals = np.zeros((lags,nfft//2+1),dtype=np.complex128)
    filterbanks = filterbank(samples,nfft=nfft,detect=False) 
    fbconj = filterbanks.conj()
    for lag in range(lags):
        print(lag)
        rvals[lag,:] = np.average(filterbanks*np.roll(fbconj,lag,axis=1),axis=0)
    return rvals

def filter_muos(samples,band=1,bw=5,decimate=False,use_window=False):
    """ Rotate baseband signal to put muos signal at 0 Hz, then filter.

    By default, return a 2x critically oversampled baseband signal.
    (Actually, at this point return full rate and leave it up to user how
    to decimate.)
    """
    # baseband signal is LSB, so MUOS offsets from 0 frequency are
    # band 4: 384-377.5 MHz.
    # band 3: 384-372.5 MHz.
    # band 2: 384-367.5 MHz.
    # band 1: 384-362.5 MHz.

    tone = 6.5*band # in MHz

    # first, get the complex baseband signal
    if not np.iscomplex(samples[0]):
        x = real_to_complex(samples,shift_band=False)
    else:
        x = samples

    # at this point, the MUOS RFI occupies the +ve frequencies of an FFT

    freq_shift = (4 + 5*(band-1))*1e6
    fsamp = 64e6
    tone_freq = np.array(-2j*np.pi*freq_shift/fsamp,dtype=np.complex64)

    x *= np.exp(tone_freq * np.arange(len(x),dtype=np.complex64))

    # at this point, the desired band is right above 0 Hz, at positive
    # frequencies, i.e. the signal is analytic

    if use_window:
        # this seems to screw things up
        raise NotImplementedError
        # TODO -- preserve normalization
        t = np.fft.fft(x*hamming(len(x)))
    else:
        t = np.fft.fft(x)

    # select bw
    nchan = int(float(bw)/64*len(t))
    if decimate:
        raise NotImplementedError
        t = t[:nchan]
    else:
        t[nchan:] = 0
    t = np.roll(t,-nchan//2)

    return np.fft.ifft(t)

def channelize(samples,nchan=64,nadv=32):
    """ Implement a moving-window filterbank preserving some of the native
    time resolution with overlap.
    """
    cmplx = np.iscomplex(samples[0])
    nsamps = nchan if cmplx else 2*nchan
    nspectra = (len(samples)-nsamps)//nadv
    window = np.hamming(nsamps)

    output = np.empty([nspectra,nchan+0 if cmplx else 1],dtype=np.complex64)
    for i in range(nspectra):
        output[i] = np.fft.rfft(samples[i*nadv:i*nadv+nsamps])
        #output[i] = pyfftw.interfaces.numpy_fft.rfft(samples[i*nadv:i*nadv+nchan*2])

    return output

def polyphase_filterbank(samples,nchan=64,nwindow=4):
    """ Implement a WOLA filterbank as an N-tap polyphase filter.
    """
    cmplx = np.iscomplex(samples[0])
    nsamps = nchan if cmplx else 2*nchan
    #window = hanning(nwindow*nsamps)
    window = hamming(nwindow*nsamps)
    #window = np.ones(nwindow*nsamps)
    norms = np.empty(nwindow)
    # this window normalization gets power normalization correct
    # normalize window
    norms = [1./np.sum(window[j*nsamps:(j+1)*nsamps]**2) for j in range(nwindow)]
    #window *= 2**0.5*(nwindow)**1.5/norm
    #window *= (2*nwindow)**0.5/(window**2).sum()
    #window *= 1./(window**2).sum()
    #print window.mean(),window.sum()
    nwindows = len(samples)//(nwindow*nsamps)-1
    nspectra = nwindows*nwindow
    output = np.empty((nspectra,nchan+int(not cmplx)),dtype=np.complex64)
    fft = np.fft.fft if cmplx else np.fft.rfft

    for i in range(nspectra):
        i0 = i*nsamps
        i1 = i0 + nwindow*nsamps
        tmp = window*samples[i0:i1]*norms[0]
        for j in range(1,nwindow):
            tmp[:nsamps] += tmp[j*nsamps:(j+1)*nsamps]*norms[j]
        output[i] = fft(tmp[:nsamps])
    #output *= 1./(2*nchan*nwindow)

    return output

def illustrative_spectrum(samples,nchan=3125):
    """ Plot PSD and cyclic spectra for 15 kHz, 120 kHz, 3.84 MHz and
    various other spectral lines."""

    # convert to real samples
    sc = real_to_complex(samples)

    # compute real power spectrum
    window = hanning(nchan)
    nspectra = len(sc)//nchan
    #psd = np.zeros(nchan,dtype=np.complex64)
    psd = np.empty((nspectra,nchan),dtype=np.float32)
    for i in range(nspectra):
        #psd += np.abs(np.fft.fft(window*sc[i*nchan:(i+1)*nchan]))
        psd[i]= np.abs(np.fft.fft(window*sc[i*nchan:(i+1)*nchan]))

    # compute shifted by 15 kHz
    dts = np.arange(len(sc),dtype=np.complex64)
    s1 = sc*np.exp((-2j*np.pi*0.5*15e3*(1./64e6))*dts)
    s2 = sc*np.exp((2j*np.pi*0.5*15e3*(1./64e6))*dts)
    #psd_15khz = np.zeros(nchan,dtype=np.complex64)
    psd_15khz = np.empty((nspectra,nchan),dtype=np.complex64)
    # do we need time twiddle factor?  This is basically TSM.
    for i in range(nspectra):
        f1 = np.fft.fft(window*s1[i*nchan:(i+1)*nchan])
        f2 = np.fft.fft(window*s2[i*nchan:(i+1)*nchan])
        psd_15khz[i] = f1*f2.conj()

    # compute shifted by 159... kHz (in principle completely uncorrelated)
    f = 159.0343343e3
    dts = np.arange(len(sc),dtype=np.complex64)
    s1 = sc*np.exp((-2j*np.pi*0.5*f*(1./64e6))*dts)
    s2 = sc*np.exp((2j*np.pi*0.5*f*(1./64e6))*dts)
    #psd_15khz = np.zeros(nchan,dtype=np.complex64)
    psd_uncorr = np.empty((nspectra,nchan),dtype=np.complex64)
    # do we need time twiddle factor?  This is basically TSM.
    for i in range(nspectra):
        f1 = np.fft.fft(window*s1[i*nchan:(i+1)*nchan])
        f2 = np.fft.fft(window*s2[i*nchan:(i+1)*nchan])
        psd_uncorr[i] = f1*f2.conj()

    return psd,psd_15khz,psd_uncorr

def do_chan(samples):
    nchunk = len(samples)//10
    chans = [channelize(samples[i*nchunk:(i+1)*nchunk],nchan=64,nadv=16) for i in range(10)]
    return chans

    # OK, arguably the best results I've gotten so far are to take these
    # channelized values, then compute the narrowband power for each chunk: 
    #ts = [chan[:,10]*chan[:,10].conj() for chan in all_chans
    # then take the Fourier transform:
    # tfs = np.asarray([np.fft.rfft(t) for t in ts])
    # then compute the spectral coherence by averaging the time chunks
    # *coherently:
    #pl.plot(np.fft.rfftfreq(len(ts[idx][1:]))*8,np.abs(tfs.mean(axis=0))[1:]/tfs[:,0].mean(),marker='.')
    # when doing this, all of the lines are only marginally resolved, so
    # one can imagine doing better with longer chunks
    # the 15 kHz line appears to be by far the best; using a window does
    # help somewhat with spectral leakage, but probably not critical


