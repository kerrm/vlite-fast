import numpy as np

class Candidate(object):

    def __init__(self,fil,line):
        self.fil = fil
        self.line = line.strip()
        toks = self.line.split()
        self.sn = float(toks[0])
        self.peak_idx = int(toks[1])
        self.peak_time = float(toks[2])
        # TODO -- is this OBO?
        self.tfilt = int(toks[3])
        self.dmi = int(toks[4])
        self.dm = float(toks[5])
        self.ngiant = int(toks[6])
        self.i0 = int(toks[7])
        self.i1 = int(toks[8])
        try:
            self.tsamp = FilReader(self.fil).header['tsamp']
        except:
            self.tsamp = 1./1280
        self.width = float(self.i1-self.i0)*self.tsamp
        self.sent_trigger = False

    def get_block(self,wmult=1,include_DM=True):
        fr = FilReader(self.fil)
        frh = fr.header
        width = self.i1-self.i0
        if include_DM:
            f0 = frh.fch1
            f1 = frh.nchans*frh.foff + frh.fch1
            dm_delay = 4.148808e3*self.dm*abs(f0**-2-f1**-2)
            dm_width = int(dm_delay/frh.tsamp)
        else:
            dm_width = 0
        start_samp = max(0,self.i0-width*wmult-dm_width)
        nsamp = self.i1+width*wmult+dm_width - start_samp
        block = fr.readBlock(start_samp,nsamp)
        return start_samp,block

    def tophat(self,block,override_tfilt=None):
        tfilt = override_tfilt or self.tfilt
        kernel = [1./2**tfilt]*2**tfilt
        if len(block.shape) == 2:
            kernel = [kernel]
        return fftconvolve(block,kernel,mode='same')

    def overlap(self,other,delta_dm=0.1,delta_w=3):
        """ Return true if DMs are close enough, the candidates overlap in time,
            and their widths aren't too discrepant.
        """
        if abs(self.dm/other.dm-1) > delta_dm:
            return 0
        w1,w2 = self.width,other.width
        if (w1 < w2):
            if w2/w1 > delta_w:
                return 0
        else:
            if w1/w2 > delta_w:
                return 0
        if self.i0 < other.i0:
            return (other.i0 < self.i1)
        return self.i0 < other.i1

    def __str__(self):
        return 'i0=%06d i1=%06d w=%3.2f sn=%3.2f dm=%3.2f'%(self.i0,self.i1,self.width*1000,self.sn,self.dm)

def coincidence (all_cands,delta_dm=0.1):

    # all_cands is a list of candidates from each beam 

    # assign each set of candidates a beam index (arbitrary)
    nbeam = len(all_cands)
    for ibeam in xrange(nbeam):
        for cand in all_cands[ibeam]:
            cand.beam = ibeam
            cand.beam_mask = np.zeros(nbeam,dtype=np.int16)

    # make a master list of all candidates
    all_cands = np.concatenate(all_cands)
    if len(all_cands)==0:
        return
    end_times = np.asarray([cand.i1 for cand in all_cands])
    a = np.argsort(end_times)
    all_cands = all_cands[a]
    end_times = end_times[a]*all_cands[0].tsamp

    # now go through each time slice (twice the maximum width) and correlate
    # candidates within it
    tslice = 1
    nslice = int(end_times[-1]/tslice) + 1
    idx0 = 0
    previous_cands = []

    for i in xrange(nslice):
        idx1 = np.searchsorted(end_times,tslice*(i+1))
        these_cands = all_cands[idx0:idx1]

        for cand in these_cands:
            # do correlation between candidates in this time slice
            for icand,ocand in enumerate(these_cands):
                overlap = cand.overlap(ocand)
                cand.beam_mask[ocand.beam] += overlap
                #ocand.beam_mask[cand.beam] += overlap
            # correlation with previous time slice (in case of straddling)
            for icand,ocand in enumerate(previous_cands):
                overlap = cand.overlap(ocand)
                cand.beam_mask[ocand.beam] += overlap
                #ocand.beam_mask[cand.beam] += overlap
        
        previous_cands = these_cands
        idx0 = idx1
    return all_cands

 
