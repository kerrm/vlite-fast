#include "process_baseband.h"

// quantities for D'Agostino normality test (see wikipedia)
#define NK float(NKURTO)
#define mu1 (-6./(NK+1))
#define mu2 ((24.*NK*(NK-2)*(NK-3))/((NK+1)*(NK+1)*(NK+3)*(NK+5)))
#define g1 (6.*(NK*NK-5*NK+2)/((NK+7)*(NK+9))*sqrt( (6.*(NK+3)*(NK+5))/(NK*(NK-2)*(NK-3)) ))
#define A (6.+(8./g1)*(2./g1 + sqrt(1. + 4./(g1*g1))))
#define Z2_1 sqrt(4.5*A)
#define Z2_2 (1-2./(9*A))
#define Z2_3 sqrt(2./(mu2*(A-4)))

#define NKb float(NFFT)
#define mu1b (-6./(NKb+1))
#define mu2b ((24.*NKb*(NKb-2)*(NKb-3))/((NKb+1)*(NKb+1)*(NKb+3)*(NKb+5)))
#define g1b (6.*(NKb*NKb-5*NKb+2)/((NKb+7)*(NKb+9))*sqrt( (6.*(NKb+3)*(NKb+5))/(NKb*(NKb-2)*(NKb-3)) ))
#define Ab (6.+(8./g1b)*(2./g1b + sqrt(1. + 4./(g1b*g1b))))
#define Z2b_1 sqrt(4.5*Ab)
#define Z2b_2 (1-2./(9*Ab))
#define Z2b_3 sqrt(2./(mu2b*(Ab-4)))

//convert unsigned char time array to float
__global__ void convertarray (cufftReal *time, unsigned char *utime, size_t n)
{
  for (int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; i += blockDim.x*gridDim.x)
  {
    if (utime[i] == 0) 
      time[i] = 0;
    else
      time[i] = (cufftReal)(utime[i])/128-1;
  }
}

__global__ void kurtosis (cufftReal *time, cufftReal *pow, cufftReal *kur)
{
  // calculate the variance (power) and kurtosis for voltage statistics in
  // relatively short windows.  Do this by first copying from global memory,
  // then using a hard-coded tree reduction.  Right now, it is set up to
  // use either 250 or 500 samples, so must be invoked with either 256 
  // or 512 threads.

  // because each thread block works on a chunk of data that's commensurate
  // with the packing of samples into the buffer, specifically with respect
  // to the two polarizations, I think we don't need to worry at all about
  // a thread block crossing the polarization.  The output will simply
  // contain the statistics for pol 0 first, then pol 1.

  volatile __shared__ float data2[256];
  volatile __shared__ float data4[256];
  unsigned int tid = threadIdx.x;
  size_t offset = blockIdx.x*NKURTO;

  if (tid < 250)
  {
    if (NKURTO==500) {
      // load up two values from global memory in this case
      data2[tid] = time[offset + tid]*time[offset + tid];
      float tmp = time[offset + tid + 250]*time[offset + tid + 250];
      data4[tid] = data2[tid]*data2[tid] + tmp*tmp;
      data2[tid] += tmp;
    }
    else {
      data2[tid] = time[offset + tid]*time[offset + tid];
      data4[tid] = data2[tid]*data2[tid];
    }
  }
  else
    data2[tid] = data4[tid] = 0;
  __syncthreads ();

  if (tid < 128)
  {
    data2[tid] += data2[tid + 128];
    data4[tid] += data4[tid + 128];
  }
  __syncthreads ();

  if (tid < 64)
  {
    data2[tid] += data2[tid + 64];
    data4[tid] += data4[tid + 64];
  }
  __syncthreads ();

  if (tid < 32)
  {
    data2[tid] += data2[tid + 32];
    data4[tid] += data4[tid + 32];
    data2[tid] += data2[tid + 16];
    data4[tid] += data4[tid + 16];
    data2[tid] += data2[tid + 8];
    data4[tid] += data4[tid + 8];
    data2[tid] += data2[tid + 4];
    data4[tid] += data4[tid + 4];
    data2[tid] += data2[tid + 2];
    data4[tid] += data4[tid + 2];
  }

  if (tid==0)
  {
    data2[tid] += data2[tid + 1];
    data4[tid] += data4[tid + 1];
    pow[blockIdx.x] = data2[0]/NKURTO;
    kur[blockIdx.x] = data4[0]/NKURTO/(pow[blockIdx.x]*pow[blockIdx.x]);
  }
}

__global__ void compute_dagostino (cufftReal* kur, cufftReal* dag, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    // I'm not sure why I have a zero check here; the only time it should
    // happen is if all of the samples are also 0.
    float dag1 = DAG_INF, dag2 = DAG_INF;
    if (kur[i] != 0.)
    {
      float t = (1-2./A)/(1.+(kur[i]-3.-mu1)*Z2_3);
      if (t > 0)
        dag1 = fabsf (Z2_1*(Z2_2 - powf (t,1./3)));
    }
    if (kur[i+n] != 0.)
    {
      float t = (1-2./A)/(1.+(kur[i+n]-3.-mu1)*Z2_3);
      if (t > 0)
        dag2 = fabsf (Z2_1*(Z2_2 - powf (t,1./3)));
    }
    // duplicate values to make bookkeeping in block_kurtosis easier
    dag[i] = dag[i+n] = fmaxf (dag1, dag2);
  }
}

// compute a filter-bank level statistic
// *** Importantly, this applies a fine-time filtering during calculation
// *** of statistic, by zero-weighting any NKURTO-sized blocks of samples
// *** that evade threshold.
__global__ void block_kurtosis (cufftReal* pow, cufftReal* kur, cufftReal* dag, cufftReal* pow_block, cufftReal* kur_block)
{
  volatile __shared__ float data2[256];
  volatile __shared__ float data4[256];
  volatile __shared__ unsigned char wt[256];
  // run with 256 threads; break it up such that we either do 5 blocks (for
  // NKURTO==500) or 10 blocks (for NKURTO=250)
  unsigned int tid = threadIdx.x;
  unsigned int warp_id = tid / 32;
  unsigned int warp_tid = tid - warp_id*32;
  if (warp_tid > 24) 
  {
    data2[tid] = 0;
    data4[tid] = 0;
    wt[tid] = 0;
  }
  else
  {
    // each thread block does 8 filterbank blocks (one for each warp)
    int idx = (blockIdx.x*8 + warp_id)*(NFFT/NKURTO) + warp_tid;
    //wt[tid] = (dag[idx]<DAG_THRESH) && (dag[idx]>-DAG_THRESH);
    // updated now that dag array is already absolute valued
    wt[tid] = dag[idx]<DAG_THRESH;
    data2[tid] = wt[tid]*pow[idx];
    data4[tid] = wt[tid]*kur[idx]*pow[idx]*pow[idx];
    if (NKURTO==250)
    {
      // if using finer time bins, add in the contribution from
      // the other pieces (see comment above)
      __syncthreads ();
      idx += 25;
      //float w = (dag[idx]<DAG_THRESH) && (dag[idx]>-DAG_THRESH);
      float w = dag[idx]<DAG_THRESH;
      data2[tid] += w*pow[idx];
      data4[tid] += w*kur[idx]*pow[idx]*pow[idx];
      wt[tid] += w;
    }
  }

  if (warp_tid > 15)
    return;

  // do sum within each warp
  data2[tid] += data2[tid + 16];
  data4[tid] += data4[tid + 16];
  wt[tid]    += wt[tid + 16];
  data2[tid] += data2[tid + 8];
  data4[tid] += data4[tid + 8];
  wt[tid]    += wt[tid + 8];
  data2[tid] += data2[tid + 4];
  data4[tid] += data4[tid + 4];
  wt[tid]    += wt[tid + 4];
  data2[tid] += data2[tid + 2];
  data4[tid] += data4[tid + 2];
  wt[tid]    += wt[tid + 2];
  data2[tid] += data2[tid + 1];
  data4[tid] += data4[tid + 1];
  wt[tid]    += wt[tid + 1];

  if (0==warp_tid)
  {
    if (wt[tid] > 0) 
    {
      float p = pow_block[blockIdx.x*8+warp_id] = data2[tid]/wt[tid];
      kur_block[blockIdx.x*8+warp_id] = data4[tid]/wt[tid]/(p*p);
    }
    else 
    {
      pow_block[blockIdx.x*8+warp_id] = 0;
      kur_block[blockIdx.x*8+warp_id] = 0;
    }
  }
}


// TODO -- this isn't quite right, because there won't necessarily be
// NFFT samples in the weighted version; since we're computing many fewer,
// don't need to precompute.  However, empirically it doesn't seem to make
// much of a difference..
__global__ void compute_dagostino2 (cufftReal* kur, cufftReal* dag, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    float dag1 = DAG_INF, dag2 = DAG_INF;
    if (kur[i] != 0.)
    {
      float t = (1-2./Ab)/(1.+(kur[i]-3.-mu1b)*Z2b_3);
      if (t > 0)
        dag1 = fabsf (Z2b_1*(Z2b_2 - powf (t,1./3)));
    }
    if (kur[i+n] != 0.)
    {
      float t = (1-2./Ab)/(1.+(kur[i+n]-3.-mu1b)*Z2b_3);
      if (t > 0)
        dag2 = fabsf (Z2b_1*(Z2b_2 - powf (t,1./3)));
    }
    dag[i] = dag[i+n] = fmaxf (dag1, dag2);
  }
}

__global__ void apply_kurtosis (
    cufftReal *in, cufftReal *out, 
    cufftReal *dag, cufftReal *dag_fb,
    cufftReal* norms)
{

  unsigned int tid = threadIdx.x;

  // D'Agostino kurtosis TS; already absolute valued and gathered
  // over the two polarizations, but is duplicated, so just use the
  // entry in the second polarization; will also make it easier if
  // we revert to independent polarizations
  //bool bad =  (dag[blockIdx.x] > DAG_THRESH) || (dag_fb[blockIdx.x/(NFFT/NKURTO)] > DAG_FB_THRESH);
  bool bad =  (dag[blockIdx.x] > DAG_THRESH);

  #ifdef DEBUG_WEIGHTS
  // if debugging, set the weights to 0 for the second half of all samples in
  // the chunk for 2nd pol and for the final eighth for the 1st pol
  int time_idx = blockIdx.x * NKURTO;
  bool c1 = time_idx > 3*(VLITE_RATE/(SEG_PER_SEC*2));
  bool c2 = (time_idx < VLITE_RATE/SEG_PER_SEC) && (time_idx > (7*VLITE_RATE/SEG_PER_SEC)/8);
  bad = c1 || c2;
  #endif

  if (bad)
  {
    // zero voltages
    if (tid < 250)
    {
      size_t offset = blockIdx.x*NKURTO;
      out[offset + tid] =  0;
      if (NKURTO==500)
        out[offset + tid + 250] =  0;
    }
  }
  else
  {
    // if copying data, copy it
    if (in != out && tid < 250)
    {
      size_t offset = blockIdx.x*NKURTO;
      out[offset + tid] = in[offset + tid];
      if (NKURTO==500)
        out[offset + tid + 250] = in[offset + tid + 250];
    }

    // add one to the filterbank block samples for weights
    if (tid==0)
    {
      atomicAdd (norms + (blockIdx.x*NKURTO)/NFFT, float(NKURTO)/NFFT);
    }
  }
}

__global__ void apply_kurtosis_fake (
    cufftReal *in, cufftReal *out, 
    cufftReal *dag, cufftReal *dag_fb,
    cufftReal* norms)
{

  unsigned int tid = threadIdx.x;

  if (in != out && tid < 250)
  {
    size_t offset = blockIdx.x*NKURTO;
    out[offset + tid] = in[offset + tid];
    if (NKURTO==500)
      out[offset + tid + 250] = in[offset + tid + 250];
  }

  // add one to the filterbank block samples for weights
  if (tid==0)
  {
    atomicAdd (norms + (blockIdx.x*NKURTO)/NFFT, float(NKURTO)/NFFT);
  }
}


__global__ void histogram ( unsigned char *utime, unsigned int* histo, size_t n)
{
  __shared__ unsigned int lhisto[512];
  lhisto[threadIdx.x] = 0;
  __syncthreads ();

  int i = threadIdx.x + blockIdx.x*blockDim.x;
  for (; i < n/2; i += blockDim.x*gridDim.x)
    atomicAdd (lhisto+utime[i], 1);
  for (; i < n; i += blockDim.x*gridDim.x)
    atomicAdd ((lhisto+256)+utime[i], 1);
  __syncthreads ();

  // MUST run with 512 threads for this global accumulation to work
  atomicAdd ( histo+threadIdx.x, lhisto[threadIdx.x]);
}

__global__ void set_frb_delays (float* frb_delays, float dm)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x; 
  if (i >= NCHAN) return;
  double freq = 0.384 - (i*0.064)/NCHAN;
  // delays are scaled by FFT timestep
  double scale = 4.15e-3*dm*SEG_PER_SEC*FFTS_PER_SEG;
  frb_delays[i] = float(scale/(freq*freq)-scale/(0.384*0.384));
}

__global__ void inject_frb ( cufftComplex *fft_out, float* frb_delays,
     int nfft_since_frb, float frb_width, float frb_amp)
{
  // NB frb_width must be in FFT time steps!

  // this is the channel; each thread does one channel for all time steps
  // and both polarizations
  int i = threadIdx.x + blockIdx.x*blockDim.x; 
  if (i >= NCHAN) return;
  
  // for now, don't try to do any interpolation, just round to the nearest
  // time index that the FRB encounters this channel

  int time_idx_lo = int(frb_delays[i]+0.5)-nfft_since_frb;
  int time_idx_hi = int(frb_delays[i]+frb_width+0.5)-nfft_since_frb;

  // if the earliest is after this chunk, return
  if (time_idx_lo >= FFTS_PER_SEG) return;

  // if the latest time precedes this chunk, return
  if (time_idx_hi < 0) return;

  // ensure indices are within data bounds
  if (time_idx_lo < 0) time_idx_lo = 0;
  if (time_idx_hi >= FFTS_PER_SEG) time_idx_hi = FFTS_PER_SEG-1;

  // otherwise, there is a portion of the FRB in this chunk, so loop over
  // the time steps that it passes through channel i
  for (int time_idx=time_idx_lo; time_idx<= time_idx_hi; time_idx++)
  {
    fft_out[time_idx*NCHAN+i].x *= frb_amp;
    fft_out[time_idx*NCHAN+i].y *= frb_amp;
  }

  // do the next polarization
  fft_out += FFTS_PER_SEG*NCHAN;

  for (int time_idx=time_idx_lo; time_idx<= time_idx_hi; time_idx++)
  {
    fft_out[time_idx*NCHAN+i].x *= frb_amp;
    fft_out[time_idx*NCHAN+i].y *= frb_amp;
  }

}

__global__ void detect_and_normalize2 (cufftComplex *fft_out, cufftReal* bp, 
        float scale)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x; 
  if (i >= NCHAN*2) return;
  if (i >= NCHAN) // advance pointer to next polarization
  {
    fft_out += FFTS_PER_SEG*NCHAN;
    bp += NCHAN;
    i -= NCHAN;
  }

  // initialize bandpass to mean of first block
  float bp_l = bp[i];
  if (0. == bp_l) {
    for (int j = i; j < FFTS_PER_SEG*NCHAN; j+= NCHAN)
      bp_l += fft_out[j].x*fft_out[j].x + fft_out[j].y*fft_out[j].y;
    bp_l /= FFTS_PER_SEG;
  }

  for (int j = i; j < FFTS_PER_SEG*NCHAN; j+= NCHAN)
  {
    // detect
    float pow = fft_out[j].x*fft_out[j].x + fft_out[j].y*fft_out[j].y;

    // update bandpass
    bp_l = scale*pow + (1-scale)*bp_l;

    // scale to bandpass and mean-subtract; this assumes the powers are 
    // chi^2_2 distributed, var(pow)=4, mean(pow)=std(pow)=2.  Therefore
    // dividing by mean will give standard deviation of 1 centred at 1.
    fft_out[j].x = pow/bp_l-1;

  }
  // write out current bandpass
  bp[i] = bp_l;
}

__global__ void detect_and_normalize3 (cufftComplex *fft_out, cufftReal* kur_weights_dev, cufftReal* bp, float scale)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x; 
  if (i >= NCHAN*2) return;
  if (i >= NCHAN) // advance pointer to next polarization
  {
    fft_out += FFTS_PER_SEG*NCHAN;
    kur_weights_dev += FFTS_PER_SEG;
    bp += NCHAN;
    i -= NCHAN;
  }

  // initialize bandpass to mean of first block
  float bp_l = bp[i];
  if (0. == bp_l) {
    int good_samples = 0;
    for (int j = i, time_idx=0; j < FFTS_PER_SEG*NCHAN; j+= NCHAN,time_idx++) {
      float w = kur_weights_dev[time_idx];
      if (0.==w)
        continue;
      good_samples++;
      bp_l += (fft_out[j].x*fft_out[j].x + fft_out[j].y*fft_out[j].y)/w;
    }
    if (0==good_samples) {
      // entire first block is bad; not sure what is best, try setting to
      // 1 and hope for the best?
      bp_l = 1;
    }
    else
      bp_l /= good_samples;
  }

  for (int j = i, time_idx=0; j < FFTS_PER_SEG*NCHAN; j+= NCHAN,time_idx++)
  {
    // detect
    //float w = kur_weights_dev[time_idx]*kur_weights_dev[time_idx];
    // NB that this formulation works because the weights are 0 or 1; if
    // we write out the expectation for the Fourier transform of the voltage
    // squared, the weights go in squared, so we normalize by the sum over
    // the weights squared, which is the same as the sum of the weights
    // (here kur_weights_dev) since they are 0 or 1
    float w = kur_weights_dev[time_idx];

    if (0.==w) {
      // if no samples are available, replace with mean bandpass
      fft_out[j].x = 0;
    }

    else {

      float pow = (fft_out[j].x*fft_out[j].x + fft_out[j].y*fft_out[j].y)/w;

      // apply a rough filter; values in excess of 11xmean shouldn't happen
      // more often than every 1.5 s, so we can clip values above this
      // without substantial distortion and possibly prevent bandpass
      // saturation; when we do, don't update the bandpass

      // TODO
      // NB this leads to a problem if we do allow in some RFI and the
      // bandpass gets stuck at a very small value.  Symptom is that the
      // output re-quantized bits are all maxval.

      if (pow > bp_l*11)
        fft_out[j].x = 10;

      else {

        // update bandpass
        bp_l = scale*pow + (1-scale)*bp_l;

        // scale to bandpass and mean-subtract; this assumes the powers are 
        // chi^2_2 distributed, var(pow)=4, mean(pow)=std(pow)=2.  Therefore
        // dividing by mean will give standard deviation of 1 centred at 1.
        fft_out[j].x = pow/bp_l-1;
      }
    }

  }
  // write out current bandpass
  bp[i] = bp_l;
}

// sum polarizations in place
__global__ void pscrunch (cufftComplex *fft_out, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    //fft_out[i].x += fft_out[i+n].x;
    fft_out[i].x = M_SQRT1_2*(fft_out[i].x + fft_out[i+n].x);
  }
}

// sum polarizations in place
__global__ void pscrunch_weights (cufftComplex *fft_out, cufftReal* kur_weights_dev, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    // this formulation excludes samples with more than 80% RFI
    float w1_f = kur_weights_dev[i/NCHAN];
    float w2_f = kur_weights_dev[(i+n)/NCHAN];
    int w1 = w1_f >= MIN_WEIGHT;
    int w2 = w2_f >= MIN_WEIGHT;
    switch (w1+w2)
    {
      case 2:
        // both samples OK, account for variance with sqrt(2)
        fft_out[i].x = M_SQRT1_2*(fft_out[i].x + fft_out[i+n].x);
        kur_weights_dev[i/NCHAN] = 0.5*(w1_f + w2_f);
        break;
      case 1:
        // only one sample OK, variance = 1
        fft_out[i].x = w1*fft_out[i].x + w2*fft_out[i+n].x;
        //kur_weights_dev[i/NCHAN] = 0.5*(w1_f*w1 + w2_f*w2);
        kur_weights_dev[i/NCHAN] = w1_f*w1 + w2_f*w2;
        break;
      case 0:
        // no good samples, average bandpass (NB isn't this just 0?)
        //fft_out[i].x = 0.5*(fft_out[i].x + fft_out[i+n].x);
        fft_out[i].x = 0.;
        kur_weights_dev[i/NCHAN] = 0;
        break;
    }
  }
}

// average time samples
// TODO -- review normalization and make sure it's correct with polarization
__global__ void tscrunch (cufftComplex *fft_out, cufftReal* fft_ave,size_t n)
{
  // loop over the output indices; calculate corresponding input index,
  // then add up the subsequent NSCRUNCH samples
  float scale = sqrt (1./NSCRUNCH);
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    // explicit calculation of indices for future reference
    ///////////////////////////////////////////////////////
    //int out_time_idx = i/NCHAN;
    //int out_chan_idx = i - out_time_idx*NCHAN;
    //int src_idx = out_time_idx * NSCRUNCH* NCHAN + out_chan_idx;
    //int src_idx = i+(NSCRUNCH-1)*out_time_idx*NCHAN;
    ///////////////////////////////////////////////////////
    int src_idx = i+(NSCRUNCH-1)*(i/NCHAN)*NCHAN;
    fft_ave[i] = 0.;
    for (int j=0; j < NSCRUNCH; ++j, src_idx += NCHAN)
    {
      fft_ave[i] += fft_out[src_idx].x;
    }
    fft_ave[i] *= scale;
  }
}

__global__ void tscrunch_weights (cufftComplex *fft_out, cufftReal* fft_ave, cufftReal* kur_weights_dev, size_t n)
{
  // loop over the output indices; calculate corresponding input index,
  // then add up the subsequent NSCRUNCH samples
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {

    // explicit calculation of indices for future reference
    ///////////////////////////////////////////////////////
    // int out_time_idx = i/NCHAN;
    // int out_chan_idx = i - out_time_idx*NCHAN;
    // int src_idx = out_time_idx * NSCRUNCH* NCHAN + out_chan_idx;
    // int src_idx = i+(NSCRUNCH-1)*out_time_idx*NCHAN;
    ///////////////////////////////////////////////////////

    // TODO -- we might not want an additional cut on MIN_WEIGHT here
    int src_idx = i+(NSCRUNCH-1)*(i/NCHAN)*NCHAN;
    fft_ave[i] = 0.;
    int wt_sum = 0;
    float wt_sumf = 0;
    for (int j=0; j < NSCRUNCH; ++j, src_idx += NCHAN)
    {
      float wt = kur_weights_dev[src_idx/NCHAN];
      if (wt < MIN_WEIGHT) continue;
      wt_sum++;
      wt_sumf += wt;
      fft_ave[i] += wt*fft_out[src_idx].x;
    }
    if (wt_sumf/NSCRUNCH >= MIN_WEIGHT)
      fft_ave[i] /= sqrt(float(wt_sum));
    else
      // this just copies the bandpass in; NB the average is needed, I'm not
      // entirely sure why
      //fft_ave[i] = fft_out[src_idx].x/NSCRUNCH;
      fft_ave[i] = 0;
  }
}

// select and digitize
__global__ void sel_and_dig_2b (
    cufftReal *fft_ave, unsigned char* fft_trim_u, size_t n, int npol)
{
  int NCHANOUT = CHANMAX-CHANMIN+1;
  //int NTIME = (n*4) / (NCHANOUT*npol); // total time samples
  int NTIME = (VLITE_RATE/SEG_PER_SEC/NFFT)/NSCRUNCH;
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    // compute index into input array
    // correct for packing of 4 (because i indexes a *byte*, not a sample)
    //int time_idx = (i*4)/(NCHANOUT);
    //int chan_idx = i*4 - time_idx*NCHANOUT;
    int time_idx = (i*4)/(NCHANOUT*npol);
    int pol_idx = (i*4 - time_idx*NCHANOUT*npol)/NCHANOUT;
    int chan_idx = i*4 - time_idx*npol*NCHANOUT - pol_idx*NCHANOUT;
    fft_trim_u[i] = 0;
    for (int j = 0; j < 4; ++j)
    {
      // I have now done an optimization of the input thresholds for the
      // approximate data format (chi^2 with 16 dof) assuming uniform
      // output.  This has about 5% more distortion than optimal output
      // with nonuniform steps, but is simpler for downstream applications.
      float tmp = fft_ave[pol_idx*NTIME*NCHAN + time_idx*NCHAN+chan_idx+CHANMIN+j];
      if (tmp < -0.6109) // do nothing, bit already correctly set
        continue;
      if (tmp < 0.3970)
        fft_trim_u[i] += 1 << 2*j;
      else if (tmp < 1.4050)
        fft_trim_u[i] += 2 << 2*j;
      else
        fft_trim_u[i] += 3 << 2*j;
    }
  }
}

// select and digitize
__global__ void sel_and_dig_4b (
    cufftReal *fft_ave, unsigned char* fft_trim_u, size_t n, int npol)
{
  int NCHANOUT = CHANMAX-CHANMIN+1;
  //int NTIME = n / (NCHANOUT*npol); // total time samples
  int NTIME = (VLITE_RATE/SEG_PER_SEC/NFFT)/NSCRUNCH;
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    // compute index into input array
    // correct for packing of 2 (because i indexes a *byte*, not a sample)
    //int time_idx = (i*2)/(NCHANOUT);
    //int chan_idx = i*2 - time_idx*NCHANOUT;
    int time_idx = (i*2)/(NCHANOUT*npol);
    int pol_idx = (i*2 - time_idx*NCHANOUT*npol)/NCHANOUT;
    int chan_idx = i*2 - time_idx*npol*NCHANOUT - pol_idx*NCHANOUT;
    // from Table 3 of Jenet & Anderson 1998
    //float tmp = fft_ave[time_idx*NCHAN+chan_idx+CHANMIN]/0.3188 + 7.5;
    float tmp = fft_ave[pol_idx*NTIME*NCHAN + time_idx*NCHAN+chan_idx+CHANMIN]/0.3188 + 7.5;
    if (tmp <= 0)
      fft_trim_u[i] = 0;
    else if (tmp >= 15)
      fft_trim_u[i] = 15;
    else
      fft_trim_u[i] = (unsigned char)(tmp);
    //tmp = fft_ave[time_idx*NCHAN+chan_idx+CHANMIN+1]/0.3188 + 7.5;
    tmp = fft_ave[pol_idx*NTIME*NCHAN + time_idx*NCHAN+chan_idx+CHANMIN+1]/0.3188 + 7.5;
    if (tmp <= 0)
      ;
    else if (tmp >= 15)
      fft_trim_u[i] += 15 << 4;
    else
      fft_trim_u[i] += (unsigned char)(tmp) << 4;
  }
}

// select and digitize
__global__ void sel_and_dig_8b (
    cufftReal *fft_ave, unsigned char* fft_trim_u, size_t n, int npol)
{
  int NCHANOUT = CHANMAX-CHANMIN+1;
  //int NTIME = n / (NCHANOUT*npol); // total time samples
  int NTIME = (VLITE_RATE/SEG_PER_SEC/NFFT)/NSCRUNCH;
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    // compute index into input array
    int time_idx = i/(NCHANOUT*npol);
    int pol_idx = (i - time_idx*NCHANOUT*npol)/NCHANOUT;
    int chan_idx = i - time_idx*npol*NCHANOUT - pol_idx*NCHANOUT;
    // from Table 3 of Jenet & Anderson 1998
    float tmp = fft_ave[pol_idx*NTIME*NCHAN + time_idx*NCHAN+chan_idx+CHANMIN]/0.02957 + 127.5;
    if (tmp <= 0)
      fft_trim_u[i] = 0;
    else if (tmp >= 255)
      fft_trim_u[i] = 255;
    else
      fft_trim_u[i] = (unsigned char) tmp;
  }
}

/*
// convert floating point to integer
__global__ void digitizearray(cufftReal *fft_ave, unsigned char* fft_ave_u, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    float tmp = fft_ave[i]/0.02957 + 127.5;
    if (tmp <= 0)
      fft_ave_u[i] = 0;
    else if (tmp >= 255)
      fft_ave_u[i] = 255;
    else
      fft_ave_u[i] = (unsigned char) tmp;
  }
}

// remove extraneous channels; in practice, this means 
__global__ void selectchannels(unsigned char* fft_ave_u, unsigned char* fft_trim_u, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    int nchan = CHANMAX-CHANMIN+1;
    int time_idx = i/nchan;
    int chan_idx = i - time_idx*nchan;
    fft_trim_u[i] = fft_ave_u[time_idx*NCHAN+chan_idx+CHANMIN];
  }
}

// detect total power and normalize the polarizations in place
// use a monolithic kernel pattern here since the total number of threads
// should be relatively small
// TODO -- this will probably be more efficient if done using lots of 
// threads and syncing; however, it doesn't seem to be a bottleneck
__global__ void detect_and_normalize (cufftComplex *fft_out, size_t ntime)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x; 
  if (i >= NCHAN*2) return;
  if (i >= NCHAN) // advance pointer to next polarization
  {
    fft_out += ntime*NCHAN;
    i -= NCHAN;
  }
  float sum1 = 0;
  float sum2 = 0;
  for (int j = i; j < ntime*NCHAN; j+= NCHAN)
  {
    float pow = fft_out[j].x*fft_out[j].x + fft_out[j].y*fft_out[j].y;
    fft_out[j].x = pow;
    sum1 += pow;
    sum2 += pow*pow;
  }
  sum1 *= 1./ntime;
  sum2 = sqrt(1./(sum2/ntime-sum1*sum1));
  for (int j = i; j < ntime*NCHAN; j+= NCHAN)
  {
    fft_out[j].x = (fft_out[j].x-sum1)*sum2;
  }
}

*/
