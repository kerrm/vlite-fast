/* Current status:

   11/14/2018

   Seems to be working fine, but I note that there are some issues with
   buffer pollution for the case of smeared pulse widths comparable to the
   buffer width.  NOT equal, I find empirically that for the data to look
   reasonable in a filterbank analysis the buffer should contain several
   periods worth of data.  This is currently checked for below.  It would
   be nice to understand why it isn't the more reasonable 
   buffer > dm_smearing but there are only so many hours in the day.

*/
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#include <cufft.h>
#include <curand.h>
#include <curand_kernel.h>
#include <helper_cuda.h>

#include "vdifio.h"

#include "dada_def.h"
#include "dada_hdu.h"
#include "ipcio.h"
#include "ascii_header.h"
#include "multilog.h"

#include "util.h"
#include "cuda_util.h"

#define DEVICE 0        //GPU. on furby: 0 = TITAN Black, 1 = GTX 780
#define NTHREAD 512
#define WRITE_DADA 0 // 0 write to file, 1 write to psrdada buffer
#define NMOMENT 256

#define VD_FRM 5032
#define VD_DAT 5000
#define VLITE_RATE 128000000
#define VLITE_FREQ 352.
#define VLITE_FRAME_RATE 25600


__global__ void init_dm_kernel (cufftComplex *fker_dev, float dm, size_t n);
__global__ void set_profile (cufftReal *fdat_dev, size_t current_sample, 
    size_t period, int skip_period, float ampl, size_t n);
__global__ void multiply_kernel (cufftComplex* dat,cufftComplex* ker, size_t n);
__global__ void swap_sideband (cufftReal *dat, size_t n);
__global__ void setup_dstate (curandState *state);
__global__ void measure_moments (cufftReal* dat, float *moments);
__global__ void add_rfi (cufftReal *dat, size_t n, curandState *d_state,
    size_t current_sample, double tsamp_in_mus);
__global__ void digitize (float* idat, uint8_t* udat, size_t n);

void usage ()
{
  fprintf(stdout,"Usage: genbase [options]\n"
	  "-t seconds to simulate (default: 5)"
	  "-n observations to simulate (default: 1)\n"
	  "-p pulse period [s; default 0.5]"
    // this is used to simulate an FRB e.g. with -p 1.0 -k 100
	  "-k skip period [int; only produce a pulse every skip_period pulses]"
	  "-a amplitude as fraction of Tsys [default 0.05]"
	  "-s scale of second polarization relative to first [default 1.0]"
	  "-r seed for random number generator [long; default=42]"
	  "-d dm [default=30; NB this is about the largest feasible]"
	  "-e write to disk rather than psrdada buffer (default: false)"
	  "-f add RFI to observations (default: false)");
}

int main(int argc, char *argv[])
{
  double tobs = 5;
  double dm = 30;
  double pulse_period = 0.5;
  float ampls[2] = {0.05,0.05};
  float poln_ratio = 1.;
  int nobs = 1;
  int arg = 0;
  long seed = 42;
  int do_add_rfi = 0;
  int write_to_dada = 1;
  int skip_period = 1;
  while ((arg = getopt(argc, argv, "hfet:n:p:a:s:d:r:k:")) != -1) {

    switch (arg) {

    case 'h':
      usage ();
      return 0;
      
    case 't':
      if (sscanf (optarg, "%lf", &tobs) != 1) {
        fprintf (stderr, "genbase: could not read obs. time from %s\n", optarg);
        return -1;
      }
      break;

    case 'n':
      if (sscanf (optarg, "%d", &nobs) != 1) {
        fprintf (stderr, "writer: could not read num obs. from %s\n", optarg);
        return -1;
      }
      break;

    case 'a':
      if (sscanf (optarg, "%f", ampls) != 1) {
        fprintf (stderr, "genbase: could not read ampl.from %s\n", optarg);
        return -1;
      }
      break;

    case 's':
      if (sscanf (optarg, "%f", &poln_ratio) != 1) {
        fprintf (stderr, "genbase: could not pol'n ratio from %s\n", optarg);
        return -1;
      }
      break;

    case 'd':
      if (sscanf (optarg, "%lf", &dm) != 1) {
        fprintf (stderr, "genbase: could not read DM from %s\n", optarg);
        return -1;
      }
      break;

    case 'p':
      if (sscanf (optarg, "%lf", &pulse_period) != 1) {
        fprintf (stderr, "genbase: could not read period from %s\n", optarg);
        return -1;
      }
      break;

    case 'r':
      if (sscanf (optarg, "%li", &seed) != 1) {
        fprintf (stderr, "genbase: could not read seed from %s\n", optarg);
        return -1;
      }
      break;

    case 'f':
      printf ("genbase: adding RFI\n" );
      do_add_rfi = 1;
      break;

    case 'e':
      printf ("genbase: writing to disk\n" );
      write_to_dada = 0;
      break;


    case 'k':
      if (sscanf (optarg, "%d", &skip_period) != 1) {
        fprintf (stderr, "writer: could not read skip period from %s\n", optarg);
        return -1;
      }
      break;
    }
  }
  printf("Skip period =%d.\n",skip_period);

  // apply any polarization scaling
  ampls[1] = ampls[0] * poln_ratio;

  // set up sample counts for given DM
  double freq = 352;
  double freq_hi = 384;
  double freq_lo = 320;
  double tsamp = 1.0/VLITE_RATE; // NB real sampling
  printf("Sampling time is %g.\n",tsamp);
  double t_dm_lo = dm/2.41e-10*(1./(freq_lo*freq_lo)-1./(freq*freq)); // mus
  double t_dm_hi = dm/2.41e-10*(1./(freq*freq)-1./(freq_hi*freq_hi)); // mus
  printf ("DM smearing time to bottom of band is %.2f ms.\n",t_dm_lo*1e-3);
  printf ("DM smearing time to top of band is %.2f ms.\n",t_dm_hi*1e-3);
  printf ("deltaDM smearing time is %.2f ms.\n",(t_dm_lo-t_dm_hi)*1e-3);
  unsigned long n_dm_samp_lo = (unsigned long) t_dm_lo*1e-6/tsamp;
  unsigned long n_dm_samp_hi = (unsigned long) t_dm_hi*1e-6/tsamp;
  printf ("DM smearing samples to bottom of band is %li.\n",n_dm_samp_lo);
  printf ("DM smearing samples to top of band is %li.\n",n_dm_samp_hi);
  n_dm_samp_lo += (n_dm_samp_lo & 1); // make it even
  n_dm_samp_hi += (n_dm_samp_hi & 1); // make it even
  {
    unsigned long tmp = n_dm_samp_lo;
    n_dm_samp_lo = n_dm_samp_hi;
    n_dm_samp_hi = tmp;
  }
  unsigned long n_dm_samp = n_dm_samp_lo + n_dm_samp_hi;
  printf("DM total samples is %li.\n",n_dm_samp);

  // pulse properties
  size_t period_in_samples = size_t(pulse_period/tsamp);
  printf ("Pulse period is %li samples.\n",period_in_samples);

  // allocate memory for 1s of data; NB this is a large buffer, but because
  // of edge effects, will discard ~0.37s of data at DM=30!
  size_t buflen = VLITE_RATE/4;
  //size_t buflen = VLITE_RATE;

  if (buflen < 2*(n_dm_samp + period_in_samples))
  {
    fprintf (stderr, "Buffer not long enough to perform dedispersion!");
    exit (EXIT_FAILURE);
  }

  // initialize GPU properties
  cudacheck (cudaSetDevice (DEVICE));
  int nsms;
  cudaDeviceGetAttribute(&nsms,cudaDevAttrMultiProcessorCount,DEVICE);

  cufftReal* fdat_dev; cudacheck (
  cudaMalloc ((void**)&fdat_dev, sizeof(cufftReal)*buflen) );

  // make a separate buffer to store the overlap; need one for each polarization
  cufftReal* fovl_dev_p0; cudacheck (
  cudaMalloc ((void**)&fovl_dev_p0, sizeof(cufftReal)*n_dm_samp) );
  cufftReal* fovl_dev_p1; cudacheck (
  cudaMalloc ((void**)&fovl_dev_p1, sizeof(cufftReal)*n_dm_samp) );
  cufftReal* fovl_dev_pols[] = {fovl_dev_p0, fovl_dev_p1};

  // allocate memory for DM kernel
  cufftComplex* fker_dev; cudacheck (
  cudaMalloc ((void**)&fker_dev, sizeof(cufftComplex)*(buflen/2+1)) );

  // allocate memory for FFT
  cufftComplex* ffft_dev; cudacheck (
  cudaMalloc ((void**)&ffft_dev, sizeof(cufftComplex)*(buflen/2+1)) );

  // allocate memory for digitized output; only do the unpolluted samples
  size_t new_samps = buflen - n_dm_samp;
  uint8_t* udat_dev; cudacheck (
  cudaMalloc ((void**)&udat_dev, new_samps) );

  // allocate memory and initialize state for curand; NB that we need a
  // state for *each* thread
  curandState *d_state;
  if (do_add_rfi)
  {
    cudaMalloc(&d_state, sizeof(curandState)*NTHREAD*32);
    setup_dstate <<<32,NTHREAD>>> (d_state);
  }

  // allocate host memory; only copy over the unpolluted samples
  // leave room for a ragged edge VDIF frame
  size_t vdif_offset = 0;
  uint8_t* udat_host_p0; cudacheck (
  cudaMallocHost ((void**)&udat_host_p0, new_samps + VD_DAT) );
  uint8_t* udat_host_p1; cudacheck (
  cudaMallocHost ((void**)&udat_host_p1, new_samps + VD_DAT) );
  uint8_t* udat_host_pols[] = {udat_host_p0, udat_host_p1};

  // initialize DM kernel; also include correction for FFT normalization
  init_dm_kernel<<<32*nsms,NTHREAD>>> (fker_dev,dm,buflen/2+1);
  cudacheck ( cudaGetLastError() );

  // set up RNG for voltage generated
  curandGenerator_t gen;
  curandcheck (curandCreateGenerator (&gen, CURAND_RNG_PSEUDO_DEFAULT) );
  curandcheck (curandSetPseudoRandomGeneratorSeed (gen, seed) );

  // set up the FFTs; NB same plan for forward and backward
  cufftHandle plan_fwd,plan_bwd;
  checkCudaErrors (cufftPlan1d (&plan_fwd,buflen,CUFFT_R2C,1));
  checkCudaErrors (cufftPlan1d (&plan_bwd,buflen,CUFFT_C2R,1));

  // set up primary VDIF header
  char hdr_buff[32];
  vdif_header* hdr = (vdif_header*) hdr_buff;
  setVDIFBitsPerSample (hdr, 8);
  setVDIFFrameBytes (hdr, VD_FRM);
  setVDIFNumChannels (hdr,1);

  // if using PSRDADA, connect to the output buffer
  key_t key = 0x40;
  multilog_t* log=NULL;
  dada_hdu_t* hdu=NULL;
  FILE *output_fp=NULL;
  if (write_to_dada)
  {
    log = multilog_open ("genbase",0);
    hdu = dada_hdu_create (log);
    dada_hdu_set_key (hdu,key);
    dada_hdu_connect (hdu);
  }
  else
  {
    output_fp = myopen("/data/kerrm/baseband_sim.uw","wb");
  }


  // set time for current set of data; will change after generating apt.
  // no. of seconds
  for (int nseg=0; nseg < nobs; ++nseg) 
  {
  printf("Working on segment %d.\n",nseg);

  // initialize overlap buffer
  printf ("Setting up a buffer of %li with an overlap of %li.\n",buflen,n_dm_samp);
  printf ("Will lose %0.2f to edge effects.\n",double(n_dm_samp)/buflen);
  fflush (stdout) ;
  size_t current_sample = 0;

  //cufftReal* new_start_p0 = fdat_dev_p0  + n_dm_samp;
  //cufftReal* new_start_p1 = fdat_dev_p1  + n_dm_samp;
  for (int ipol=0; ipol < 2; ++ipol) {
    curandcheck (curandGenerateNormal (
        gen, (float*) fovl_dev_pols[ipol], n_dm_samp, 0, 1) );
    set_profile<<<32*nsms, NTHREAD>>> (

        fovl_dev_pols[ipol], current_sample, period_in_samples, 
        skip_period, 1+ampls[ipol], n_dm_samp);
    cudacheck ( cudaGetLastError() );
  }
  current_sample += n_dm_samp;

  if (write_to_dada)
  {
    // connect to DADA buffer and set current system time for epoch
    dada_hdu_lock_write (hdu);
    setVDIFFrameTime (hdr, time (NULL) );

    // write psrdada output header values
    char* ascii_hdr = ipcbuf_get_next_write (hdu->header_block);
    dadacheck (ascii_header_set (ascii_hdr, "NAME", "%s", "B0833-45" ) );
    dadacheck (ascii_header_set (ascii_hdr, "NCHAN", "%d", 1) );
    dadacheck (ascii_header_set (ascii_hdr, "BANDWIDTH", "%lf", 64.0) );
    dadacheck (ascii_header_set (ascii_hdr, "CFREQ", "%lf", 352.0) );
    dadacheck (ascii_header_set (ascii_hdr, "NPOL", "%d", 2) );
    dadacheck (ascii_header_set (ascii_hdr, "NBIT", "%d", 8) );
    dadacheck (ascii_header_set (ascii_hdr, "RA", "%lf", 0.87180) );
    dadacheck (ascii_header_set (ascii_hdr, "DEC", "%lf", 0.72452) );
    // NB psrdada format has TSAMP in microseconds
    //dadacheck (ascii_header_set (ascii_hdr, "TSAMP", "%lf", tsamp*1e6) );
    // set up epoch appropriately -- first, make a tm struct for VDIF epoch
    struct tm tm_epoch = {0};
    int vdif_epoch = getVDIFEpoch (hdr);
    tm_epoch.tm_year = 100 + vdif_epoch/2;
    tm_epoch.tm_mon = 6*(vdif_epoch%2);
    time_t epoch_seconds = mktime (&tm_epoch) + getVDIFFrameEpochSecOffset (hdr);
    struct tm* utc_time = gmtime (&epoch_seconds);
    char dada_utc[64];
    strftime (dada_utc, 64, DADA_TIMESTR, utc_time);
    printf("UTC START: %s\n",dada_utc);
    dadacheck (ascii_header_set (ascii_hdr, "UTC_START", "%s", dada_utc) );
    printf("%s",ascii_hdr);
    ipcbuf_mark_filled (hdu->header_block, 4096);
  }

  double sec_to_sim = tobs;
  size_t end_sample = size_t(sec_to_sim/tsamp);
  size_t current_frame = 0;
  int frame_seconds = 0;

  // sanity checks on voltage levels
  //float moments[2] = {0,0};
  //float* moments_dev;
  //cudaMalloc ((void**)&moments_dev, sizeof(float)*2);

  while (current_sample < end_sample)
  {
    
    for (int ipol=0; ipol < 2; ++ipol) {

      // copy overlap from previous input
      cudacheck (cudaMemcpy (
          fdat_dev, fovl_dev_pols[ipol], 
          n_dm_samp*sizeof(cufftReal), cudaMemcpyDeviceToDevice) );

      // generate input to fill non-overlap region; generate real samps
      curandcheck (curandGenerateNormal (
          gen, (float*) fdat_dev+n_dm_samp, new_samps, 0, 1) );

      // set pulse profile
      set_profile<<<32*nsms, NTHREAD>>> (
          fdat_dev+n_dm_samp, current_sample, period_in_samples,
          skip_period,1.+ampls[ipol], new_samps);
      cudacheck ( cudaGetLastError() );

      // copy input for next overlap to overlap buffer
      cudacheck (cudaMemcpy (
            fovl_dev_pols[ipol], fdat_dev+buflen-n_dm_samp, 
          n_dm_samp*sizeof(cufftReal), cudaMemcpyDeviceToDevice) );

      // forward transform the input
      cufftcheck (cufftExecR2C (plan_fwd, fdat_dev, ffft_dev) );

      // multiply by DM kernel
      multiply_kernel <<<32*nsms, NTHREAD>>> (ffft_dev, fker_dev, buflen/2+1);
      cudacheck ( cudaGetLastError() );

      // inverse transform
      cufftcheck (cufftExecC2R (plan_bwd, ffft_dev, fdat_dev) );

      // change to same sideband sense as VLITE
      swap_sideband <<<32*nsms, NTHREAD>>> (fdat_dev, buflen);
      cudacheck ( cudaGetLastError() );

      /*
      // an optional sanity check on the moments; last calculation showed
      // they nicedly followed a standard normal distribution
      moments[0] = 0;
      moments[1] = 0;
      cudacheck (cudaMemcpy (moments_dev, moments, 2*sizeof(float),
          cudaMemcpyHostToDevice));

      measure_moments <<< buflen/256, 256>>> (fdat_dev, moments_dev);
      cudacheck ( cudaGetLastError() );
      cudacheck (cudaMemcpy (moments, moments_dev, 2*sizeof(float), cudaMemcpyDeviceToHost) );
      moments[0] *= double(256)/buflen;
      moments[1] *= double(256)/buflen;
      printf ("moment2 %.6f %.6f\n", moments[0], sqrt(moments[0]));
      printf ("moment4 %.6f\n", moments[1]);
      */

      if (do_add_rfi)
      {
        // this sample offset means that RFI "phase" is referenced to the
        // very first sample written out; makes analysis easier
        add_rfi <<<32, NTHREAD>>> (fdat_dev, buflen, d_state, 
            current_sample-n_dm_samp-n_dm_samp_lo, tsamp*1e6);
        cudacheck ( cudaGetLastError() );
      }

      // digitize to 8-bit uints; simultaneously select only valid samples
      digitize <<<32*nsms, NTHREAD>>> (
          (float*)(fdat_dev+n_dm_samp_lo), udat_dev, new_samps);
      cudacheck ( cudaGetLastError() );

      // copy to host
      cudacheck (cudaMemcpy (udat_host_pols[ipol] + vdif_offset, udat_dev, new_samps, cudaMemcpyDeviceToHost) );
    } // end loop over polarizations

    current_sample += new_samps;

    // write to psrdada buffer or file
    size_t nframes = (new_samps + vdif_offset)/VD_DAT;

    if (write_to_dada)
    {
      for (size_t iframe = 0; iframe < nframes; ++iframe)
      {
        // update VDIF header
        if (current_frame == VLITE_FRAME_RATE)
        {
          frame_seconds ++;
          setVDIFFrameSecond (hdr, getVDIFFrameSecond (hdr) + 1);
          current_frame = 0;
        }
        setVDIFFrameNumber (hdr, current_frame);
        for (int ipol = 0; ipol < 2; ++ipol)
        {
          setVDIFThreadID(hdr, ipol);
          ipcio_write (hdu->data_block,hdr_buff,32);
          ipcio_write (hdu->data_block,(char*)(udat_host_pols[ipol] + VD_DAT*iframe), VD_DAT);
        }
        current_frame++;
      }
    }
    else
    {
      for (size_t iframe = 0; iframe < nframes; ++iframe)
      {
        // update VDIF header
        if (current_frame == VLITE_FRAME_RATE)
        {
          frame_seconds ++;
          setVDIFFrameSecond (hdr, getVDIFFrameSecond (hdr) + 1);
          current_frame = 0;
        }
        setVDIFFrameNumber (hdr, current_frame);
        for (int ipol = 0; ipol < 2; ++ipol)
        {
          setVDIFThreadID(hdr, ipol);
          fwrite(hdr_buff,1,32,output_fp);
          fwrite(udat_host_pols[ipol]+VD_DAT*iframe,1,VD_DAT,output_fp);
        }
        current_frame++;
      }
    }

    // copy remainder to beginning of buffer for next time
    // NB this is always a number <5032 samples, i.e. one frame
    // TODO -- figure this out and what to do RE polarization
    // MTK -- obviously I don't now know what the above TODO means
    size_t tocopy = new_samps + vdif_offset - nframes*VD_DAT;
    if (tocopy > 0)
    {
      cudacheck (cudaMemcpy (udat_host_p0, udat_host_p0 + nframes*VD_DAT, tocopy, cudaMemcpyHostToHost) );
      cudacheck (cudaMemcpy (udat_host_p1, udat_host_p1 + nframes*VD_DAT, tocopy, cudaMemcpyHostToHost) );
      vdif_offset = tocopy;
    }
   
  } // end loop over samples

  if (hdu)
    dada_hdu_unlock_write (hdu);
  if (output_fp)
    fclose (output_fp);

  } // end loop over observations

  // this cleanup a bit trivial at end of program
  curandDestroyGenerator (gen);
  cudaFree (fdat_dev);
  cudaFree (fovl_dev_p0);
  cudaFree (fovl_dev_p1);
  cudaFree (fker_dev);
  cudaFree (ffft_dev);
  cudaFree (udat_dev);
  cudaFreeHost (udat_host_p1);
  cudaFreeHost (udat_host_p0);

}



// set up DM kernel
__global__ void init_dm_kernel(cufftComplex *ker, float dm, size_t n)
{
  // i is the index into the array, and the FFT is arranged with frequencies
  // 0/N, 1/N, ...  (real to complex FFT, only +ve freqs)
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    // NB hardcoded bw and freq for now
    double freq = (64.*double(i))/double(n);
    double freq0 = 320.;
    double arg = (2*M_PI*dm/2.41e-10)*freq*freq/(freq0*freq0*(freq0+freq));
    double rcos,rsin;
    sincos(arg, &rsin, &rcos);
    ker[i].x = rcos/(2*(n-1));
    ker[i].y = rsin/(2*(n-1));

    // make a slightly more realistic bandpass; this has a relatively fast,
    // asymmetric taper on each side as well as a modest ramp
    freq *= 1./64;
    double scale = 1-exp(-(freq*freq)/(0.05*0.05));
    scale -= exp(-((1-freq)*(1-freq))/(0.10*0.10));
    scale *= (1+0.20*freq);
    ker[i].x *= scale;
    ker[i].y *= scale;
  }
}

__global__ void set_profile(cufftReal *dat, size_t current_sample, 
        size_t period_in_samples, int skip_period, float ampl, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    size_t sample = current_sample + i;
    /*
    // to emulate complex mixing, every odd sample should be multiplied by
    // -i, which is specific to the centre frequency / bandwidth of this app
    if (sample & 1)
    {
      float datx = dat[i].x;
      dat[i].x = dat[i].y;
      dat[i].y = -datx;
    }
    */
    // calculate integer phase for use in skip period
    int phasei = sample/period_in_samples;
    float phasef = float(sample-phasei*period_in_samples)/period_in_samples;
    if ((phasef < 0.03) && ((phasei%skip_period)==0))
    //if (phasef < 0.03)
    {
      //float tmp =1-abs(phase/0.025-1);
      //float amp = 1+tmp*tmp;
      //dat[i] *= amp;
      dat[i] *= ampl;
    }
  }
}

__global__ void multiply_kernel(cufftComplex* dat,cufftComplex* ker, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    float datx = dat[i].x;
    dat[i].x = datx*ker[i].x - dat[i].y*ker[i].y;
    dat[i].y = datx*ker[i].y + dat[i].y*ker[i].x;
  }
}

__global__ void measure_moments (cufftReal* dat, float *moments)
{
  unsigned int tid = threadIdx.x;
  size_t offset = blockIdx.x*NMOMENT;

  // general plan of work: do explicit sum within each warp
  volatile __shared__ float data2[256];
  volatile __shared__ float data4[256];

  data2[tid] = dat[offset + tid]*dat[offset + tid];
  data4[tid] = data2[tid]*data2[tid];

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
    atomicAdd (moments + 0, data2[0] / NMOMENT );
    atomicAdd (moments + 1, data4[0] / NMOMENT );
  }
}

__global__ void swap_sideband(cufftReal* dat, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    if (i & 1)
      dat[i] = -dat[i];
  }
}

__global__ void setup_dstate(curandState *state)
{

  int idx = threadIdx.x+blockDim.x*blockIdx.x;
  curand_init(idx+1233456, 0, 0, &state[idx]);
}

// Add roughly 1 mus of RFI every 10 mus of data.
__global__ void add_rfi(cufftReal* dat, size_t n, curandState *d_state, size_t current_sample, double tsamp_in_mus)
{
  curandState *state = &d_state[threadIdx.x + blockIdx.x*NTHREAD];
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {

    float phase = fmodf((i+current_sample) * (tsamp_in_mus/11.3),1);
    if (phase < 0.1)
    {
      // Quick and dirty, add a random uniform signal
      dat[i] += 5.*(curand_uniform (state) - 0.5);
    }
  }
}

__global__ void digitize(float* idat, uint8_t* udat, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    // add an extra 2 here for overhead in case we make it bright
    //float tmp = idat[i]/0.02957/2 + 127.5;
    // this normalization appears to be more consistent with the VLITE
    // digitizers, which have a mean of 128
    float tmp = idat[i]/0.02957/2 + 128.5;
    if (tmp <= 0)
      udat[i] = 0;
    else if (tmp >= 255)
      udat[i] = 255;
    else
      udat[i] = (uint8_t) tmp;
  }
}

