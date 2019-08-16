#ifndef __PROCESS_BASEBAND_H__
#define __PROCESS_BASEBAND_H__

#include <cufft.h>

// GPU options
#define DEVICE 0
#define NTHREAD 512
#define PROFILE 0
#define CUDA_PROFILE_STOP(X,Y,Z) {\
cudaEventRecord (Y,0);\
cudaEventSynchronize (Y);\
cudaEventElapsedTime (Z,X,Y);}


#define VD_FRM 5032
#define VD_DAT 5000
#define VLITE_RATE 128000000
#define VLITE_FRAME_RATE 25600
#define NFFT 12500 // number of samples to filterbank
//#define NFFT 8192 // number of samples to filterbank
#define NCHAN (NFFT/2+1) //output filterbank channels, including DC
//#define NCHAN 4097 //output filterbank channels, including DC
#define NSCRUNCH 8 // time scrunch factor
//#define NSCRUNCH 5 // time scrunch factor
#define SEG_PER_SEC 10 // break each second of data up into chunks
//#define SEG_PER_SEC 5 // break each second of data up into chunks
#define FFTS_PER_SEG VLITE_RATE/SEG_PER_SEC/NFFT // filterbanks for a single pol

// write to multiple output buffers so we don't have to block on a
// completed observation; switch between them in ring fashion
//#define NOUTBUFF 2

// RFI excision and statistics recording options
#define NKURTO 500
#define DOHISTO 0
#define WRITE_KURTO 0
// TS should be normally distributed, so take 3-sigma as a threshold
#define DAG_THRESH 3.0
// these are more skew because they integrate over longer time, only
// excise the really bad ones
#define DAG_FB_THRESH 5.0
#define DAG_INF DAG_THRESH + DAG_FB_THRESH + 1
//#define DEBUG_WEIGHTS 1
# // exclude samples with more than 80% RFI
#define MIN_WEIGHT 0.2

//static volatile int CHANMIN=2411; // minimum output channel (counting DC), from low
//static volatile int CHANMIN=1; // minimum output channel (counting DC), from low
//static volatile int CHANMAX=6250; // maximum output channel (counting DC), from low
//__device__ int CHANMIN_d=1;
//__device__ int CHANMAX_d=1;
//#define CHANMIN 11 // minimum output channel (counting DC), from low
#define CHANMIN 2155 // minimum output channel (counting DC), from low/
#define CHANMAX 6250 // maximum output channel (counting DC), from low
//#define CHANMIN 1 // minimum output channel (counting DC), from low
//#define CHANMAX 4096 // maximum output channel (counting DC), from low

__global__ void convertarray (cufftReal *, unsigned char *, size_t);
__global__ void kurtosis (cufftReal *, cufftReal *, cufftReal *);
__global__ void compute_dagostino (cufftReal *, cufftReal* ,size_t);
__global__ void compute_dagostino2 (cufftReal *, cufftReal* ,size_t);
__global__ void block_kurtosis (cufftReal*, cufftReal*, cufftReal*, cufftReal*, cufftReal*);
__global__ void apply_kurtosis (cufftReal *, cufftReal *, cufftReal *, cufftReal*, cufftReal*);
__global__ void apply_kurtosis_fake (cufftReal *, cufftReal *, cufftReal *, cufftReal*, cufftReal*);
//__global__ void detect_and_normalize (cufftComplex *, size_t);
//
__global__ void set_frb_delays (float*, float);
__global__ void inject_frb ( cufftComplex *, float* , int , float , float );
__global__ void detect_and_normalize2 (cufftComplex *, cufftReal*, float);
__global__ void detect_and_normalize3 (cufftComplex *, cufftReal*, cufftReal*, float);
__global__ void histogram ( unsigned char *, unsigned int*, size_t);
__global__ void pscrunch (cufftComplex *, size_t);
__global__ void pscrunch_weights (cufftComplex *, cufftReal*, size_t);
__global__ void tscrunch (cufftComplex *, cufftReal *, size_t);
__global__ void tscrunch_weights (cufftComplex *, cufftReal* , cufftReal* , size_t);
__global__ void sel_and_dig_2b (cufftReal*, unsigned char*, size_t, int);
__global__ void sel_and_dig_4b (cufftReal*, unsigned char*, size_t, int);
__global__ void sel_and_dig_8b (cufftReal*, unsigned char*, size_t, int);

#endif
