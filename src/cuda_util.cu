#include "cuda_util.h"
#include <stdio.h>

void cudacheck (cudaError_t err)
{
  if (err!=cudaSuccess)
  {
    fprintf (stderr, "CUDA ERROR! %s\n",cudaGetErrorString(err));
    //exit(0);
    throw 20;
  }
}

void cufftcheck(cufftResult err){
  if (err!=CUFFT_SUCCESS)
  {
    fprintf (stderr, "cufft ERROR! No. = %d\n", err);
    //exit(0);
    throw 20;
  }
}

void curandcheck(curandStatus_t err){
  if (err!=CURAND_STATUS_SUCCESS)
  {
    fprintf (stderr, "curand ERROR! No. = %d\n", err);
    //exit(0);
    throw 20;
  }
}

void print_cuda_properties (int device_idx, FILE* stream) {
  cudaDeviceProp prop;
  cudaGetDeviceProperties (&prop,device_idx);
  int nsms;
  cudaDeviceGetAttribute (&nsms,cudaDevAttrMultiProcessorCount,device_idx);
  fprintf (stream, "%s:\n",prop.name);
  fprintf (stream, "maxThreadsPerBlock= %d\n",prop.maxThreadsPerBlock);
  fprintf (stream, "maxThreadsDim[0]= %d\n",prop.maxThreadsDim[0]);
  fprintf (stream, "maxThreadsDim[1]= %d\n",prop.maxThreadsDim[1]);
  fprintf (stream, "maxThreadsDim[2]= %d\n",prop.maxThreadsDim[2]);
  fprintf (stream, "maxGristream, dSize[0]= %d\n",prop.maxGridSize[0]);
  fprintf (stream, "maxGridSize[1]= %d\n",prop.maxGridSize[1]);
  fprintf (stream, "maxGridSize[2]= %d\n",prop.maxGridSize[2]);
  fprintf (stream, "multiProcessors= %d\n",nsms);
}
