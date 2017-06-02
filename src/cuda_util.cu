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
