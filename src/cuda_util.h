#include <cufft.h>
#include <curand.h>
#include <stdio.h>

void cudacheck (cudaError_t err);
void cufftcheck (cufftResult err);
void curandcheck (curandStatus_t err);
void print_cuda_properties (int device_idx, FILE* stream);
