#include <cufft.h>
#include <curand.h>

void cudacheck (cudaError_t err);
void cufftcheck (cufftResult err);
void curandcheck (curandStatus_t err);
