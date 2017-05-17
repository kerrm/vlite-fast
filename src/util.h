#include <cufft.h>
#include <curand.h>
#include <unistd.h>
#include <stdio.h>

void cudacheck (cudaError_t err);
void cufftcheck (cufftResult err);
void curandcheck (curandStatus_t err);
void dadacheck (int rcode);
FILE* myopen (const char* fname, const char* mode="rb", bool do_remove=false);
//
//output functions, modified from sigproc.h
void send_coords (double raj, double dej, double az, double za,FILE *output);
void send_double (const char *name, double double_precision, FILE *output);
void send_float (const char *name,float floating_point, FILE *output);
void send_int (const char *name, int integer, FILE *output);
void send_long (const char *name, long integer, FILE *output);
void send_string (const char *string, FILE *output);
int check_name (char*);
