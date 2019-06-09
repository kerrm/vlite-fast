#include <unistd.h>
#include <stdio.h>

void dadacheck (int rcode);
FILE* myopen (const char* fname, const char* mode="rb", bool do_remove=false, size_t bufsize=0);
//
//output functions, modified from sigproc.h
void send_coords (double raj, double dej, double az, double za,FILE *output);
void send_double (const char *name, double double_precision, FILE *output);
void send_float (const char *name,float floating_point, FILE *output);
void send_int (const char *name, int integer, FILE *output);
void send_long (const char *name, long integer, FILE *output);
void send_string (const char *string, FILE *output);
int check_name (char*);
int check_id (char*);
double coord_dist (double ra1, double ra2, double de1, double de2);
int check_coords (double raj, double dej, double tol=0.01);
void send_email(char* src_name, char* host_name);
