#include "util.h"
#include <string.h>
#include <stdlib.h>

void dadacheck (int rcode)
{
  if (rcode < 0)
  {
    fprintf (stderr, "dadacheck failed\n");
    throw 20;
  }
}


FILE* myopen(const char* fname, const char* mode, bool do_remove, size_t bufsize)
{
  if (do_remove && access(fname,F_OK)!=-1 )
    remove (fname);
  FILE* fp = fopen (fname, mode);
  if (NULL==fp)
  {
    fprintf (stderr, "Unable to open %s, aborting.\n",fname);
    exit (EXIT_FAILURE);
  }
  if (bufsize)
    setvbuf (fp, NULL, _IOFBF, bufsize);
  return (fp);
}

void myfwrite_byte (FILE* fp, void* buff, size_t nelem)
{
  if (fwrite (buff,1,nelem,fp) != nelem)
  {
    fprintf (stderr, "fwrite failed to write %ld bytes.\n",nelem);
    exit (EXIT_FAILURE);
  }
}

void myfwrite_float (FILE* fp, float* buff, size_t nelem)
{
  if (fwrite ((void*)buff,sizeof(float),nelem,fp) != nelem)
  {
    fprintf (stderr, "fwrite failed to write %ld floats.\n",nelem);
    exit (EXIT_FAILURE);
  }
}

void send_string(const char *string, FILE *output)
{
  int len;
  len=strlen(string);
  fwrite(&len, sizeof(int), 1, output);
  fwrite(string, sizeof(char), len, output);
}

void send_float(const char *name,float floating_point, FILE *output)
{
  send_string(name,output);
  fwrite(&floating_point,sizeof(float),1,output);
}

void send_double (const char *name, double double_precision, FILE *output)
{
  send_string(name,output);
  fwrite(&double_precision,sizeof(double),1,output);
}

void send_int(const char *name, int integer, FILE *output)
{
  send_string(name,output);
  fwrite(&integer,sizeof(int),1,output);
}

void send_long(const char *name, long integer, FILE *output)
{
  send_string(name,output);
  fwrite(&integer,sizeof(long),1,output);
}

void send_coords(double raj, double dej, double az, double za, _IO_FILE *output)
{
  if ((raj != 0.0) || (raj != -1.0)) send_double("src_raj",raj,output);
  if ((dej != 0.0) || (dej != -1.0)) send_double("src_dej",dej,output);
  if ((az != 0.0)  || (az != -1.0))  send_double("az_start",az,output);
  if ((za != 0.0)  || (za != -1.0))  send_double("za_start",za,output);
}

int check_name (char* src)
{
  return (strstr (src, "B0329+54") != NULL ||
          strstr (src, "B0531+21") != NULL || 
          strstr (src, "B0833-45") != NULL ||
          strstr (src, "B0950+08") != NULL || 
          strstr (src, "B1133+16") != NULL ||
          strstr (src, "B1237+25") != NULL ||
          strstr (src, "B1642-03") != NULL ||
          strstr (src, "B1749-28") != NULL || 
          strstr (src, "B1929+10") != NULL || 
          strstr (src, "J0341+5711") != NULL ||
          strstr (src, "J1713+0747") != NULL ||
          strstr (src, "J1909-3744") != NULL);
}
