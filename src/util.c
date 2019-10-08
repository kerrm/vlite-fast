#include "util.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void dadacheck (int rcode)
{
  if (rcode < 0)
  {
    fprintf (stderr, "dadacheck failed\n");
    fflush (stderr);
    exit (EXIT_FAILURE);
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
          strstr (src, "J0332+54") != NULL || 
          strstr (src, "B0531+21") != NULL || 
          strstr (src, "J0534+22") != NULL || 
          strstr (src, "B2319+60") != NULL || 
          strstr (src, "J2321+6024") != NULL || 
          strstr (src, "B0833-45") != NULL ||
          strstr (src, "J0835-45") != NULL ||
          strstr (src, "B1237+25") != NULL ||
          //strstr (src, "B0950+08") != NULL || 
          //strstr (src, "B1133+16") != NULL ||
          //strstr (src, "B1237+25") != NULL ||
          //strstr (src, "B1642-03") != NULL ||
          //strstr (src, "B1749-28") != NULL || 
          //strstr (src, "B1929+10") != NULL || 
          //strstr (src, "3C147") != NULL || 
          //strstr (src, "3C48") != NULL || 
          //strstr (src, "J0341+5711") != NULL ||
          //strstr (src, "J1713+0747") != NULL ||
          strstr (src, "R2") != NULL ||
          strstr (src, "R3") != NULL);

}

int check_id (char* src)
{
  return (strstr (src, "18B-405") != NULL ||
          strstr (src, "19A-331") != NULL ||
          strstr (src, "SC1046")  != NULL
         );
}

double coord_dist (double ra1, double ra2, double de1, double de2)
{
  // NB everything is in radians!
  double dde = (de2-de1);
  double dra = (ra2-ra1)*cos(de1);
  return sqrt (dde*dde + dra*dra);
}

int check_coords (double raj, double dej, double tol)
{
  // NB everything is in radians!
  // positions of interest

  // position 1 -- arr2
  if (coord_dist(1.14479055, raj, 1.28572588, dej) < tol)
    return 1;
  //
  // position 1 -- arr3
  if (coord_dist(0.5110324, raj, 1.14737945, dej) < tol)
    return 1;

  // position 2 -- XTE 1809-197
  if (coord_dist(4.755373, raj, -0.344372, dej) < tol)
    return 1;

  return 0;
}

void send_email(char* src_name, char* host_name)
{
  // only send an email if executing on vlite-difx3; a better way to do this
  // is obviously with a daemon running on vlite-nrl or something, but OK
  if (strstr (host_name, "vlite-difx3") == NULL)
    return;
  static time_t t[5] = {0,0,0,0,0};
  int idx = 4;
  if (strstr (src_name, "R2") != NULL)
    idx = 0;
  else if (strstr (src_name, "3C147") != NULL)
    idx = 1;
  else if (strstr (src_name, "3C48") != NULL)
    idx = 2;
  else if (strstr (src_name, "R3") != NULL)
    idx = 3;
  // Lump all other sources into a single indx
  time_t new_time;
  time (&new_time);
  if ((new_time - t[idx]) > 1800)
  {
    // send email
    t[idx] = new_time;
    char cmd[256];
    snprintf (cmd,255,"echo '' | mail -s 'Source %s now being observed.' matthew.kerr@gmail.com,shining.surya.d8@gmail.com",src_name);
    system (cmd);
  }
}
