#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>

#include "util.h"
#include "vdifio.h"

#include "dada_def.h"
#include "dada_hdu.h"
#include "ipcio.h"
#include "ascii_header.h"
#include "multilog.h"

#define VD_FRM 5032
#define VD_DAT 5000
#define VLITE_RATE 128000000
#define VLITE_FREQ 352.
#define VLITE_FRAME_RATE 25600


void usage ()
{
  fprintf(stdout,
      "Read a recorded baseband dump into a psrdada buffer "
      "(for testing / post-processing.\n"
      "Usage: readbase [filename]\n");
}

int main(int argc, char *argv[])
{
  char input_file[256];
  // check for mandatory input
  if (argc != 2)
    fprintf (stderr, "Must specify input file!.\n");
  snprintf (input_file, 256, argv[1]);

  //char hdr_buff[32];
  char input_buff[5032];
  //vdif_header* hdr = (vdif_header*) hdr_buff;

  // read first frame to get timestamp TODO?
  //size_t nread = fread (input_buff, 1, VD_FRM, input_fp);
  //if (nread == VD_FRM)
  //  ipcio_write (hdu->data_block,input_buff,VD_FRM);

  // connect to the output buffer
  key_t key = 0x40;
  multilog_t* log = multilog_open ("readbase",0);
  multilog_add (log, stderr);
  dada_hdu_t* hdu = dada_hdu_create (log);
  dada_hdu_set_key (hdu,key);
  dadacheck (dada_hdu_connect (hdu));

  // connect to DADA buffer and set current system time for epoch
  dadacheck (dada_hdu_lock_write (hdu));

  // write psrdada output header values
  char* ascii_hdr = ipcbuf_get_next_write (hdu->header_block);
  dadacheck (ascii_header_set (ascii_hdr, "NAME", "%s", "B0833-45" ) );
  dadacheck (ascii_header_set (ascii_hdr, "NCHAN", "%d", 1) );
  dadacheck (ascii_header_set (ascii_hdr, "BANDWIDTH", "%lf", 64) );
  dadacheck (ascii_header_set (ascii_hdr, "CFREQ", "%lf", 352) );
  dadacheck (ascii_header_set (ascii_hdr, "NPOL", "%d", 2) );
  dadacheck (ascii_header_set (ascii_hdr, "NBIT", "%d", 8) );
  dadacheck (ascii_header_set (ascii_hdr, "RA", "%lf", 0.87180) );
  dadacheck (ascii_header_set (ascii_hdr, "DEC", "%lf", 0.72452) );
  /*
  struct tm tm_epoch = {0};
  int vdif_epoch = getVDIFEpoch (hdr);
  tm_epoch.tm_year = 100 + vdif_epoch/2;
  tm_epoch.tm_mon = 6*(vdif_epoch%2)-1;
  time_t epoch_seconds = mktime (&tm_epoch) + getVDIFFrameEpochSecOffset (hdr);
  struct tm* utc_time = gmtime (&epoch_seconds);
  char dada_utc[64];
  strftime (dada_utc, 64, DADA_TIMESTR, utc_time);
  printf("UTC START: %s\n",dada_utc);
  dadacheck (ascii_header_set (ascii_hdr, "UTC_START", "%s", dada_utc) );
  */
  printf("%s",ascii_hdr);
  dadacheck (ipcbuf_mark_filled (hdu->header_block, 4096));

  FILE* input_fp = fopen (input_file, "rb");
  if (input_fp==NULL)
  {
    fprintf (stderr, "Unable to open %s, aborting.\n",input_file);
    exit (EXIT_FAILURE);
  }

  //size_t bytes_written = 0;
  while (1) {
    size_t nread = fread (input_buff, 1, VD_FRM, input_fp);
    if (nread != VD_FRM) 
      break;
    ipcio_write (hdu->data_block,input_buff,VD_FRM);
    //bytes_written += VD_FRM;
    //if (bytes_written % (VD_FRM*VLITE_FRAME_RATE) == 0)
    //  printf ("Read 1 s of data.\n");
  }
  dada_hdu_unlock_write (hdu);
  fclose (input_fp);
}


