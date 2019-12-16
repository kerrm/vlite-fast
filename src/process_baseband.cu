// system support
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <time.h>
#include <fcntl.h>

#include "vdifio.h"

// psrdada support
#include "dada_def.h"
#include "dada_hdu.h"
#include "ipcio.h"
#include "ascii_header.h"
#include "multilog.h"

// local support
#include "process_baseband.h"
#include "util.h"
#include "cuda_util.h"

// from Julia's code
extern "C" {
#include "utils.h"
#include "def.h"
#include "executor.h"
#include "multicast.h"
}

// rt profiling
#define RT_PROFILE 1

static volatile int NBIT = 2;
static FILE* logfile_fp = NULL;
static FILE* fb_fp = NULL;
static FILE* fb_kurto_fp = NULL;
static multilog_t* mlog = NULL;
static dada_hdu_t* hdu_in = NULL;
static dada_hdu_t* hdu_out = NULL;
static dada_hdu_t* hdu_co = NULL;
static int mc_control_sock = 0;


void usage ()
{
  fprintf (stdout,"Usage: process [options]\n"
	  "-k hexadecimal shared memory key for input (default: 40)\n"
	  "-K hexadecimal shared memory key for output (default: 0=disabled)\n"
	  "-C hexadecimal shared memory key for coadder (default: 0=disabled)\n"
//	  "-p listening port number (default: %zu; if 0, disable)\n"
	  "-o print logging messages to stdout (as well as logfile)\n"
	  "-w output filterbank data (0=no sources, 1=listed sources, 2=all sources [def])\n"
	  "-b reduce output to b bits (2 [def], 4, and 8 are supported))\n"
	  "-P number of output polarizations (1=AA+BB, 2=AA,BB; 4 not implemented)\n"
	  "-r RFI-excision mode (0=no excision, 1=only excision, 2=[default]write both data sets)\n"
	  "-i Inject an FRB at DM=80 periodicially in the output data.\n"
    "-s process a single observation, then quit (used for debugging)\n"
    "-p process a single recorded observation, then quit (used for profiling))\n"
    "-g run on specified GPU\n");
	  //"-m retain MUOS band\n"
//	  ,(uint64_t)READER_SERVICE_PORT);
}

void cleanup (void)
{
  fprintf (stderr,"called cleanup! [PROCESS_BASEBAND]\n");
  fflush (stderr);
  if (fb_fp) fclose (fb_fp);
  if (fb_kurto_fp) fclose (fb_kurto_fp);
  fprintf (stderr,"h1\n");
  fflush (stderr);
  if (hdu_in != NULL)
  {
    dada_hdu_disconnect (hdu_in);
    dada_hdu_destroy (hdu_in);
  }
  fprintf (stderr,"h2\n");
  fflush (stderr);
  if (hdu_out != NULL)
  {
    dada_hdu_disconnect (hdu_out);
    dada_hdu_destroy (hdu_out);
  }
  fprintf (stderr,"h3\n");
  fflush (stderr);
  if (hdu_co)
  {
    dada_hdu_disconnect (hdu_co);
    dada_hdu_destroy (hdu_co);
  }
  fprintf (stderr,"h4\n");
  fflush (stderr);
  if (mc_control_sock > 0)
    shutdown (mc_control_sock, 2);
  fprintf (stderr,"h5\n");
  fflush (stderr);
  multilog (mlog, LOG_INFO, "Completed shutdown [PROCESS_BASEBAND].\n");
  //multilog_close (mlog);
  fprintf (stderr,"h6\n");
  fflush (stderr);
  if (logfile_fp) fclose (logfile_fp);
  fprintf (stderr,"h7\n");
  fflush (stderr);
}

void exit_handler (void) {
  fprintf (stderr, "exit handler called\n");
  fflush (stderr);
  cleanup ();
}

void sigint_handler (int dummy) {
  cleanup ();
  exit (EXIT_SUCCESS);
}


void change_extension (const char* in, char* out, const char* oldext, const char* newext)
{
  char* substr;
  strcpy (out, in);
  substr = strstr (out, oldext);
  if (NULL==substr)
    strcat (out,newext);
  else
    snprintf (substr, 256-strlen(out)-1, newext);  
}

/*
TODO generalize above to allow changing a substring, particularly to change
"muos" to "muos_kur" and then make this the default for output.
*/

// TODO -- check that this is consistent with sigproc
int write_psrdada_header (dada_hdu_t* hdu, char* inchdr, vdif_header* vdhdr, int npol, char* fb_file)
{

  // handle values from incoming header; set defaults if not available
  int station_id = 0;
  ascii_header_get (inchdr, "STATIONID", "%d", &station_id);
  double ra = 0;
  ascii_header_get (inchdr, "RA", "%lf", &ra);
  double dec = 0;
  ascii_header_get (inchdr, "DEC", "%lf", &dec);
  char name[OBSERVATION_NAME_SIZE];
  ascii_header_get (inchdr, "NAME", "%s", name);
  name[OBSERVATION_NAME_SIZE-1] = '\0';
  double start = 0;
  ascii_header_get (inchdr, "SCANSTART", "%lf", &start);

  // update the time with the actual data start, since we have discarded
  // some data to reach a 1s boundary
  time_t epoch_seconds = vdif_to_unixepoch (vdhdr);
  struct tm utc_time;
  gmtime_r (&epoch_seconds, &utc_time);
  char dada_utc[DADA_TIMESTR_LENGTH];
  strftime (dada_utc, DADA_TIMESTR_LENGTH, DADA_TIMESTR, &utc_time);

  // initialize observation parameters for filterbank
  // NB the data are upper sideband, so negative channel bandwidth
  double chbw = -64./NCHAN;
  double tsamp = double(NFFT)/VLITE_RATE*NSCRUNCH*1e6; // NB in mus
  int nchan = CHANMAX-CHANMIN+1;
  double bw = nchan*chbw;
  double freq0 = 384.;
  double freq = freq0 + 0.5*(CHANMIN+CHANMAX-1)*chbw;

  fprintf (stderr, "before lock\n");
  dadacheck (dada_hdu_lock_write (hdu));
  fprintf (stderr, "after lock\n");
  char* ascii_hdr = ipcbuf_get_next_write (hdu->header_block);
  fprintf (stderr, "after next write\n");
  dadacheck (ascii_header_set (ascii_hdr, "STATIONID", "%d", station_id));
  dadacheck (ascii_header_set (ascii_hdr, "BEAM", "%d", station_id));
  dadacheck (ascii_header_set (ascii_hdr, "RA", "%lf", ra));
  dadacheck (ascii_header_set (ascii_hdr, "DEC", "%lf", dec));
  dadacheck (ascii_header_set (ascii_hdr, "NAME", "%s", name));
  dadacheck (ascii_header_set (ascii_hdr, "SCANSTART", "%lf", start));
  dadacheck (ascii_header_set (ascii_hdr, "NCHAN", "%d", nchan) );
  dadacheck (ascii_header_set (ascii_hdr, "BANDWIDTH", "%lf", bw) );
  dadacheck (ascii_header_set (ascii_hdr, "CFREQ", "%lf", freq) );
  dadacheck (ascii_header_set (ascii_hdr, "NPOL", "%d", npol) );
  dadacheck (ascii_header_set (ascii_hdr, "NBIT", "%d", NBIT) );
  dadacheck (ascii_header_set (ascii_hdr, "TSAMP", "%lf", tsamp) );
  dadacheck (ascii_header_set (ascii_hdr, "UTC_START", "%s", dada_utc) );
  double epc = epoch_seconds;
  dadacheck (ascii_header_set (ascii_hdr, "UNIXEPOCH", "%lf", epc) );
  // also record the VDIF MJD info, this is useful for finding
  // transients in the baseband stream.
  dadacheck (ascii_header_set (ascii_hdr, "VDIF_MJD", "%d", 
      getVDIFFrameMJD (vdhdr)) );
  dadacheck (ascii_header_set (ascii_hdr, "VDIF_SEC", "%lu", 
      getVDIFFrameMJDSec (vdhdr)) );

  if (fb_file)
    dadacheck (ascii_header_set (ascii_hdr, "SIGPROC_FILE", "%s", fb_file) );
  multilog (hdu->log, LOG_INFO, "%s",ascii_hdr);
  ipcbuf_mark_filled (hdu->header_block, 4096);
  return 0;
}

int  clear_psrdada_buffer (dada_hdu_t* hdu)
{
  fprintf (stderr, "locking read");
  dadacheck (dada_hdu_lock_read (hdu));
  uint64_t hdr_size = 0;
  fprintf (stderr, "clearing headers");
  for (uint64_t i = 0; i < ipcbuf_get_nfull (hdu->header_block); ++i)
  {
    ipcbuf_get_next_read (hdu->header_block,&hdr_size);
    dadacheck (ipcbuf_mark_cleared (hdu->header_block));
  }
  ipcbuf_t* db = (ipcbuf_t*)hdu->data_block;
  fprintf (stderr, "clearing datas");
  for (uint64_t i = 0; i < ipcbuf_get_nfull (db); ++i)
  {
    ipcbuf_get_next_read (db,&hdr_size);
    dadacheck (ipcbuf_mark_cleared (db));
  }
  fprintf (stderr, "cleared data\n");
  dadacheck (dada_hdu_unlock_read (hdu));
  return 0;
}

int write_sigproc_header (FILE* output_fp, char* inchdr, vdif_header* vdhdr,int npol)
{
  // handle values from incoming header; set defaults if not available
  int station_id = 0;
  ascii_header_get (inchdr, "STATIONID", "%d", &station_id);
  double ra = 0;
  ascii_header_get (inchdr, "RA", "%lf", &ra);
  double dec = 0;
  ascii_header_get (inchdr, "DEC", "%lf", &dec);
  char name[OBSERVATION_NAME_SIZE];
  ascii_header_get (inchdr, "NAME", "%s", name);
  name[OBSERVATION_NAME_SIZE-1] = '\0';

  double chbw = -64./NCHAN;
  double tsamp = double(NFFT)/VLITE_RATE*NSCRUNCH;
  double mjd = getVDIFFrameDMJD (vdhdr, VLITE_FRAME_RATE);
  // write out a sigproc header
  send_string ("HEADER_START",output_fp);
  send_string ("source_name",output_fp);
  send_string (name,output_fp);
  send_int ("barycentric",0,output_fp);
  send_int ("telescope_id",station_id,output_fp);
  // VLA coordinates are in radians; convert to sigproc
  float hh = (180/M_PI)*(24./360)*ra;
  float mm = (hh-int(hh))*60;
  float ss = (mm-int(mm))*60;
  float sigproc_ra = int(hh)*1e4 + int(mm)*1e2 + ss;
  send_double ("src_raj",sigproc_ra,output_fp);
  float dd = (180/M_PI)*fabs(dec);
  mm = (dd-int(dd))*60;
  ss = (mm-int(mm))*60;
  float sigproc_dec = int(dd)*1e4 + int(mm)*1e2 + ss;
  if (dec < 0) dec = -dec;
  send_double ("src_dej",sigproc_dec,output_fp);
  send_int ("data_type",1,output_fp);
  send_double ("fch1",384+(CHANMIN-0.5)*chbw,output_fp);
  send_double ("foff",chbw,output_fp);//negative foff, fch1 is highest freq
  send_int ("nchans",CHANMAX-CHANMIN+1,output_fp);
  send_int ("nbits",NBIT,output_fp);
  send_double ("tstart",mjd,output_fp);
  send_double ("tsamp",tsamp,output_fp);//[sec]
  send_int ("nifs",npol,output_fp);
  send_string ("HEADER_END",output_fp);
  return station_id;
}

void get_cofbfile (char* fbfile, ssize_t fbfile_len, char* inchdr, vdif_header* vdhdr)
{
  // Open up filterbank file using timestamp and antenna
  int station_id = 99;
  //ascii_header_get (inchdr, "STATIONID", "%d", &station_id);
  char currt_string[128];
  time_t epoch_seconds = vdif_to_unixepoch(vdhdr);
  struct tm utc_time;
  gmtime_r (&epoch_seconds, &utc_time);
  strftime (currt_string,sizeof(currt_string), "%Y%m%d_%H%M%S", &utc_time);
  *(currt_string+15) = 0;
  if (CHANMIN < 2411)
    snprintf (fbfile,fbfile_len,"%s/%s_muos_ea%02d.fil",DATADIR,currt_string,station_id);
  else
    snprintf (fbfile,fbfile_len,"%s/%s_ea%02d.fil",DATADIR,currt_string,station_id);
}

void get_fbfile (char* fbfile, ssize_t fbfile_len, char* inchdr, vdif_header* vdhdr)
{
  // Open up filterbank file using timestamp and antenna
  int station_id = 0;
  ascii_header_get (inchdr, "STATIONID", "%d", &station_id);
  char currt_string[128];
  time_t epoch_seconds = vdif_to_unixepoch(vdhdr);
  struct tm utc_time;
  gmtime_r (&epoch_seconds, &utc_time);
  strftime (currt_string,sizeof(currt_string), "%Y%m%d_%H%M%S", &utc_time);
  *(currt_string+15) = 0;
  if (CHANMIN < 2411)
    snprintf (fbfile,fbfile_len,"%s/%s_muos_ea%02d.fil",DATADIR,currt_string,station_id);
  else
    snprintf (fbfile,fbfile_len,"%s/%s_ea%02d.fil",DATADIR,currt_string,station_id);
}

void check_buffer (dada_hdu_t* hdu, multilog_t* log)
{
  ipcbuf_t* buf = (ipcbuf_t*) hdu->data_block;
  uint64_t m_nbufs = ipcbuf_get_nbufs (buf);
  uint64_t m_full_bufs = ipcbuf_get_nfull (buf);
  if (m_full_bufs == (m_nbufs - 1))
  {
    fprintf (stderr,"failed buffer check\n");
    fflush (stderr);
    dadacheck (dada_hdu_unlock_write (hdu));
    multilog (mlog, LOG_ERR,
        "Only one free buffer left!  Aborting output.\n");
    exit (EXIT_FAILURE);
  }
}

void check_ipcio_write (dada_hdu_t* hdu, char* buf, size_t to_write, multilog_t* log)
{
  size_t written = ipcio_write (hdu->data_block,buf,to_write);
  if (written != to_write)
  {
    fprintf (stderr, "failed ipcio write\n"); 
    fflush (stderr);
    multilog (mlog, LOG_ERR, "Tried to write %lu bytes to psrdada buffer but only wrote %lu.", to_write, written);
    exit (EXIT_FAILURE);
  }
}

int main (int argc, char *argv[])
{
  // register SIGINT handling
  signal (SIGINT, sigint_handler);

  // register exit function
  atexit (exit_handler);

  int exit_status = EXIT_SUCCESS;
  key_t key_in = 0x40;
  key_t key_out = 0x0;
  key_t key_co = 0x0;
  //uint64_t port = READER_SERVICE_PORT;
  int stdout_output = 0;
  int write_fb = 2;
  int npol = 1;
  int do_inject_frb = 0;
  int RFI_MODE=2;
  int profile_pass = 0;  
  int single_pass = 0;  
  int gpu_id = 0;
  size_t maxn = 0;

  int arg = 0;
  while ((arg = getopt(argc, argv, "hik:K:C:omw:b:P:r:stg:")) != -1) {

    switch (arg) {

    case 'h':
      usage ();
      return 0;
      
    case 'k':
      if (sscanf (optarg, "%x", &key_in) != 1) {
        fprintf (stderr, "writer: could not parse key from %s\n", optarg);
        return -1;
      }
      break;
    case 'C':
      if (sscanf (optarg, "%x", &key_co) != 1) {
        fprintf (stderr, "writer: could not parse key from %s\n", optarg);
        return -1;
      }
      break;

    case 'K':
      if (sscanf (optarg, "%x", &key_out) != 1) {
        fprintf (stderr, "writer: could not parse key from %s\n", optarg);
        return -1;
      }
      break;

      /*
    case 'p':
      if (sscanf (optarg, "%zu", &port) != 1) {
        fprintf (stderr, "writer: could not parse port from %s\n", optarg);
        return -1;
      }
      break;
      */

    case 'o':
      stdout_output = 1;
      break;

    case 'i':
      do_inject_frb = 1;
      break;

    /*
    case 'm':
      CHANMIN = 1;
      break;
    */

    case 'w':
      if (sscanf (optarg, "%d", &write_fb) != 1) {
        fprintf (stderr, "writer: could not parse write mode %s\n", optarg);
        return -1;
      }
      break;

    case 'b':
      if (sscanf (optarg, "%d", &NBIT) != 1) {
        fprintf (stderr, "writer: could not parse number of bits %s\n", optarg);
        return -1;
      }
      if (!(NBIT==2 || NBIT==4 || NBIT==8)) {
        fprintf (stderr, "Unsupported NBIT!\n");
        return -1;
      }
      break;

    case 'r':
      if (sscanf (optarg, "%d", &RFI_MODE) != 1) {
        fprintf (stderr, "writer: could not parse RFI mode %s\n", optarg);
        return -1;
      }
      if ( (RFI_MODE < 0) || (RFI_MODE > 2) )
      {
        fprintf (stderr, "Unsupported RFI mode!\n");
        return -1;
      }
      break;

    case 'P':
      if (sscanf (optarg, "%d", &npol) != 1) {
        fprintf (stderr, "writer: could not parse number of pols %s\n", optarg);
        return -1;
      }
      if (!(npol==1 || npol==2)) {
        fprintf (stderr, "Unsupported npol!\n");
        return -1;
      }
      break;

    case 's':
      single_pass = 1;
      break;

    case 't':
      profile_pass = 1;
      break;

    case 'g':
      if (sscanf (optarg, "%d", &gpu_id) != 1) {
        fprintf (stderr, "writer: could not parse GPU id %s\n", optarg);
        return -1;
      }
      if (!(gpu_id==0 || gpu_id==1)) {
        fprintf (stderr, "Unsupported GPU id!\n");
        return -1;
      }
      break;

    }
  }

  cudacheck (cudaSetDevice (gpu_id));
  printf ("Setting CUDA device to %d.\n",gpu_id);
  int nsms;
  cudaDeviceGetAttribute (&nsms,cudaDevAttrMultiProcessorCount,gpu_id);

  struct timespec ts_1ms = get_ms_ts (1);
  struct timespec ts_10s = get_ms_ts (10000);

  #if RT_PROFILE
  cudaEvent_t rt_start, rt_stop; //,start_total,stop_total;
  cudaEventCreate (&rt_start);
  cudaEventCreate (&rt_stop);

  float measured_time=0;
  float read_time=0, proc_time=0, write_time=0, rt_elapsed=0;
  #endif

  #if PROFILE
  // support for measuring run times of parts
  cudaEvent_t start,stop;
  cudaEventCreate (&start);
  cudaEventCreate (&stop);
  float alloc_time=0, hdr_time=0, read_time=0, todev_time=0, 
      convert_time=0, kurtosis_time=0, fft_time=0, histo_time=0,
      normalize_time=0, tscrunch_time=0, pscrunch_time=0, 
      digitize_time=0, write_time=0, flush_time=0, misc_time=0,
      elapsed=0;
  #endif
  // measure full run time time
  cudaEvent_t obs_start, obs_stop;
  cudaEventCreate (&obs_start);
  cudaEventCreate (&obs_stop);

  multilog_t* mlog = multilog_open ("process_baseband",0);
  char logfile[128];
  time_t currt = time (NULL);
  struct tm tmpt; 
  gmtime_r (&currt, &tmpt);
  char currt_string[32];
  strftime (currt_string,sizeof(currt_string), "%Y%m%d_%H%M%S", &tmpt);
  currt_string[15] = 0;
  char hostname[MAXHOSTNAME];
  gethostname (hostname,MAXHOSTNAME);
  pid_t pid = getpid();
  snprintf (logfile,128,
      "%s/%s_%s_process_%06d.log",LOGDIR,currt_string,hostname,pid);
  logfile_fp = myopen (logfile, "w");
  multilog_add (mlog, logfile_fp);
  if (stdout_output)
    multilog_add (mlog, stdout);

  // log configuration parameters
  char incoming_hdr[4096]; // borrow this buffer
  strncpy (incoming_hdr, argv[0], sizeof(incoming_hdr)-1);
  for (int i = 1; i < argc; ++i) {
    strcat (incoming_hdr, " ");
    strncat (incoming_hdr, argv[i], 4095-strlen(incoming_hdr));
  }
  multilog (mlog, LOG_INFO, "[PROCESS_BASEBAND] invoked with: \n%s\n", incoming_hdr);
  if (key_out)
    multilog (mlog, LOG_INFO, "Will write to psrdada buffer.\n");

  // sanity checks on configuration parameters
  if (NKURTO != 250 && NKURTO != 500) {
    multilog (mlog, LOG_ERR, "Only NKURTO==250 or 500 supported.\n");
    exit (EXIT_FAILURE);
  }

  // connect to input buffer
  hdu_in = dada_hdu_create (mlog);
  dada_hdu_set_key (hdu_in,key_in);
  if (dada_hdu_connect (hdu_in) != 0) {
    multilog (mlog, LOG_ERR, 
        "Unable to connect to incoming PSRDADA buffer!\n");
    exit (EXIT_FAILURE);
  }

  // connect to output buffer (optional)
  if (key_out)
  {
      hdu_out = dada_hdu_create (mlog);
      dada_hdu_set_key (hdu_out,key_out);
      if (dada_hdu_connect (hdu_out) != 0) {
        multilog (mlog, LOG_ERR, 
            "Unable to connect to outgoing PSRDADA buffer!\n");
        exit (EXIT_FAILURE);
      }
  }
  if (key_co)
  {
    hdu_co = dada_hdu_create (mlog);
    dada_hdu_set_key(hdu_co, key_co);
    if (dada_hdu_connect (hdu_co) != 0) {
     multilog (mlog, LOG_ERR, 
       "Unable to connect to Coadding PSRDADA buffer key=%x!\n",key_co);
     exit (EXIT_FAILURE);
    }
  }

  #if PROFILE
  cudaEventRecord(start,0);
  #endif
  
  // allocate memory for 1 second; in the future, consider keeping
  // multiple buffers to allow for out-of-order data; this implementation
  // will handle dropped data
  unsigned char* udat_hst; cudacheck (
  cudaMallocHost ((void**)&udat_hst, 2*VLITE_RATE) );
  // device memory for 8-bit samples
  unsigned char* udat_dev; cudacheck (
  cudaMalloc ((void**)&udat_dev,2*VLITE_RATE/SEG_PER_SEC) );

  #if DOHISTO
  // memory for sample histograms
  unsigned int* histo_dev; cudacheck (
  cudaMalloc ((void**)&histo_dev,2*256*sizeof(unsigned int)) );
  unsigned int* histo_hst; cudacheck (
  cudaMallocHost ((void**)&histo_hst,2*256*sizeof(unsigned int)) );
  #endif

  // voltage samples in a chunk, 2*VLITE_RATE/SEG_PER_SEC (both pols)
  size_t samps_per_chunk = 2*VLITE_RATE/SEG_PER_SEC;
  // FFTs per processing chunk (2 per pol)
  int fft_per_chunk = 2*FFTS_PER_SEG;

  cufftHandle plan;
  cufftcheck (cufftPlan1d (&plan,NFFT,CUFFT_R2C,fft_per_chunk));

  // memory for FFTs
  cufftReal* fft_in; cudacheck (
  cudaMalloc ((void**)&fft_in,sizeof(cufftReal)*samps_per_chunk) );
  cufftComplex* fft_out; cudacheck (
  cudaMalloc ((void**)&fft_out,sizeof(cufftComplex)*fft_per_chunk*NCHAN) );
  // if RFI_MODE==1, just use the same buffers; otherwise, duplicate
  cufftReal* fft_in_kur = fft_in;
  if (2 == RFI_MODE)
    cudacheck ( cudaMalloc (
      (void**)&fft_in_kur, sizeof(cufftReal)*samps_per_chunk) );
  cufftComplex* fft_out_kur = fft_out;
  if (2 == RFI_MODE)
    cudacheck ( cudaMalloc (
      (void**)&fft_out_kur, sizeof(cufftComplex)*fft_per_chunk*NCHAN) );

  // device memory for kurtosis statistics, uses NKURTO samples, both pols
  // NKURTO must be commensurate with samps_per_chunk/2, i.e. samples per
  // chunk in a pol
  size_t nkurto_per_chunk = samps_per_chunk / NKURTO;

  // extra factor of 2 to store both power and kurtosis statistics
  // storage for high time resolution and filterbank block ("fb") scales
  cufftReal *pow_dev(NULL), *kur_dev(NULL), 
            *pow_fb_dev(NULL), *kur_fb_dev(NULL);
  if (RFI_MODE) {
    cudacheck (
      cudaMalloc ((void**)&pow_dev,2*sizeof(cufftReal)*nkurto_per_chunk) );
    cudacheck (
      cudaMalloc ((void**)&pow_fb_dev,2*sizeof(cufftReal)*fft_per_chunk) );
    kur_dev = pow_dev + nkurto_per_chunk;
    kur_fb_dev = pow_fb_dev + fft_per_chunk;
  }

  // store D'Agostino statistic for thresholding
  // only using one per pol now, but keep memory size for both pols
  // to make life easier; the values are duplicated
  cufftReal* dag_dev=NULL; if (RFI_MODE) cudacheck (
  cudaMalloc ((void**)&dag_dev,sizeof(cufftReal)*nkurto_per_chunk) );
  cufftReal* dag_fb_dev=NULL; if (RFI_MODE) cudacheck (
  cudaMalloc ((void**)&dag_fb_dev,sizeof(cufftReal)*fft_per_chunk) );

  // store a set of to re-normalize voltages after applying kurtosis
  cufftReal* kur_weights_dev=NULL; if (RFI_MODE) cudacheck (
  cudaMalloc ((void**)&kur_weights_dev,sizeof(cufftReal)*fft_per_chunk) );

  #if WRITE_KURTO
  // storage on host if writing out kurtosis statistics
  cufftReal* kur_hst; cudacheck (
  cudaMallocHost ((void**)&kur_hst,2*sizeof(cufftReal)*nkurto_per_chunk));
  cufftReal* kur_fb_hst; cudacheck (
  cudaMallocHost ((void**)&kur_fb_hst,2*sizeof(cufftReal)*fft_per_chunk) );
  cufftReal* kur_weights_hst; cudacheck (
  cudaMallocHost ((void**)&kur_weights_hst,sizeof(cufftReal)*fft_per_chunk) );
  #endif


  // NB, reduce by a further factor of 2 if pscrunching
  int polfac = npol==1?2:1;
  int scrunch = (fft_per_chunk*NCHAN)/(polfac*NSCRUNCH);
  cufftReal* fft_ave; cudacheck (
  cudaMalloc ((void**)&fft_ave,sizeof(cufftReal)*scrunch) );
  cufftReal* fft_ave_kur=fft_ave;
  if (2 == RFI_MODE)
    cudacheck ( cudaMalloc (
        (void**)&fft_ave_kur,sizeof(cufftReal)*scrunch) );

  // error check that NBIT is commensurate with trimmed array size
  int trim = (fft_per_chunk*(CHANMAX-CHANMIN+1))/(polfac*NSCRUNCH);
  if (trim % (8/NBIT) != 0) {
    multilog (mlog, LOG_ERR, 
        "Selected channel and bit scheme is not commensurate!.\n");
    exit (EXIT_FAILURE);
  }

  // reduce array size by packing of samples into byte
  trim /= (8/NBIT);
  unsigned char* fft_trim_u_dev; cudacheck (
  cudaMalloc ((void**)&fft_trim_u_dev,trim) );
  unsigned char* fft_trim_u_hst;
  // we only need to malloc this if we need to buffers; otherwise will
  // use the big internal buffer
  if (2 == RFI_MODE) cudacheck (
    cudaMallocHost ((void**)&fft_trim_u_hst,trim) );

  unsigned char* fft_trim_u_kur_dev = fft_trim_u_dev;
  if (2 == RFI_MODE)
    cudacheck (cudaMalloc ((void**)&fft_trim_u_kur_dev,trim) );
  unsigned char* fft_trim_u_kur_hst = NULL;
  //if (2 == RFI_MODE)
    //cudacheck (cudaMallocHost ((void**)&fft_trim_u_kur_hst,trim) );

  // memory for a 10s buffer of output filterbank data
  int output_buf_sec = 10;
  int output_buf_seg_size = trim;
  int output_buf_size = output_buf_seg_size*output_buf_sec*SEG_PER_SEC;
  unsigned char* output_buf_mem;
  cudacheck (cudaMallocHost ((void**)&output_buf_mem,output_buf_size) );
  unsigned char* output_buf_cur = output_buf_mem;

  // memory for running bandpass correction; 2 pol * NCHAN
  cufftReal* bp_dev; cudacheck (
  cudaMalloc ((void**)&bp_dev,sizeof(cufftReal)*NCHAN*2));
  cudacheck (cudaMemset (bp_dev, 0, sizeof(cufftReal)*NCHAN*2));
  cufftReal* bp_kur_dev = bp_dev;
  if (2 == RFI_MODE )
  {
    cudacheck (cudaMalloc (
        (void**)&bp_kur_dev,sizeof(cufftReal)*NCHAN*2));
    cudacheck (cudaMemset (bp_kur_dev, 0, sizeof(cufftReal)*NCHAN*2));
  }

  // memory for FRB injection
  float* frb_delays_dev = NULL;
  if (do_inject_frb)
  {
    cudacheck (cudaMalloc ((void**)&frb_delays_dev, sizeof(float)*NCHAN));
    // DM of 80 fits almost exactly in 1s of data across the whole band
    set_frb_delays <<< NCHAN/NTHREAD+1, NTHREAD >>> (frb_delays_dev, 80);
    cudacheck (cudaGetLastError () );

    /*
    float* frb_delays_hst = NULL;
    cudacheck (cudaMallocHost ((void**)&frb_delays_hst, sizeof(float)*NCHAN));
    cudacheck (cudaMemcpy (

          frb_delays_hst,frb_delays_dev,sizeof(float)*NCHAN,
          cudaMemcpyDeviceToHost) );
    for (int ichan=0; ichan < 6250; ichan += 100)
      fprintf (stdout, "frb_delay %d = %.6f\n", ichan, frb_delays_hst[ichan]);
    */
  }


  #if PROFILE
  CUDA_PROFILE_STOP(start,stop,&alloc_time)
  #endif

  // constants for bandpass normalization: tsamp/tsmooth, giving a time
  // constant of about tsmooth secs.
  double tsmooth = 1;
  double tsamp = double(NFFT)/VLITE_RATE*NSCRUNCH;
  float bp_scale = tsamp/tsmooth;

  // allocate a frame buffer on the stack
  char frame_buff[5032];
  vdif_header* vdhdr = (vdif_header*) frame_buff;
  uint8_t* dat = (uint8_t*)(frame_buff + 32);

  /*
  // connect to control socket
  Connection conn;
  conn.sockoptval = 1; //release port immediately after closing connection
  if (port) {
    if (serve (port, &conn) < 0) {
      multilog (mlog, LOG_ERR,
          "Failed to create control socket on port %d.\n", port);
      exit (EXIT_FAILURE);
    }
    fcntl (conn.rqst, F_SETFL, O_NONBLOCK); // set up for polling
  }
  char cmd_buff[32];
  */

  // connect to multicast control socket
  mc_control_sock = open_mc_socket (mc_vlitegrp, MC_READER_PORT,
      (char*)"Control Socket [Process Baseband]", NULL, mlog);
  char mc_control_buff[32];
  fcntl (mc_control_sock, F_SETFL, O_NONBLOCK); // set up for polling
  int cmds[5] = {0,0,0,0,0};

  int quit = 0;
  int inject_frb_now = 0;

  // TEMP for profiling
  if (profile_pass)
  {
    char cmd[256];
    snprintf (cmd, 255, "/home/vlite-master/mtk/src/readbase /home/vlite-master/mtk/baseband/20150302_151845_B0329+54_ev_0005.out.uw"); 
    multilog (mlog, LOG_INFO, "%s\n", cmd);
    popen (cmd, "r");
    nanosleep (&ts_10s,NULL);
  }

  // This point start to loop over observations.
  while(true)
  {

  if (quit)
    break;

  // TODO -- probably want to make this have a timeout to avoid pegging the CPU
  // TODO -- make writer issue a START over multicast so we don't need to block
  //         except, we don't want all writers talking to all PBs!
  get_cmds (cmds, mc_control_sock, mc_control_buff, 32, logfile_fp);
  if (cmds[2]) // CMD_QUIT
    break;
  //if (!cmds[0]) // CMD_START
  //  continue;

  dadacheck (dada_hdu_lock_read (hdu_in) );

  // this will block until a header has been created by writer, 
  // signalling the start of an observation
  // TODO -- this is encapsulated in dada_hdu_open
  multilog (mlog, LOG_INFO, "Waiting for DADA header.\n");
  fflush (logfile_fp);
  uint64_t hdr_size = 0;
  char* ascii_hdr = ipcbuf_get_next_read (hdu_in->header_block,&hdr_size);
  multilog (mlog, LOG_INFO, "Received header with size %d.\n", hdr_size);

  if (ascii_hdr == NULL) {
    // this could be an error condition, or it could be a result of 
    // shutting down the data acquisition system; check to see if a
    // CMD_QUIT has been issued in order to log the appropriate outcome

    if (test_for_cmd (CMD_QUIT, mc_control_sock, mc_control_buff, 32, logfile_fp))
      multilog (mlog, LOG_INFO, "Received CMD_QUIT, indicating data taking is ceasing.  Exiting.\n");
    else
    {
      // did not receive CMD_QUIT, so this is an abnormal state
      multilog (mlog, LOG_ERR, "PSRDADA read failed unexpectedly.  Terminating.\n");
      exit (EXIT_FAILURE);
    }
    break;
  }

  cudaEventRecord(obs_start,0);

  multilog (mlog, LOG_INFO, "psrdada header:\n%s",ascii_hdr);
  multilog (mlog, LOG_INFO, "Beginning new observation.\n\n");
  fflush (logfile_fp);
  memcpy (incoming_hdr,ascii_hdr,hdr_size);
  dadacheck (ipcbuf_mark_cleared (hdu_in->header_block));

  // July 17, 2019.  Changes in writer ensure all arriving data are
  // contiguous and aligned to 1-s boundaries in the buffer. Read first
  // frame and verify this, and use for metadata.
  ssize_t nread = ipcio_read (hdu_in->data_block,frame_buff,VD_FRM);
  if (nread != VD_FRM)
  {
    multilog (mlog, LOG_ERR, "Problem reading first bloody frame!  Bailing.");
    exit (EXIT_FAILURE);
  }
  int thread = getVDIFThreadID (vdhdr) != 0;
  int frame = getVDIFFrameNumber (vdhdr);
  if ((frame != 0) && (thread != 0))
  {
    multilog (mlog, LOG_ERR, "Incoming data were not aligned!");
    exit (EXIT_FAILURE);
  }

  // Open up filterbank file at appropriate time
  char fbfile[256], cofbfile[256];
  get_fbfile (fbfile, 256, incoming_hdr, vdhdr);
  get_cofbfile (cofbfile, 256, incoming_hdr, vdhdr);
  char fbfile_kur[256], cofbfile_kur[256];
  change_extension (fbfile, fbfile_kur, ".fil", "_kur.fil");
  change_extension (cofbfile, cofbfile_kur, ".fil", "_kur.fil");

  // use same scheme for ancillary data files
  #if DOHISTO
  char histofile[256] = "";
  change_extension (fbfile, histofile, ".fil", ".histo");
  #endif
  #if WRITE_KURTO
  char kurtofile[256] = "";
  change_extension (fbfile, kurtofile, ".fil", ".kurto");
  char block_kurtofile[256] = "";
  change_extension (fbfile, block_kurtofile, ".fil", ".block_kurto");
  char weightfile[256] = "";
  change_extension (fbfile, weightfile, ".fil", ".weights");
  #endif

  // cache file name in case we need to change it to /dev/null
  char heimdall_file[256] = "";
  char coheimdall_file[256] = "";
  strncpy (heimdall_file,RFI_MODE?fbfile_kur:fbfile,255);
  strncpy (coheimdall_file,RFI_MODE?cofbfile_kur:cofbfile,255);

  // check for write to /dev/null and udpate file names if necessary
  int write_to_null = 0==write_fb;
  if (write_fb == 1)
  {
    char name[OBSERVATION_NAME_SIZE];
    ascii_header_get (incoming_hdr, "NAME", "%s", name);
    name[OBSERVATION_NAME_SIZE-1] = '\0';
    double ra = 0;
    ascii_header_get (incoming_hdr, "RA", "%lf", &ra);
    double dec = 0;
    ascii_header_get (incoming_hdr, "DEC", "%lf", &dec);
	  char datasetId[EXECUTOR_DATASETID_SIZE];
    ascii_header_get (incoming_hdr, "DATAID", "%s", &datasetId);

    if (check_coords (ra, dec) )
    {
      multilog (mlog, LOG_INFO,
          "Source %s matches target coords, recording filterbank data.\n",
          name);
      send_email (name, hostname);
    }
    else if (check_name (name))
    {
      multilog (mlog, LOG_INFO,
          "Source %s matches target list, recording filterbank data.\n",
          name);
      send_email (name, hostname);
    }
    else if (check_id (datasetId))
    {
      multilog (mlog, LOG_INFO,
          "DATAID %s matches list, recording filterbank data.\n",
          datasetId);
      send_email (datasetId, hostname);
    }
    else {
      write_to_null = 1;
      multilog (mlog, LOG_INFO,
          "Source %s not on target list, disabling filterbank data.\n",
          name);
    }
  }
  if (write_to_null)
  {
    multilog (mlog, LOG_INFO,
        "Filterbank output disabled.  Would have written to %s.\n",fbfile);
    strncpy (fbfile,"/dev/null",255); 
    #if DOHISTO
    strncpy (histofile,"/dev/null",255); 
    #endif
    if (RFI_MODE)
      strncpy (fbfile_kur,"/dev/null",255); 
    #if WRITE_KURT0
      strncpy (kurtofile,"/dev/null",255); 
      strncpy (block_kurtofile,"/dev/null",255); 
      strncpy (weightfile,"/dev/null",255); 
    #endif
  }

  //FILE *fb_fp,*fb_kurto_fp=NULL;
  uint64_t fb_bytes_written = 0;
  // this buffer size is correct for 100ms of default data
  if (RFI_MODE == 0 || RFI_MODE == 2 )
  {
    fb_fp = myopen (fbfile, "wb", true, trim);
    multilog (mlog, LOG_INFO,
        "Writing no-RFI-excision filterbanks to %s.\n",fbfile);
  }
  else
  {
    fb_fp = myopen (fbfile_kur, "wb", trim);
    multilog (mlog, LOG_INFO,
        "Writing RFI-excision filterbanks to %s.\n",fbfile_kur);
  }
  if (RFI_MODE == 2)
  {
    fb_kurto_fp = myopen (fbfile_kur, "wb", trim);
    multilog (mlog, LOG_INFO,
        "Writing RFI-excision filterbanks to %s.\n",fbfile_kur);
  }

  #if DOHISTO
  FILE *histo_fp = myopen (histofile, "wb", true);
  multilog (mlog, LOG_INFO, "Writing histograms to %s.\n",histofile);
  #endif

  #if WRITE_KURTO
  FILE *kurto_fp = myopen (kurtofile, "wb", true);
  multilog (mlog, LOG_INFO, "Writing kurtosis to %s.\n",kurtofile);
  FILE *block_kurto_fp = myopen (block_kurtofile, "wb", true);
  multilog (
       mlog, LOG_INFO, "Writing block kurtosis to %s.\n",block_kurtofile);
  FILE *weight_fp = myopen (weightfile, "wb", true);
  multilog (mlog, LOG_INFO, "Writing weights to %s.\n",weightfile);
  #endif
  fflush (logfile_fp);

  #if PROFILE
  cudaEventRecord(start,0);
  #endif

  // write out psrdada header; NB this locks the hdu for writing
  if (key_out)
  {
    write_psrdada_header (hdu_out, incoming_hdr, vdhdr, npol, heimdall_file);
    fprintf (stderr, "write psrdada header\n");
  }
  if (key_co) {
    write_psrdada_header(hdu_co, incoming_hdr, vdhdr, npol, coheimdall_file); 
    fprintf(stderr, "Write Coadd header\n");
  }

  // write out a sigproc header
  int station_id = write_sigproc_header (fb_fp, incoming_hdr, vdhdr, npol);
  if (2 == RFI_MODE)
    write_sigproc_header (fb_kurto_fp, incoming_hdr, vdhdr, npol);

  #if PROFILE
  CUDA_PROFILE_STOP(start,stop,&hdr_time)
  #endif

  // set up sample# trackers and copy first frame into 1s-buffer
  int current_sec = getVDIFFrameSecond(vdhdr);
  int current_thread = getVDIFThreadID (vdhdr) != 0;
  multilog (mlog, LOG_INFO, "Starting sec=%d, thread=%d\n",current_sec,current_thread);

  double integrated = 0.0;
  int integrated_sec = 0;
  output_buf_cur = output_buf_mem;

  #if RT_PROFILE
  cudaEventRecord(rt_start,0);
  #endif 
  measured_time = 0;


  while(true) // loop over data packets
  {
    int thread = getVDIFThreadID (vdhdr) != 0;
    int sample = getVDIFFrameNumber (vdhdr);
    int second = getVDIFFrameSecond (vdhdr);

    // this segment will continue loop until one full second is read
    if (second==current_sec)
    {

      #if PROFILE
      cudaEventRecord (start,0);
      #endif 

      //#if RT_PROFILE
      //cudaEventRecord(rt_start,0);
      //#endif 

      // copy previous frame to contiguous 1-s buffer
      size_t idx = VLITE_RATE*thread + sample*VD_DAT;
      memcpy (udat_hst + idx,dat,VD_DAT);

      // read the next frame into the buffer
      ssize_t nread = ipcio_read (hdu_in->data_block,frame_buff,VD_FRM);

      #if PROFILE
      CUDA_PROFILE_STOP(start,stop,&elapsed)
      read_time += elapsed;
      #endif

      //#if RT_PROFILE
      //CUDA_PROFILE_STOP(rt_start,rt_stop,&rt_elapsed)
      //read_time += rt_elapsed;
      //#endif

      if (nread < 0 ) {
        multilog (mlog, LOG_ERR, "Error on nread=%d.\n",nread);
        exit (EXIT_FAILURE);
      }

      // THIS IS THE PRIMARY SUCCESSFUL EXIT FROM THE PACKET LOOP HERE:
      // nread==0 --> break

      if (nread != VD_FRM) {
        if (nread != 0)
          multilog (mlog, LOG_INFO,
            "Packet size=%ld, expected %ld.  Aborting this observation.\n",
            nread,VD_FRM);
        break; // potential end of observation, or error
      }

      continue; // back to top of loop over packets
    }

    if (second - current_sec > 1)
    {
      // this should only happen when observing over epoch changes,
      // leap seconds, or when there are major data drops.  When that
      // happens, just restart the code.
      multilog (mlog, LOG_ERR, "Major data skip!  (%d vs. %d; thread = %d) Aborting this observation.\n", second,current_sec,thread);
      exit_status = EXIT_FAILURE;
      quit = 1;
      break;
    }

    // check for a QUIT command once every second
    if (test_for_cmd (CMD_QUIT, mc_control_sock, mc_control_buff, 32, logfile_fp))
    {
      multilog (mlog, LOG_INFO, "Received CMD_QUIT, indicating data taking is ceasing.  Exiting.\n");
      quit = 1;
      break;
    }
    if (quit) {
      multilog (mlog, LOG_INFO, "Received CMD_QUIT.  Exiting!.\n");
      // this will break out of data loop, write out file, then break out
      // of observation loop to exit program
      break;
    }

    // keep log on disk up to date
    fflush (logfile_fp);

    // check for FRB injection conditions
    inject_frb_now = do_inject_frb && (current_sec%60==0);
    if (inject_frb_now)
      multilog (mlog, LOG_INFO,
        "Injecting an FRB with integrated = %.2f!!!.\n", integrated);

    // check that both buffers are full
    current_sec = second;

    // do dispatch -- break into chunks to fit in GPU memory; this is
    // currently 100 milliseconds
    for (int iseg = 0; iseg < SEG_PER_SEC; iseg++)
    {

      #if PROFILE
      cudaEventRecord(start,0);
      #endif 


      // need to copy pols separately
      for (int ipol=0; ipol < 2; ++ipol)
      {
        unsigned char* dst = udat_dev + ipol*samps_per_chunk/2;
        unsigned char* src = udat_hst + ipol*VLITE_RATE + iseg*samps_per_chunk/2;
        cudacheck (cudaMemcpy (dst,src,samps_per_chunk/2,cudaMemcpyHostToDevice) );
      }
      #if PROFILE
      CUDA_PROFILE_STOP(start,stop,&elapsed)
      todev_time += elapsed;
      #endif

      ////// HISTOGRAM //////
      #if DOHISTO
      #if PROFILE
      cudaEventRecord (start,0);
      #endif 
      // compute sample histogram
      // memset implicitly barriers between blocks, necessary
      cudacheck (cudaMemset (histo_dev,0,2*256*sizeof(unsigned int)));
      histogram <<<nsms*32,512>>> (udat_dev,histo_dev,samps_per_chunk);
      cudacheck (cudaGetLastError () );
      #if PROFILE
      CUDA_PROFILE_STOP(start,stop,&elapsed)
      histo_time += elapsed;
      #endif
      #endif

      ////// RT_PROFILE
      //#if RT_PROFILE
      //cudaEventRecord (rt_start,0);
      //#endif
      ////// CONVERT UINTS to FLOATS //////
      #if PROFILE
      cudaEventRecord (start,0);
      #endif 
      convertarray <<<nsms*32,NTHREAD>>> (fft_in,udat_dev,samps_per_chunk);
      cudacheck (cudaGetLastError () );
      #if PROFILE
      CUDA_PROFILE_STOP(start,stop,&elapsed)
      convert_time += elapsed;
      #endif

      ////// CALCULATE KURTOSIS STATISTICS //////
      if (RFI_MODE)
      {
        #if PROFILE
        cudaEventRecord (start,0);
        #endif

        // calculate high time resolution kurtosis (250 or 500 samples)
        kurtosis <<<nkurto_per_chunk, 256>>> (
            fft_in,pow_dev,kur_dev);
        cudacheck (cudaGetLastError () );

        // compute the thresholding statistic
        // NB now modified to combine polarizations
        compute_dagostino <<<nsms*32,NTHREAD>>> (
            kur_dev,dag_dev,nkurto_per_chunk/2);
        cudacheck (cudaGetLastError () );

        // calculate coarser kurtosis (for entire filterbank sample, e.g. 12500 samples)
        // NB this relies on results of previous D'Agostino calculation
        block_kurtosis <<<fft_per_chunk/8,256>>> (
            pow_dev,kur_dev,dag_dev,pow_fb_dev,kur_fb_dev);
        cudacheck (cudaGetLastError () );
        // NB now modified to combine polarizations
        compute_dagostino2 <<<nsms*32,NTHREAD>>> (
            kur_fb_dev,dag_fb_dev,fft_per_chunk/2);
        cudacheck (cudaGetLastError () );

        cudacheck (cudaMemset (
            kur_weights_dev,0,sizeof(cufftReal)*fft_per_chunk) );
        // (1) NB that fft_in_kur==fft_in if not writing both streams
        // (2) original implementation had a block for each pol; keeping
        //     that, but two blocks now access the same Dagostino entry
        //if (integrated >= 0.1)
        if (1)
        {
          apply_kurtosis <<<nkurto_per_chunk, 256>>> (
              fft_in,fft_in_kur,dag_dev,dag_fb_dev,kur_weights_dev);
          cudacheck (cudaGetLastError () );
        }
        else
        {
          apply_kurtosis_fake <<<nkurto_per_chunk, 256>>> (
              fft_in,fft_in_kur,dag_dev,dag_fb_dev,kur_weights_dev);
          cudacheck (cudaGetLastError () );
        }

        #if WRITE_KURTO
        cudacheck (cudaMemcpy (
            kur_weights_hst, kur_weights_dev, 
            sizeof(cufftReal)*fft_per_chunk, cudaMemcpyDeviceToHost) );
        #endif

        #if PROFILE
        CUDA_PROFILE_STOP(start,stop,&elapsed)
        kurtosis_time += elapsed;
        #endif
      }

      ////// PERFORM FFTs //////
      #if PROFILE
      cudaEventRecord (start,0);
      #endif 
      cufftcheck (cufftExecR2C (plan,fft_in,fft_out) );
      if (RFI_MODE > 1)
        cufftcheck (cufftExecR2C (plan,fft_in_kur,fft_out_kur) );
      #if PROFILE
      CUDA_PROFILE_STOP(start,stop,&elapsed)
      fft_time += elapsed;
      #endif

      ////// INJECT FRB AS REQUESTED //////
      if (inject_frb_now > 0)
      {
        // NB that inject_frb_now is only reset every 1s, so we also use
        // it to keep track of how many segments have elapsed since the
        // FRB time, since this loop is over 100ms chunks which will be
        // < dispersed FRB width, typically

        float frb_width = 2e-3*SEG_PER_SEC*FFTS_PER_SEG; //2ms 
        float frb_amp = 1.05; // gives a single-antenna S/N of about 25-30 for 2ms
        int nfft_since_frb = (inject_frb_now-1)*FFTS_PER_SEG;
        inject_frb <<< NCHAN/NTHREAD+1,NTHREAD >>> (fft_out, 
            frb_delays_dev, nfft_since_frb, frb_width, frb_amp);
        cudacheck ( cudaGetLastError () );
        if (RFI_MODE > 1)
        {
          inject_frb <<< NCHAN/NTHREAD+1,NTHREAD >>> (fft_out_kur, 
              frb_delays_dev, nfft_since_frb, frb_width, frb_amp);
          cudacheck ( cudaGetLastError () );
        }
        inject_frb_now += 1;
      }

      ////// NORMALIZE BANDPASS //////
      #if PROFILE
      cudaEventRecord(start,0);
      #endif 
      if (RFI_MODE == 0 || RFI_MODE == 2)
      {
        detect_and_normalize2 <<<(NCHAN*2)/NTHREAD+1,NTHREAD>>> (
            fft_out,bp_dev,bp_scale);
        cudacheck ( cudaGetLastError () );
      }
      if (RFI_MODE == 1 || RFI_MODE == 2)
      {
        detect_and_normalize3 <<<(NCHAN*2)/NTHREAD+1,NTHREAD>>> (
            fft_out_kur,kur_weights_dev,bp_kur_dev,bp_scale);
        cudacheck ( cudaGetLastError () );
      }
      #if PROFILE
      CUDA_PROFILE_STOP(start,stop,&elapsed)
      normalize_time += elapsed;
      #endif

      ////// [OPTIONALLY] ADD POLARIZATIONS //////
      #if PROFILE
      cudaEventRecord (start,0);
      #endif 
      maxn = (fft_per_chunk*NCHAN)/polfac;
      if (npol==1) {
        if (RFI_MODE == 0 || RFI_MODE == 2)
        {
          pscrunch <<<nsms*32,NTHREAD>>> (fft_out,maxn);
          cudacheck ( cudaGetLastError () );
        }
        if (RFI_MODE == 1 || RFI_MODE == 2)
        {
          pscrunch_weights <<<nsms*32,NTHREAD>>> (
              fft_out_kur,kur_weights_dev,maxn);
          cudacheck ( cudaGetLastError () );
        }
      }
      #if PROFILE
      CUDA_PROFILE_STOP(start,stop,&elapsed)
      pscrunch_time += elapsed;
      #endif

      ////// AVERAGE TIME DOMAIN //////
      #if PROFILE
      cudaEventRecord( start,0);
      #endif 
      maxn /= NSCRUNCH;
      if (RFI_MODE == 0 || RFI_MODE == 2)
      {
        tscrunch <<<nsms*32,NTHREAD>>> (fft_out,fft_ave,maxn);
        cudacheck ( cudaGetLastError () );
      }
      if (RFI_MODE == 1 || RFI_MODE == 2)
      {
        tscrunch_weights <<<nsms*32,NTHREAD>>> (
            fft_out_kur,fft_ave_kur,kur_weights_dev,maxn);
        cudacheck ( cudaGetLastError () );
      }
      #if PROFILE
      CUDA_PROFILE_STOP(start,stop,&elapsed)
      tscrunch_time += elapsed;
      #endif

      ////// TRIM CHANNELS AND DIGITIZE //////
      #if PROFILE
      cudaEventRecord (start,0);
      #endif 
      maxn = (CHANMAX-CHANMIN+1)*(maxn/NCHAN)/(8/NBIT);
      switch (NBIT)
      {
        case 2:
          sel_and_dig_2b <<<nsms*32,NTHREAD>>> (
              fft_ave,fft_trim_u_dev,maxn, npol);
          if (RFI_MODE > 1)
            sel_and_dig_2b <<<nsms*32,NTHREAD>>> (
                fft_ave_kur,fft_trim_u_kur_dev,maxn, npol);
          break;
        case 4:
          sel_and_dig_4b <<<nsms*32,NTHREAD>>> (
              fft_ave,fft_trim_u_dev,maxn, npol);
          if (RFI_MODE > 1)
            sel_and_dig_4b <<<nsms*32,NTHREAD>>> (
                fft_ave_kur,fft_trim_u_kur_dev,maxn, npol);
          break;
        case 8:
          sel_and_dig_8b <<<nsms*32,NTHREAD>>> (
              fft_ave,fft_trim_u_dev,maxn, npol);
          if (RFI_MODE > 1)
          sel_and_dig_8b <<<nsms*32,NTHREAD>>> (
              fft_ave_kur,fft_trim_u_kur_dev,maxn, npol);
          break;
        default:
          sel_and_dig_2b <<<nsms*32,NTHREAD>>> (
              fft_ave,fft_trim_u_dev,maxn, npol);
          if (RFI_MODE > 1)
            sel_and_dig_2b <<<nsms*32,NTHREAD>>> (
                fft_ave_kur,fft_trim_u_kur_dev,maxn, npol);
          break;
      }
      cudacheck ( cudaGetLastError () );
      #if PROFILE
      CUDA_PROFILE_STOP(start,stop,&elapsed)
      digitize_time += elapsed;
      #endif

      #if PROFILE
      cudaEventRecord (start,0);
      #endif 

      // copy filterbanked data back to host; use big buffer to avoid
      // a second copy; NB that if we are only recording a single RFI
      // excision mode, the _kur buffer points to same place.  And if we
      // are recording both, then we output the kurtosis.  So we can
      // just always record the kurtosis to the output buffer, and if
      // we are recording both, copy the second to the small buffer
      fft_trim_u_kur_hst = output_buf_cur;
      cudacheck (cudaMemcpy (
          fft_trim_u_kur_hst,fft_trim_u_kur_dev,maxn,cudaMemcpyDeviceToHost) );
      if (2 == RFI_MODE)
        cudacheck (cudaMemcpy ( fft_trim_u_hst,
            fft_trim_u_dev,maxn,cudaMemcpyDeviceToHost) );
      output_buf_cur += output_buf_seg_size; 

      #if DOHISTO
      // copy histograms
      cudacheck (cudaMemcpy (
            histo_hst,histo_dev,2*256*sizeof(unsigned int),
            cudaMemcpyDeviceToHost) );
      #endif

      #if WRITE_KURTO
      // copy kurtosis
      cudacheck (cudaMemcpy (
            kur_hst,kur_dev,2*nkurto_per_chunk*sizeof(cufftReal),
            cudaMemcpyDeviceToHost) );
      cudacheck (cudaMemcpy (
            kur_fb_hst,kur_fb_dev,2*fft_per_chunk*sizeof(cufftReal),
            cudaMemcpyDeviceToHost) );
      #endif

      #if PROFILE
      CUDA_PROFILE_STOP(start,stop,&elapsed)
      misc_time += elapsed;
      #endif

      //#if RT_PROFILE
      //CUDA_PROFILE_STOP(rt_start,rt_stop,&rt_elapsed)
      //proc_time += rt_elapsed;
      //#endif

      // finally, push the filterbanked time samples onto psrdada buffer
      // and/or write out to sigproc

      #if PROFILE
      cudaEventRecord (start,0);
      #endif 

      //#if RT_PROFILE
      //cudaEventRecord (rt_start,0);
      //#endif

      if (key_co)
      {
        char* outbuff = RFI_MODE==2? (char*)fft_trim_u_kur_hst:
                                     (char*)fft_trim_u_hst;
        check_buffer (hdu_co, mlog);
        check_ipcio_write (hdu_co, outbuff, maxn, mlog);
      }

      //#if RT_PROFILE
      //CUDA_PROFILE_STOP(rt_start,rt_stop,&rt_elapsed)
      //write_time += rt_elapsed;
      //#endif

      // TODO -- tune this I/O.  The buffer size is set to 8192, but
      // according to fstat the nfs wants a block size of 1048576! Each
      // 100ms of data is 65536 with the current parameters.  So optimally
      // we would buffer in memory for the full 1 second before a write.
      // However, a simple improvement will be reducing the write calls by
      // a factor of 8 by either changing the buffer size or using the

      // TODO -- add error checking for these writes

      fwrite (fft_trim_u_hst,1,maxn,fb_fp);
      if (RFI_MODE == 2)
        fwrite (fft_trim_u_kur_hst,1,maxn,fb_kurto_fp);
      fb_bytes_written += maxn;

      #if DOHISTO
      fwrite (histo_hst,sizeof(unsigned int),512,histo_fp);
      #endif
      #if WRITE_KURTO
      fwrite (kur_hst,sizeof(cufftReal),2*nkurto_per_chunk,kurto_fp);
      fwrite (kur_fb_hst,sizeof(cufftReal),2*fft_per_chunk,block_kurto_fp);
      fwrite (kur_weights_hst,sizeof(cufftReal),fft_per_chunk,weight_fp);
      #endif
      #if PROFILE
      CUDA_PROFILE_STOP(start,stop,&elapsed)
      write_time += elapsed;
      #endif

      integrated += 1./double(SEG_PER_SEC);

    } // end loop over data segments

    integrated_sec += 1;
    #if RT_PROFILE
    if (integrated_sec%10==0) {
    #else
    if (integrated_sec%30==0) {
    #endif
      // TMP
      #if RT_PROFILE
      CUDA_PROFILE_STOP(rt_start,rt_stop,&rt_elapsed)
      measured_time += rt_elapsed;
      //fprintf (stderr, "data  sec=%d\n",10);
      //fprintf (stderr, "read  sec=%.2f\n",read_time*1e-3);
      //fprintf (stderr, "proc  sec=%.2f\n",proc_time*1e-3);
      //fprintf (stderr, "write sec=%.2f\n",write_time*1e-3);
      if ((measured_time*1e-3-integrated) > 0.5)
        multilog (mlog, LOG_ERR, "Measured time exceeding integrated time: measured: %.2f, integrated: %.2f\n",measured_time*1e-3,integrated);
      cudaEventRecord (rt_start,0);
      read_time=0;proc_time=0;write_time=0;
      #else
      fprintf (stderr, "integrated sec=%d\n",integrated_sec);
      #endif 
    }
    if (integrated_sec >= output_buf_sec)
    {
      output_buf_cur = output_buf_mem;
      int to_write = output_buf_seg_size*SEG_PER_SEC;
      if (integrated_sec == output_buf_sec)
        to_write = output_buf_size;

      // write out either the full buffer or 1s
      if(key_out) {
        check_buffer (hdu_out, mlog);
        check_ipcio_write (hdu_out,(char*)output_buf_cur, to_write, mlog);
      }
    }

  } // end loop over packets

  if (key_out)
  {
    fprintf (stderr, "process_baseband: before dada_hdu_unlock_write\n");
    dadacheck (dada_hdu_unlock_write (hdu_out));
    fprintf (stderr, "process_baseband: after dada_hdu_unlock_write\n");
    fflush (stderr);
  }

  if (key_co)
  {
    multilog(mlog, LOG_INFO, "process_baseband: coadd before dada_hdu_unlock_write\n");
    dadacheck( dada_hdu_unlock_write(hdu_co));
    multilog(mlog, LOG_INFO, "process_baseband: coadd after dada_hdu_unlock_write\n");
  }

  dadacheck (dada_hdu_unlock_read (hdu_in));

  #if PROFILE
  cudaEventRecord(start,0);
  #endif 

  // close files
  if (fb_fp) {fclose (fb_fp); fb_fp = NULL;}
  if (fb_kurto_fp) {fclose (fb_kurto_fp); fb_kurto_fp = NULL;}
  #if DOHISTO
  if (histo_fp) {fclose (histo_fp); histo_fp = NULL;}
  #endif
  #if WRITE_KURTO
  if (kurto_fp) {fclose (kurto_fp); kurto_fp = NULL;}
  if (block_kurto_fp) {fclose (block_kurto_fp); block_kurto_fp = NULL;}
  if (weight_fp) {fclose (weight_fp); weight_fp = NULL;}
  #endif
  uint64_t samps_written = (fb_bytes_written*(8/NBIT))/(CHANMAX-CHANMIN+1);
  multilog (mlog, LOG_INFO, "Wrote %.2f MB (%.2f s) to %s\n",
      fb_bytes_written*1e-6,samps_written*tsamp,fbfile);

  float obs_time;
  CUDA_PROFILE_STOP (obs_start,obs_stop,&obs_time);
  multilog (mlog, LOG_INFO, "Proc Time...%.3f\n", obs_time*1e-3);

  #if PROFILE
  CUDA_PROFILE_STOP(start,stop,&flush_time)
  float sub_time = hdr_time + read_time + todev_time + 
      histo_time + convert_time + kurtosis_time + fft_time + 
      normalize_time + pscrunch_time + tscrunch_time + digitize_time + 
      write_time + flush_time + misc_time;
  multilog (mlog, LOG_INFO, "Alloc Time..%.3f\n", alloc_time*1e-3);
  multilog (mlog, LOG_INFO, "Read Time...%.3f\n", read_time*1e-3);
  multilog (mlog, LOG_INFO, "Copy To Dev.%.3f\n", todev_time*1e-3);
  multilog (mlog, LOG_INFO, "Histogram...%.3f\n", histo_time*1e-3);
  multilog (mlog, LOG_INFO, "Convert.....%.3f\n", convert_time*1e-3);
  multilog (mlog, LOG_INFO, "Kurtosis....%.3f\n", kurtosis_time*1e-3);
  multilog (mlog, LOG_INFO, "FFT.........%.3f\n", fft_time*1e-3);
  multilog (mlog, LOG_INFO, "Normalize...%.3f\n", normalize_time*1e-3);
  multilog (mlog, LOG_INFO, "Pscrunch....%.3f\n", pscrunch_time*1e-3);
  multilog (mlog, LOG_INFO, "Tscrunch....%.3f\n", tscrunch_time*1e-3);
  multilog (mlog, LOG_INFO, "Digitize....%.3f\n", digitize_time*1e-3);
  multilog (mlog, LOG_INFO, "Write.......%.3f\n", write_time*1e-3);
  multilog (mlog, LOG_INFO, "Flush.......%.3f\n", flush_time*1e-3);
  multilog (mlog, LOG_INFO, "Misc........%.3f\n", misc_time*1e-3);
  multilog (mlog, LOG_INFO, "Sum of subs.%.3f\n", sub_time*1e-3);

  // reset values for next loop
  hdr_time=read_time=todev_time=convert_time=kurtosis_time=fft_time=0;
  histo_time=normalize_time=tscrunch_time=pscrunch_time=digitize_time=0;
  write_time=flush_time=misc_time=elapsed=0;

  #endif

  if (profile_pass || single_pass)
    break;

  } // end loop over observations

  // free memory
  if (udat_hst) cudaFreeHost (udat_hst);
  if (udat_dev) cudaFree (udat_dev);
  if (fft_in) cudaFree (fft_in);
  if (fft_out) cudaFree (fft_out);
  if (fft_ave) cudaFree (fft_ave);
  if (fft_trim_u_dev) cudaFree (fft_trim_u_dev);
  if (fft_trim_u_hst) cudaFreeHost (fft_trim_u_hst);
  if (output_buf_mem ) cudaFreeHost (output_buf_mem );
  if (RFI_MODE == 2)
  {
    if (fft_out_kur) cudaFree (fft_out_kur);
    if (fft_ave_kur) cudaFree (fft_ave_kur);
    if (fft_trim_u_kur_dev) cudaFree (fft_trim_u_kur_dev);
    if (fft_trim_u_kur_hst) cudaFreeHost (fft_trim_u_kur_hst);
  }
  if (RFI_MODE)
  {
    if (kur_dev) cudaFree (kur_dev);
    if (kur_fb_dev) cudaFree (kur_fb_dev);
    if (kur_weights_dev) cudaFree (kur_weights_dev);
  }
  #if DOHISTO
  if (histo_dev) cudaFree (histo_dev);
  if (histo_hst) cudaFreeHost (histo_hst);
  #endif
  #if WRITE_KURTO
  if (kur_hst) cudaFreeHost (kur_hst);
  if (kur_fb_hst) cudaFreeHost (kur_fb_hst);
  if (kur_weights_hst) cudaFreeHost (kur_weights_hst);
  #endif

  // these are on the stack, don't think we want this
  //#if PROFILE
  //cudaEventDestroy (start);
  //cudaEventDestroy (stop);
  //#endif
  //cudaEventDestroy (obs_start);
  //cudaEventDestroy (obs_stop);

  return exit_status;
    
}

