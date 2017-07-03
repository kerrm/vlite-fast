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
}

static volatile int NBIT = 2;
static FILE* logfile_fp = NULL;


void usage ()
{
  fprintf (stdout,"Usage: process [options]\n"
	  "-k hexadecimal shared memory key for input (default: 40)\n"
	  "-K hexadecimal shared memory key for output (default: 0=disabled)\n"
	  "-p listening port number (default: %zu; if 0, disable)\n"
	  "-o print logging messages to stdout (as well as logfile)\n"
	  "-w output filterbank data (0=no sources, 1=listed sources, 2=all sources [def])\n"
	  "-b reduce output to b bits (2 [def], 4, and 8 are supported))\n"
	  "-P number of output polarizations (1=AA+BB, 2=AA,BB; 4 not implemented)\n"
	  "-r RFI-excision mode (0=no excision, 1=only excision, 2=[default]write both data sets)\n"
    "-s process a single observation, then quit (used for debugging)\n"
	  //"-m retain MUOS band\n"
	  ,(uint64_t)READER_SERVICE_PORT);
}

void sigint_handler (int dummy) {
  fclose (logfile_fp);
  // TODO other cleanup?
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

// TODO -- make sure this is consistent with sigproc, specifically for time
int write_psrdada_header (dada_hdu_t* hdu, char* inchdr, int npol, char* fb_file)
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
  char dada_utc[64];
  ascii_header_get (inchdr, "UTC_START", "%s", dada_utc);
  dada_utc[63] = '\0';

  // initialize observation parameters for filterbank
  // NB the data are upper sideband, so negative channel bandwidth
  double chbw = -64./NCHAN;
  double tsamp = double(NFFT)/VLITE_RATE*NSCRUNCH*1e6; // NB in mus
  int nchan = CHANMAX-CHANMIN+1;
  double bw = nchan*chbw;
  double freq0 = 384.;
  double freq = freq0 + 0.5*(CHANMIN+CHANMAX-1)*chbw;

  dada_hdu_lock_write (hdu);
  char* ascii_hdr = ipcbuf_get_next_write (hdu->header_block);
  dadacheck (ascii_header_set (ascii_hdr, "STATIONID", "%d", station_id));
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
  if (fb_file)
    dadacheck (ascii_header_set (ascii_hdr, "SIGPROC_FILE", "%s", fb_file) );
  multilog (hdu->log, LOG_INFO, "%s",ascii_hdr);
  ipcbuf_mark_filled (hdu->header_block, 4096);
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
  return 0;
}


void get_fbfile (char* fbfile, ssize_t fbfile_len, char* inchdr, vdif_header* vdhdr)
{
  // Open up filterbank file using timestamp and antenna
  int station_id = 0;
  ascii_header_get (inchdr, "STATIONID", "%d", &station_id);
  char currt_string[128];
  struct tm tm_epoch = {0};
  int vdif_epoch = getVDIFEpoch (vdhdr);
  tm_epoch.tm_year = 100 + vdif_epoch/2;
  tm_epoch.tm_mon = 6*(vdif_epoch%2);
  time_t epoch_seconds = mktime (&tm_epoch) + getVDIFFrameEpochSecOffset (vdhdr);
  struct tm* utc_time = gmtime (&epoch_seconds);
  strftime (currt_string,sizeof(currt_string), "%Y%m%d_%H%M%S", utc_time);
  *(currt_string+15) = 0;
  if (CHANMIN < 2411)
    snprintf (fbfile,fbfile_len,"%s/%s_muos_ea%02d.fil",DATADIR,currt_string,station_id);
  else
    snprintf (fbfile,fbfile_len,"%s/%s_ea%02d.fil",DATADIR,currt_string,station_id);
}

int main (int argc, char *argv[])
{
  // register SIGINT handling
  signal (SIGINT, sigint_handler);

  int exit_status = EXIT_SUCCESS;
  key_t key_in = 0x40;
  key_t key_out = 0x0;
  uint64_t port = READER_SERVICE_PORT;
  int stdout_output = 0;
  int write_fb = 2;
  int npol = 1;
  int RFI_MODE=2;
  int major_skip = 0;
  int single_pass = 0;  

  int arg = 0;
  while ((arg = getopt(argc, argv, "hk:K:p:omw:b:P:r:s")) != -1) {

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

    case 'K':
      if (sscanf (optarg, "%x", &key_out) != 1) {
        fprintf (stderr, "writer: could not parse key from %s\n", optarg);
        return -1;
      }
      break;

    case 'p':
      if (sscanf (optarg, "%zu", &port) != 1) {
        fprintf (stderr, "writer: could not parse port from %s\n", optarg);
        return -1;
      }
      break;

    case 'o':
      stdout_output = 1;
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

    }
  }

  cudacheck (cudaSetDevice (DEVICE));
  int nsms;
  cudaDeviceGetAttribute (&nsms,cudaDevAttrMultiProcessorCount,DEVICE);

  struct timespec ts_1ms = get_ms_ts (1);
  struct timespec ts_10s = get_ms_ts (10000);

  #if PROFILE
  // support for measuring run times of parts
  cudaEvent_t start,stop,start_total,stop_total;
  cudaEventCreate (&start);
  cudaEventCreate (&stop);
  cudaEventCreate (&start_total);
  cudaEventCreate (&stop_total);
  cudaEventRecord (start_total,0);
  float alloc_time=0, hdr_time=0, read_time=0, todev_time=0, 
      convert_time=0, kurtosis_time=0, fft_time=0, histo_time=0,
      normalize_time=0, tscrunch_time=0, pscrunch_time=0, 
      digitize_time=0, write_time=0, flush_time=0, misc_time=0,
      elapsed=0, total_time=0;
  #endif
  // measure full run time time
  cudaEvent_t obs_start, obs_stop;
  cudaEventCreate (&obs_start);
  cudaEventCreate (&obs_stop);

  multilog_t* log = multilog_open ("process_baseband",0);
  char logfile[128];
  time_t currt = time (NULL);
  struct tm tmpt; 
  localtime_r (&currt, &tmpt);
  char currt_string[32];
  strftime (currt_string,sizeof(currt_string), "%Y%m%d_%H%M%S", &tmpt);
  currt_string[15] = 0;
  char hostname[MAXHOSTNAME];
  gethostname (hostname,MAXHOSTNAME);
  pid_t pid = getpid();
  snprintf (logfile,128,
      "%s/%s_%s_process_%06d.log",LOGDIR,currt_string,hostname,pid);
  logfile_fp = myopen (logfile, "w");
  multilog_add (log, logfile_fp);
  if (stdout_output)
    multilog_add (log, stdout);

  // log configuration parameters
  char incoming_hdr[4096]; // borrow this buffer
  strncpy (incoming_hdr, argv[0], sizeof(incoming_hdr)-1);
  for (int i = 1; i < argc; ++i) {
    strcat (incoming_hdr, " ");
    strncat (incoming_hdr, argv[i], 4095-strlen(incoming_hdr));
  }
  multilog (log, LOG_INFO, "invoked with: \n%s\n", incoming_hdr);
  if (key_out)
    multilog (log, LOG_INFO, "Will write to %d psrdada buffers.\n",
      NOUTBUFF);

  // sanity checks on configuration parameters
  if (NKURTO != 250 && NKURTO != 500) {
    multilog (log, LOG_ERR, "Only NKURTO==250 or 500 supported.\n");
    exit (EXIT_FAILURE);
  }

  // connect to input buffer
  dada_hdu_t* hdu_in = dada_hdu_create (log);
  dada_hdu_set_key (hdu_in,key_in);
  if (dada_hdu_connect (hdu_in) != 0) {
    multilog (log, LOG_ERR, 
        "Unable to connect to incoming PSRDADA buffer!\n");
    exit (EXIT_FAILURE);
  }

  // connect to output buffer (optional)
  dada_hdu_t* hdu_out[NOUTBUFF] = {NULL};
  FILE* heimdall_fp[NOUTBUFF] = {NULL};
  int active_buffer = 0;
  if (key_out) {
    // connect to the output buffer(s)
    for (int i=0; i < NOUTBUFF; i++)
    {
      hdu_out[i] = dada_hdu_create (log);
      dada_hdu_set_key (hdu_out[i],key_out+i*100);
      if (dada_hdu_connect (hdu_out[i]) != 0) {
        multilog (log, LOG_ERR, 
            "Unable to connect to outgoing PSRDADA buffer #%d!\n",i+1);
        exit (EXIT_FAILURE);
      }
    }
  }

  #if PROFILE
  cudaEventRecord(start,0);
  #endif
  
  // allocate memory for 1 second; in the future, consider keeping
  // multiple buffers to allow for out-of-order data; this implementation
  // will handle dropped data
  unsigned char* udat; cudacheck (
  cudaMallocHost ((void**)&udat, 2*VLITE_RATE) );
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
  // FFTs per processing chunk
  int fft_per_chunk = samps_per_chunk / NFFT;

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

  // device memory for kurtosis statistics, uses NKURTO samples
  // NKURTO must be commensurate with samps_per_chunk/2
  size_t nkurto_per_chunk = samps_per_chunk / NKURTO;

  // extra factor of 2 to store both power and kurtosis statistics
  // storage for high time resolution and filterbank block ("fb") scales
  cufftReal* kur_dev=NULL; if (RFI_MODE) cudacheck (
  cudaMalloc ((void**)&kur_dev,2*sizeof(cufftReal)*nkurto_per_chunk) );
  cufftReal* kur_fb_dev=NULL; if (RFI_MODE) cudacheck (
  cudaMalloc ((void**)&kur_fb_dev,2*sizeof(cufftReal)*fft_per_chunk) );

  // store D'Agostino statistic for thresholding
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
    multilog (log, LOG_ERR, 
        "Selected channel and bit scheme is not commensurate!.\n");
    exit (EXIT_FAILURE);
  }

  // reduce array size by packing of samples into byte
  trim /= (8/NBIT);
  unsigned char* fft_trim_u; cudacheck (
  cudaMalloc ((void**)&fft_trim_u,trim) );
  unsigned char* fft_trim_u_hst; cudacheck (
  cudaMallocHost ((void**)&fft_trim_u_hst,trim) );

  unsigned char* fft_trim_u_kur=fft_trim_u;
  if (2 == RFI_MODE)
    cudacheck (cudaMalloc ((void**)&fft_trim_u_kur,trim) );
  unsigned char* fft_trim_u_kur_hst=fft_trim_u_hst;
  if (2 == RFI_MODE)
    cudacheck (cudaMallocHost ((void**)&fft_trim_u_kur_hst,trim) );

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

  // connect to control socket
  Connection conn;
  conn.sockoptval = 1; //release port immediately after closing connection
  if (port) {
    if (serve (port, &conn) < 0) {
      multilog (log, LOG_ERR,
          "Failed to create control socket on port %d.\n", port);
      exit (EXIT_FAILURE);
    }
    fcntl (conn.rqst, F_SETFL, O_NONBLOCK); // set up for polling
  }
  char cmd_buff[32];

  int quit = 0;

  // TEMP for profiling
  if (single_pass)
  {
    char cmd[256];
    snprintf (cmd, 255, "/home/vlite-master/mtk/src/readbase /home/vlite-master/mtk/baseband/20150302_151845_B0329+54_ev_0005.out.uw"); 
    multilog (log, LOG_INFO, "%s\n", cmd);
    popen (cmd, "r");
    nanosleep (&ts_10s,NULL);
  }

  // This point start to loop over observations.
  while(true)
  {

  if (quit)
    break;

  fflush (logfile_fp);
  dadacheck (dada_hdu_lock_read (hdu_in) );
  multilog (log, LOG_INFO, "Waiting for DADA header.\n");
  fflush (logfile_fp);

  // this will block until a header has been created by writer, signalling the start of
  // an observation
  uint64_t hdr_size = 0;
  char* ascii_hdr = ipcbuf_get_next_read (hdu_in->header_block,&hdr_size);
  multilog (log, LOG_INFO, "Received header with size %d.\n", hdr_size);

  if (ascii_hdr == NULL) {
    // this could be an error condition, or it could be a result of 
    // shutting down the data acquisition system; check to see if a
    // CMD_QUIT has been issued in order to log the appropriate outcome
    ssize_t npoll_bytes = 0;
    if (port) npoll_bytes = read (conn.rqst, cmd_buff, 32);
    if (npoll_bytes >  0) {
      for (int ib = 0; ib < npoll_bytes; ib++) {
        printf("Read command character %c.\n",cmd_buff[ib]);
        if (cmd_buff[ib] == CMD_QUIT) {
          multilog (log, LOG_INFO,
              "Received CMD_QUIT, indicating data taking is ceasing.  Exiting.\n");
          quit = 1;
          break;
        }
      }
      if (quit) break;
    }
    // did not receive CMD_QUIT, so this is an abnormal state
    multilog (log, LOG_ERR, "PSRDADA read failed unexpectedly.  Terminating.\n");
    exit_status = EXIT_FAILURE;
    break;
  }

  cudaEventRecord(obs_start,0);

  multilog (log, LOG_INFO, "psrdada header:\n%s",ascii_hdr);
  multilog (log, LOG_INFO, "Beginning new observation.\n\n");
  fflush (logfile_fp);
  memcpy (incoming_hdr,ascii_hdr,hdr_size);
  dadacheck (ipcbuf_mark_cleared (hdu_in->header_block));

  int bad_packets = 0;
  int skip_observation = 1;
  int frames_in_order[2] = {0,0};
  int last_sample[2] = {-1,-1};

  // consume from the buffer until we reach 1-s boundary

  // change this to consume until we reach 1-s boundary AND the previous
  // 1s of data are contiguous; this will hopefully prevent problems with
  // packet loss; if we haven't achieved a good chunk after 10 "seconds"
  // worth of packets, abort the observation
  for (int i=0; i < VLITE_FRAME_RATE*20; ++i)
  {
    ssize_t nread = ipcio_read (hdu_in->data_block,frame_buff,VD_FRM);
    if (nread != VD_FRM)
    {
      bad_packets++;
      continue;
    }
    if (nread < 0)
    {
      multilog (log, LOG_ERR, "Problem with nread! Skipping observation.\n");
      skip_observation = 1;
      break;
    }
    if (0==i)
      multilog (log, LOG_INFO,"Starting trim at frame %d and second %d.\n",
          getVDIFFrameNumber(vdhdr),getVDIFFrameSecond(vdhdr));

    int thread = getVDIFThreadID (vdhdr) != 0;
    int sample = getVDIFFrameNumber (vdhdr);

    // at the beginning of a second, check for contiguity
    if (0==sample)
    {
      // we have passed through and have acquired a full second of data
      if ((VLITE_FRAME_RATE-1)==frames_in_order[0] &&
          (VLITE_FRAME_RATE-1)==frames_in_order[1])
      {
        skip_observation = 0;
        break;
      }
      else
        frames_in_order[thread] = 0;
    }
    else if (sample==last_sample[thread]+1)
      frames_in_order[thread]++;

    last_sample[thread] = sample;

  }
  multilog (log, LOG_INFO, "hdu_in=%p, hdu_in->data_block=%p\n",hdu_in,hdu_in->data_block);

  if (bad_packets > 0)
    multilog (log, LOG_INFO,
      "Exiting buffer consumption with %d bad packets.\n",bad_packets);

  if (skip_observation)
  {
    multilog (log, LOG_INFO, "Skipping observation.\n");
    dadacheck (dada_hdu_unlock_read (hdu_in));
    continue;
  }

  multilog (log, LOG_INFO, "hdu_in=%p, hdu_in->data_block=%p\n",hdu_in,hdu_in->data_block);
  nanosleep (&ts_1ms,NULL);
  nanosleep (&ts_1ms,NULL);
  nanosleep (&ts_1ms,NULL);
  nanosleep (&ts_1ms,NULL);
  nanosleep (&ts_1ms,NULL);
  multilog (log, LOG_INFO, "hdu_in=%p, hdu_in->data_block=%p\n",hdu_in,hdu_in->data_block);

  // Open up filterbank file at appropriate time
  char fbfile[256];
  multilog (log, LOG_INFO, "hdu_in=%p, hdu_in->data_block=%p\n",hdu_in,hdu_in->data_block);
  get_fbfile (fbfile, 256, incoming_hdr, vdhdr);
  multilog (log, LOG_INFO, "hdu_in=%p, hdu_in->data_block=%p\n",hdu_in,hdu_in->data_block);
  char fbfile_kur[256];
  change_extension (fbfile, fbfile_kur, ".fil", "_kur.fil");

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
  multilog (log, LOG_INFO, "hdu_in=%p, hdu_in->data_block=%p\n",hdu_in,hdu_in->data_block);

  // check for write to /dev/null
  int write_to_null =  0==write_fb;
  if (write_fb == 1) {
    char name[OBSERVATION_NAME_SIZE];
    ascii_header_get (incoming_hdr, "NAME", "%s", name);
    name[OBSERVATION_NAME_SIZE-1] = '\0';
    if (check_name (name)) {
      multilog (log, LOG_INFO,
          "Source %s matches target list, recording filterbank data.\n",
          name);
    }
    else {
      multilog (log, LOG_INFO,
          "Source %s not on target list, disabling filterbank data.\n",
          name);
      write_to_null = 1;
    }
  }
  if (write_to_null)
  {
    multilog (log, LOG_INFO,
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
  multilog (log, LOG_INFO, "hdu_in=%p, hdu_in->data_block=%p\n",hdu_in,hdu_in->data_block);

  FILE *fb_fp,*fb_kurto_fp=NULL;
  uint64_t fb_bytes_written = 0;
  if (RFI_MODE == 0 || RFI_MODE == 2 )
  {
    fb_fp = myopen (fbfile, "wb", true);
    multilog (log, LOG_INFO,
        "Writing no-RFI-excision filterbanks to %s.\n",fbfile);
  }
  else
  {
    fb_fp = myopen (fbfile_kur, "wb", true);
    multilog (log, LOG_INFO,
        "Writing RFI-excision filterbanks to %s.\n",fbfile_kur);
  }
  if (RFI_MODE == 2)
  {
    fb_kurto_fp = myopen (fbfile_kur, "wb", true);
    multilog (log, LOG_INFO,
        "Writing RFI-excision filterbanks to %s.\n",fbfile_kur);
  }
  multilog (log, LOG_INFO, "hdu_in=%p, hdu_in->data_block=%p\n",hdu_in,hdu_in->data_block);

  #if DOHISTO
  FILE *histo_fp = myopen (histofile, "wb", true);
  multilog (log, LOG_INFO, "Writing histograms to %s.\n",histofile);
  #endif

  #if WRITE_KURTO
  FILE *kurto_fp = myopen (kurtofile, "wb", true);
  multilog (log, LOG_INFO, "Writing kurtosis to %s.\n",kurtofile);
  FILE *block_kurto_fp = myopen (block_kurtofile, "wb", true);
  multilog (
       log, LOG_INFO, "Writing block kurtosis to %s.\n",block_kurtofile);
  FILE *weight_fp = myopen (weightfile, "wb", true);
  multilog (log, LOG_INFO, "Writing weights to %s.\n",weightfile);
  #endif
  fflush (logfile_fp);
  multilog (log, LOG_INFO, "hdu_in=%p, hdu_in->data_block=%p\n",hdu_in,hdu_in->data_block);

  #if PROFILE
  cudaEventRecord(start,0);
  #endif

  if (key_out)
    //Write psrdada output header values
    write_psrdada_header (hdu_out[active_buffer], incoming_hdr, npol, 
      RFI_MODE?fbfile_kur:fbfile);

  // write out a sigproc header
  write_sigproc_header (fb_fp, incoming_hdr, vdhdr, npol);
  if (2 == RFI_MODE)
    write_sigproc_header (fb_kurto_fp, incoming_hdr, vdhdr, npol);

  #if PROFILE
  CUDA_PROFILE_STOP(start,stop,&hdr_time)
  #endif

  // launch heimdall
  if (key_out) {
    if (heimdall_fp[active_buffer])
    {
      // have used previous buffer, make sure processing is finished
      // NB this runs the risk of a data skip if it hasn't, but will 
      // check for that error condition later
      if (pclose (heimdall_fp[active_buffer]) != 0)
      {
        multilog (log, LOG_ERR, "heimdall exit status nonzero\n");
        // TODO -- exit on such an error?
      }
      heimdall_fp[active_buffer] = NULL;
    }
    char heimdall_cmd[256];
    snprintf (heimdall_cmd, 255, "/home/vlite-master/mtk/bin/heimdall -nsamps_gulp 62500 -gpu_id 0 -dm 2 1000 -boxcar_max 64 -group_output -zap_chans 0 190 -k %x -v", key_out + active_buffer*100); 
    multilog (log, LOG_INFO, "%s\n", heimdall_cmd);
    heimdall_fp[active_buffer] = popen (heimdall_cmd, "r");
  }

  // set up sample# trackers and copy first frame into 1s-buffer
  int current_sec = getVDIFFrameSecond(vdhdr);
  last_sample[0] = last_sample[1] = -1;
  int current_thread = getVDIFThreadID (vdhdr) != 0;
  multilog (log, LOG_INFO, "Starting sec=%d, thread=%d\n",current_sec,current_thread);

  bool skipped_data_announced = false;

  while(true) // loop over data packets
  {
    int thread = getVDIFThreadID (vdhdr) != 0;
    int sample = getVDIFFrameNumber (vdhdr);
    int second = getVDIFFrameSecond (vdhdr);

    if (second==current_sec || major_skip)
    {

      #if PROFILE
      cudaEventRecord (start,0);
      #endif 

      size_t idx = VLITE_RATE*thread + sample*VD_DAT;
      memcpy (udat + idx,dat,VD_DAT);

      // check for skipped data
      if (sample!=last_sample[thread] + 1)
      {
        if (!skipped_data_announced) {
          multilog (log, LOG_INFO, "Found skipped data in current second!\n");
          skipped_data_announced = true;
        }
        // set the intervening memory to zeros
        // TODO -- not sure this is correct, and not sure what to set
        // memory to, if anything; probably should just abort
        //idx = VLITE_RATE*thread + (last_sample[thread]+1)*VD_DAT;
        //memset (udat+idx,0,VD_DAT*(sample-last_sample[thread]-1));
      }
      last_sample[thread] = sample;

      // read the next frame into the buffer
      ssize_t nread = ipcio_read (hdu_in->data_block,frame_buff,VD_FRM);

#if PROFILE
      CUDA_PROFILE_STOP(start,stop,&elapsed)
      read_time += elapsed;
#endif

      if (nread < 0 ) {
        multilog (log, LOG_ERR, "Error on nread.\n");
        exit (EXIT_FAILURE);
      }
      if (nread != VD_FRM) {
        if (nread != 0)
          multilog (log, LOG_INFO,
            "Packet size=%ld, expected %ld.  Aborting this observation.\n",
            nread,VD_FRM);
        break;
      }

      continue;

    }

    if (second - current_sec > 1)
    {
      // this should only happen when observing over epoch changes,
      // leap seconds, or when there are major data drops.  When that
      // happens, just restart the code.
      multilog (log, LOG_ERR,
        "Major data skip!  (%d vs. %d; thread = %d) Aborting this observation.\n",
         second,current_sec,thread);
         major_skip = 1;
         break;
    }

    // every 1s, check for a QUIT command
    ssize_t npoll_bytes = 0;
    if (port) npoll_bytes = read (conn.rqst, cmd_buff, 32);
    if (npoll_bytes >  0) {
      for (int ib = 0; ib < npoll_bytes; ib++) {
        printf("Read command character %c.\n",cmd_buff[ib]);
        if (cmd_buff[ib] == CMD_QUIT) {
          quit = 1;
          break;
        }
      }
    }
    if (quit) {
      multilog (log, LOG_INFO, "Received CMD_QUIT.  Exiting!.\n");
      // this will break out of data loop, write out file, then break out of observation
      // loop to exit program
      break;
    }
    // keep log on disk up to date
    fflush (logfile_fp);

    // check that both buffers are full
    skipped_data_announced = false;
    current_sec = second;

    // need to dispatch the current buffer; zero fill any discontinuity
    // at the end; TODO is this redundant?  Already checking for skipped
    // data; I guess at this point we might have more of an idea about the
    // magnitude of a skip
    for (int i = 0; i < 2; ++i)
    {
      int target = VLITE_FRAME_RATE - 1;
      if (last_sample[i] !=  target)
      {
        multilog (log, LOG_ERR,
          "Found skipped data in dispatch!  Aborting this observation.\n");
        //size_t idx = VLITE_RATE*i + last_sample[i]*VD_DAT;
        //memset(udat+idx,0,(target-last_sample[i])*VD_DAT);
        break;
      }
      last_sample[i] = -1;
    }

    // do dispatch -- break into chunks to fit in GPU memory
    for (int iseg = 0; iseg < SEG_PER_SEC; iseg++)
    {

      #if PROFILE
      cudaEventRecord(start,0);
      #endif 
      // need to copy pols separately
      for (int ipol=0; ipol < 2; ++ipol)
      {
        unsigned char* dst = udat_dev + ipol*samps_per_chunk/2;
        unsigned char* src = udat + ipol*VLITE_RATE + iseg*samps_per_chunk/2;
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

        kurtosis <<<nkurto_per_chunk, 256>>> (
            fft_in,kur_dev,kur_dev+nkurto_per_chunk);
        cudacheck (cudaGetLastError () );

        // compute the thresholding statistic
        compute_dagostino <<<nsms*32,NTHREAD>>> (
            kur_dev+nkurto_per_chunk,dag_dev,nkurto_per_chunk);
        cudacheck (cudaGetLastError () );

        // do the same for block level
        block_kurtosis <<<fft_per_chunk/8,256>>> (
            kur_dev,kur_dev+nkurto_per_chunk,dag_dev,
            kur_fb_dev,kur_fb_dev+fft_per_chunk);
        cudacheck (cudaGetLastError () );
        compute_dagostino2 <<<nsms*32,NTHREAD>>> (
            kur_fb_dev+fft_per_chunk,dag_fb_dev,fft_per_chunk);
        cudacheck (cudaGetLastError () );

        cudacheck (cudaMemset (
            kur_weights_dev,0,sizeof(cufftReal)*fft_per_chunk) );
        // NB that fft_in_kur==fft_in if not writing both streams
        apply_kurtosis <<<nkurto_per_chunk, 256>>> (
            fft_in,fft_in_kur,dag_dev,dag_fb_dev,kur_weights_dev);
        cudacheck (cudaGetLastError () );

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

      ////// NORMALIZE BANDPASS //////
      #if PROFILE
      cudaEventRecord(start,0);
      #endif 
      if (RFI_MODE == 0 || RFI_MODE == 2)
      {
        detect_and_normalize2 <<<(NCHAN*2)/NTHREAD+1,NTHREAD>>> (
            fft_out,bp_dev,bp_scale,fft_per_chunk/2);
        cudacheck ( cudaGetLastError () );
      }
      if (RFI_MODE == 1 || RFI_MODE == 2)
      {
        detect_and_normalize3 <<<(NCHAN*2)/NTHREAD+1,NTHREAD>>> (
            fft_out_kur,kur_weights_dev,bp_kur_dev,bp_scale,
            fft_per_chunk/2);
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
      size_t maxn = (fft_per_chunk*NCHAN)/polfac;
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
              fft_ave,fft_trim_u,maxn, npol);
          if (RFI_MODE > 1)
            sel_and_dig_2b <<<nsms*32,NTHREAD>>> (
                fft_ave_kur,fft_trim_u_kur,maxn, npol);
          break;
        case 4:
          sel_and_dig_4b <<<nsms*32,NTHREAD>>> (
              fft_ave,fft_trim_u,maxn, npol);
          if (RFI_MODE > 1)
            sel_and_dig_4b <<<nsms*32,NTHREAD>>> (
                fft_ave_kur,fft_trim_u_kur,maxn, npol);
          break;
        case 8:
          sel_and_dig_8b <<<nsms*32,NTHREAD>>> (
              fft_ave,fft_trim_u,maxn, npol);
          if (RFI_MODE > 1)
          sel_and_dig_8b <<<nsms*32,NTHREAD>>> (
              fft_ave_kur,fft_trim_u_kur,maxn, npol);
          break;
        default:
          sel_and_dig_2b <<<nsms*32,NTHREAD>>> (
              fft_ave,fft_trim_u,maxn, npol);
          if (RFI_MODE > 1)
            sel_and_dig_2b <<<nsms*32,NTHREAD>>> (
                fft_ave_kur,fft_trim_u_kur,maxn, npol);
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

      // copy filterbanked data back to host
      cudacheck (cudaMemcpy (
            fft_trim_u_hst,fft_trim_u,maxn,cudaMemcpyDeviceToHost) );
      if (RFI_MODE > 1)
        cudacheck (cudaMemcpy ( fft_trim_u_kur_hst,
            fft_trim_u_kur,maxn,cudaMemcpyDeviceToHost) );

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

      // finally, push the filterbanked time samples onto psrdada buffer
      // and/or write out to sigproc

      #if PROFILE
      cudaEventRecord (start,0);
      #endif 

      if (key_out)
      {
        if (RFI_MODE == 2)
          ipcio_write (hdu_out[active_buffer]->data_block,
              (char *)fft_trim_u_kur_hst,maxn);
        else
          ipcio_write (hdu_out[active_buffer]->data_block,
              (char *)fft_trim_u_hst,maxn);
      }

      fwrite (fft_trim_u_hst,1,maxn,fb_fp);
      fb_bytes_written += maxn;
      if (RFI_MODE == 2)
        fwrite (fft_trim_u_kur_hst,1,maxn,fb_kurto_fp);

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

    }
  } // end loop over packets

  #if PROFILE
  cudaEventRecord(start,0);
  #endif 

  dadacheck (dada_hdu_unlock_read (hdu_in));
  if (key_out)
  {
    dadacheck (dada_hdu_unlock_write (hdu_out[active_buffer]));
    // switch to next output buffer
    if (NOUTBUFF == ++active_buffer)
      active_buffer = 0;
  }

  // close files
  fclose (fb_fp);
  if (fb_kurto_fp) fclose (fb_kurto_fp);
  #if DOHISTO
  fclose (histo_fp);
  #endif
  #if WRITE_KURTO
  fclose (kurto_fp);
  fclose (block_kurto_fp);
  fclose (weight_fp);
  #endif
  uint64_t samps_written = (fb_bytes_written*(8/NBIT))/(CHANMAX-CHANMIN+1);
  multilog (log, LOG_INFO, "Wrote %.2f MB (%.2f s) to %s\n",
      fb_bytes_written*1e-6,samps_written*tsamp,fbfile);

  float obs_time;
  CUDA_PROFILE_STOP(obs_start,obs_stop,&obs_time);
  multilog (log, LOG_INFO, "Proc Time...%.3f\n", obs_time*1e-3);

  // TODO -- this can't happen, because we've already disconnected
  // from the input buffer!  Moreover, it seems to be an infinite
  // loop; take it out for now
  if (0) {
  // If there has been a major skip, need to continue to read from
  // buffer so it doesn't fill up
  //if (major_skip) {
    multilog (log, LOG_INFO, "Reading from buffer to clear major skip.\n");
    major_skip = 0;
    while (true) {
      ssize_t nread = ipcio_read (hdu_in->data_block,frame_buff,VD_FRM);
      if (nread < 0 ) {
        multilog (log, LOG_ERR, "Error on nread on major skip clear.\n");
        // TODO -- exit more gracefully here?
        return (EXIT_FAILURE);
      }
      if (nread != VD_FRM) {
        if (nread != 0)
          multilog (log, LOG_INFO,
            "Packet size=%ld, expected %ld on major skip clear.\n",
            nread,VD_FRM);
        break;
      }
    }
  }

  // TODO -- if it was a scan observation and not long enough to complete a full cycle,, delete the file
  /*
  if (fb_bytes_written == 0) {
    remove (fname);
  }
  */

  #if PROFILE
  CUDA_PROFILE_STOP(start,stop,&flush_time)
  CUDA_PROFILE_STOP(start_total,stop_total,&total_time)
  float sub_time = alloc_time + hdr_time + read_time + todev_time + 
      histo_time + convert_time + kurtosis_time + fft_time + 
      normalize_time + pscrunch_time + tscrunch_time + digitize_time + 
      write_time + flush_time + misc_time;
  multilog (log, LOG_INFO, "Alloc Time..%.3f\n", alloc_time*1e-3);
  multilog (log, LOG_INFO, "Read Time...%.3f\n", read_time*1e-3);
  multilog (log, LOG_INFO, "Copy To Dev.%.3f\n", todev_time*1e-3);
  multilog (log, LOG_INFO, "Histogram...%.3f\n", histo_time*1e-3);
  multilog (log, LOG_INFO, "Convert.....%.3f\n", convert_time*1e-3);
  multilog (log, LOG_INFO, "Kurtosis....%.3f\n", kurtosis_time*1e-3);
  multilog (log, LOG_INFO, "FFT.........%.3f\n", fft_time*1e-3);
  multilog (log, LOG_INFO, "Normalize...%.3f\n", normalize_time*1e-3);
  multilog (log, LOG_INFO, "Pscrunch....%.3f\n", pscrunch_time*1e-3);
  multilog (log, LOG_INFO, "Tscrunch....%.3f\n", tscrunch_time*1e-3);
  multilog (log, LOG_INFO, "Digitize....%.3f\n", digitize_time*1e-3);
  multilog (log, LOG_INFO, "Write.......%.3f\n", write_time*1e-3);
  multilog (log, LOG_INFO, "Flush.......%.3f\n", flush_time*1e-3);
  multilog (log, LOG_INFO, "Misc........%.3f\n", misc_time*1e-3);
  multilog (log, LOG_INFO, "Sum of subs.%.3f\n", sub_time*1e-3);
  multilog (log, LOG_INFO, "Total run...%.3f\n", total_time*1e-3);
  #endif

  if (single_pass)
    break;

  } // end loop over observations

  // close sockets and files
  shutdown (conn.rqst, 2);

  // clean up psrdada data structures
  dada_hdu_disconnect (hdu_in);
  dada_hdu_destroy (hdu_in);
  if (key_out) {
    for (int i=0; i < NOUTBUFF; i++)
    {
      dada_hdu_disconnect (hdu_out[i]);
      dada_hdu_destroy (hdu_out[i]);
    }
  }

  //free memory
  cudaFreeHost (udat);
  cudaFree (udat_dev);
  cudaFree (fft_in);
  cudaFree (fft_out);
  cudaFree (fft_ave);
  cudaFree (fft_trim_u);
  cudaFreeHost (fft_trim_u_hst);
  if (RFI_MODE == 2)
  {
    cudaFree (fft_out_kur);
    cudaFree (fft_ave_kur);
    cudaFree (fft_trim_u_kur);
    cudaFreeHost (fft_trim_u_kur_hst);
  }
  if (RFI_MODE)
  {
    cudaFree (kur_dev);
    cudaFree (kur_fb_dev);
    cudaFree (kur_weights_dev);
  }
  #if DOHISTO
  cudaFree (histo_dev);
  cudaFreeHost (histo_hst);
  #endif
  #if WRITE_KURTO
  cudaFreeHost (kur_hst);
  cudaFreeHost (kur_fb_hst);
  cudaFreeHost (kur_weights_hst);
  #endif

  #if PROFILE
  cudaEventDestroy (start);
  cudaEventDestroy (stop);
  cudaEventDestroy (start_total);
  cudaEventDestroy (stop_total);
  #endif
  cudaEventDestroy (obs_start);
  cudaEventDestroy (obs_stop);

  multilog (log, LOG_INFO, "Completed shutdown.");
  multilog_close (log);
  fclose (logfile_fp);

  return exit_status;
    
}

// quantities for D'Agostino normality test (see wikipedia)
#define NK float(NKURTO)
#define mu1 (-6./(NK+1))
#define mu2 ((24.*NK*(NK-2)*(NK-3))/((NK+1)*(NK+1)*(NK+3)*(NK+5)))
#define g1 (6.*(NK*NK-5*NK+2)/((NK+7)*(NK+9))*sqrt( (6.*(NK+3)*(NK+5))/(NK*(NK-2)*(NK-3)) ))
#define A (6.+(8./g1)*(2./g1 + sqrt(1. + 4./(g1*g1))))
#define Z2_1 sqrt(4.5*A)
#define Z2_2 (1-2./(9*A))
#define Z2_3 sqrt(2./(mu2*(A-4)))

#define NKb float(NFFT)
#define mu1b (-6./(NKb+1))
#define mu2b ((24.*NKb*(NKb-2)*(NKb-3))/((NKb+1)*(NKb+1)*(NKb+3)*(NKb+5)))
#define g1b (6.*(NKb*NKb-5*NKb+2)/((NKb+7)*(NKb+9))*sqrt( (6.*(NKb+3)*(NKb+5))/(NKb*(NKb-2)*(NKb-3)) ))
#define Ab (6.+(8./g1b)*(2./g1b + sqrt(1. + 4./(g1b*g1b))))
#define Z2b_1 sqrt(4.5*Ab)
#define Z2b_2 (1-2./(9*Ab))
#define Z2b_3 sqrt(2./(mu2b*(Ab-4)))

//convert unsigned char time array to float
__global__ void convertarray (cufftReal *time, unsigned char *utime, size_t n)
{
  for (int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; i += blockDim.x*gridDim.x)
  {
    if (utime[i] == 0) 
      time[i] = 0;
    else
      time[i] = (cufftReal)(utime[i])/128-1;
  }
}

__global__ void kurtosis (cufftReal *time, cufftReal *pow, cufftReal *kur)
{
  // calculate the variance (power) and kurtosis for voltage statistics in
  // relatively short windows.  Do this by first copying from global memory,
  // then using a hard-coded tree reduction.  Right now, it is set up to
  // use either 250 or 500 samples, so must be invoked with either 256 
  // or 512 threads.

  // because each thread block works on a chunk of data that's commensurate
  // with the packing of samples into the buffer, specifically with respect
  // to the two polarizations, I think we don't need to worry at all about
  // a thread block crossing the polarization.  The output will simply
  // contain the statistics for pol 0 first, then pol 1.

  volatile __shared__ float data2[256];
  volatile __shared__ float data4[256];
  unsigned int tid = threadIdx.x;
  size_t offset = blockIdx.x*NKURTO;

  if (tid < 250)
  {
    if (NKURTO==500) {
      // load up two values from global memory in this case
      data2[tid] = time[offset + tid]*time[offset + tid];
      float tmp = time[offset + tid + 250]*time[offset + tid + 250];
      data4[tid] = data2[tid]*data2[tid] + tmp*tmp;
      data2[tid] += tmp;
    }
    else {
      data2[tid] = time[offset + tid]*time[offset + tid];
      data4[tid] = data2[tid]*data2[tid];
    }
  }
  else
    data2[tid] = data4[tid] = 0;
  __syncthreads ();

  if (tid < 128)
  {
    data2[tid] += data2[tid + 128];
    data4[tid] += data4[tid + 128];
  }
  __syncthreads ();

  if (tid < 64)
  {
    data2[tid] += data2[tid + 64];
    data4[tid] += data4[tid + 64];
  }
  __syncthreads ();

  if (tid < 32)
  {
    data2[tid] += data2[tid + 32];
    data4[tid] += data4[tid + 32];
    data2[tid] += data2[tid + 16];
    data4[tid] += data4[tid + 16];
    data2[tid] += data2[tid + 8];
    data4[tid] += data4[tid + 8];
    data2[tid] += data2[tid + 4];
    data4[tid] += data4[tid + 4];
    data2[tid] += data2[tid + 2];
    data4[tid] += data4[tid + 2];
  }

  if (tid==0)
  {
    data2[tid] += data2[tid + 1];
    data4[tid] += data4[tid + 1];
    pow[blockIdx.x] = data2[0]/NKURTO;
    kur[blockIdx.x] = data4[0]/NKURTO/(pow[blockIdx.x]*pow[blockIdx.x]);
  }
}

__global__ void compute_dagostino (cufftReal* kur, cufftReal* dag, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    if (kur[i] != 0.)
    {
      float t = (1-2./A)/(1.+(kur[i]-3.-mu1)*Z2_3);
      dag[i] = Z2_1*(Z2_2 - powf (t,1./3));
    }
    else
      dag[i] = 0.;
  }
}

// TODO -- this isn't quite right, because there won't necessarily be
// NFFT samples in the weighted version; since we're computing many fewer,
// don't need to precompute.  However, empirically it doesn't seem to make
// much of a difference..
__global__ void compute_dagostino2 (cufftReal* kur, cufftReal* dag, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    if (kur[i] != 0.)
    {
      float t = (1-2./Ab)/(1.+(kur[i]-3.-mu1b)*Z2b_3);
      dag[i] = Z2b_1*(Z2b_2 - powf (t,1./3));
    }
    else
      dag[i] = DAG_FB_THRESH + 1;
  }
}

// compute a filter-bank level statistic
// TODO -- can use the weighted sum here rather than the atomicAdds in
// apply kurtosis
__global__ void block_kurtosis (cufftReal* pow, cufftReal* kur, cufftReal* dag, cufftReal* pow_block, cufftReal* kur_block)
{
  volatile __shared__ float data2[256];
  volatile __shared__ float data4[256];
  volatile __shared__ unsigned char wt[256];
  // run with 256 threads; break it up such that we either do 5 blocks (for
  // NKURTO==500) or 10 blocks (for NKURTO=250)
  unsigned int tid = threadIdx.x;
  unsigned int warp_id = tid / 32;
  unsigned int warp_tid = tid - warp_id*32;
  if (warp_tid > 24) 
  {
    data2[tid] = 0;
    data4[tid] = 0;
    wt[tid] = 0;
  }
  else
  {
    // each thread block does 8 filterbank blocks (one for each warp)
    int idx = (blockIdx.x*8 + warp_id)*(NFFT/NKURTO) + warp_tid;
    wt[tid] = (dag[idx]<DAG_THRESH) && (dag[idx]>-DAG_THRESH);
    data2[tid] = wt[tid]*pow[idx];
    data4[tid] = wt[tid]*kur[idx]*pow[idx]*pow[idx];
    if (NKURTO==250)
    {
      __syncthreads ();
      idx += 25;
      float w = (dag[idx]<DAG_THRESH) && (dag[idx]>-DAG_THRESH);
      data2[tid] += w*pow[idx];
      data4[tid] += w*kur[idx]*pow[idx]*pow[idx];
      wt[tid] += w;
    }
  }

  if (warp_tid > 15)
    return;

  // do sum within each warp
  data2[tid] += data2[tid + 16];
  data4[tid] += data4[tid + 16];
  wt[tid]    += wt[tid + 16];
  data2[tid] += data2[tid + 8];
  data4[tid] += data4[tid + 8];
  wt[tid]    += wt[tid + 8];
  data2[tid] += data2[tid + 4];
  data4[tid] += data4[tid + 4];
  wt[tid]    += wt[tid + 4];
  data2[tid] += data2[tid + 2];
  data4[tid] += data4[tid + 2];
  wt[tid]    += wt[tid + 2];
  data2[tid] += data2[tid + 1];
  data4[tid] += data4[tid + 1];
  wt[tid]    += wt[tid + 1];

  if (0==warp_tid)
  {
    if (wt[tid] > 0) 
    {
      float p = pow_block[blockIdx.x*8+warp_id] = data2[tid]/wt[tid];
      kur_block[blockIdx.x*8+warp_id] = data4[tid]/wt[tid]/(p*p);
    }
    else 
    {
      pow_block[blockIdx.x*8+warp_id] = 0;
      kur_block[blockIdx.x*8+warp_id] = 0;
    }
  }
}

__global__ void apply_kurtosis (cufftReal *in, cufftReal *out, cufftReal *dag, cufftReal *dag_fb, cufftReal* norms)
{

  unsigned int tid = threadIdx.x;

  // D'Agostino kurtosis TS
  float d = dag[blockIdx.x];
  float d_fb = dag_fb[blockIdx.x/(NFFT/NKURTO)];
  #ifdef DEBUG_WEIGHTS
  // if debugging, set the weights to 0 for the second half of all samples in
  // the chunk for 2nd pol and for the final eighth for the 1st pol
  int time_idx = blockIdx.x * NKURTO;
  bool c1 = time_idx > 3*(VLITE_RATE/(SEG_PER_SEC*2));
  bool c2 = (time_idx < VLITE_RATE/SEG_PER_SEC) && (time_idx > (7*VLITE_RATE/SEG_PER_SEC)/8);
  d = d_fb = c1||c2?DAG_THRESH + 1:0;
  #endif
  if ((d > DAG_THRESH) || (d < -DAG_THRESH) || (d_fb > DAG_FB_THRESH) || (d < -DAG_FB_THRESH)) {
  //if ((d > DAG_THRESH) || (d < -DAG_THRESH)) {
    // zero voltages
    size_t offset = blockIdx.x*NKURTO;
    if (tid < 250) {
      out[offset + tid] =  0;
      if (NKURTO==500)
        out[offset + tid + 250] =  0;
    }
  }
  else
  {
    // if copying data, copy it
    if (in != out && tid < 250) {
      size_t offset = blockIdx.x*NKURTO;
      out[offset + tid] = in[offset + tid];
      if (NKURTO==500)
        out[offset + tid + 250] = in[offset + tid + 250];
    }

    // add one to the filterbank block samples
    if (tid==0)
      atomicAdd (norms + (blockIdx.x*NKURTO)/NFFT, float(NKURTO)/NFFT);
  }
}


__global__ void histogram ( unsigned char *utime, unsigned int* histo, size_t n)
{
  __shared__ unsigned int lhisto[512];
  lhisto[threadIdx.x] = 0;
  __syncthreads ();

  int i = threadIdx.x + blockIdx.x*blockDim.x;
  for (; i < n/2; i += blockDim.x*gridDim.x)
    atomicAdd (lhisto+utime[i], 1);
  for (; i < n; i += blockDim.x*gridDim.x)
    atomicAdd ((lhisto+256)+utime[i], 1);
  __syncthreads ();

  // MUST run with 512 threads for this global accumulation to work
  atomicAdd ( histo+threadIdx.x, lhisto[threadIdx.x]);
}

__global__ void detect_and_normalize2 (cufftComplex *fft_out, cufftReal* bp, 
        float scale, size_t ntime)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x; 
  if (i >= NCHAN*2) return;
  if (i >= NCHAN) // advance pointer to next polarization
  {
    fft_out += ntime*NCHAN;
    bp += NCHAN;
    i -= NCHAN;
  }

  // initialize bandpass to mean of first block
  float bp_l = bp[i];
  if (0. == bp_l) {
    for (int j = i; j < ntime*NCHAN; j+= NCHAN)
      bp_l += fft_out[j].x*fft_out[j].x + fft_out[j].y*fft_out[j].y;
    bp_l /= ntime;
  }

  for (int j = i; j < ntime*NCHAN; j+= NCHAN)
  {
    // detect
    float pow = fft_out[j].x*fft_out[j].x + fft_out[j].y*fft_out[j].y;

    // update bandpass
    bp_l = scale*pow + (1-scale)*bp_l;

    // scale to bandpass and mean-subtract; this assumes the powers are 
    // chi^2_2 distributed, var(pow)=4, mean(pow)=std(pow)=2.  Therefore
    // dividing by mean will give standard deviation of 1 centred at 1.
    fft_out[j].x = pow/bp_l-1;

  }
  // write out current bandpass
  bp[i] = bp_l;
}

__global__ void detect_and_normalize3 (cufftComplex *fft_out, cufftReal* kur_weights_dev, cufftReal* bp, float scale, size_t ntime)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x; 
  if (i >= NCHAN*2) return;
  if (i >= NCHAN) // advance pointer to next polarization
  {
    fft_out += ntime*NCHAN;
    kur_weights_dev += ntime;
    bp += NCHAN;
    i -= NCHAN;
  }

  // initialize bandpass to mean of first block
  float bp_l = bp[i];
  if (0. == bp_l) {
    int good_samples = 0;
    for (int j = i, time_idx=0; j < ntime*NCHAN; j+= NCHAN,time_idx++) {
      float w = kur_weights_dev[time_idx];
      if (0.==w)
        continue;
      good_samples++;
      bp_l += (fft_out[j].x*fft_out[j].x + fft_out[j].y*fft_out[j].y)/w;
    }
    if (0==good_samples) {
      // entire first block is bad; not sure what is best, try setting to
      // 1 and hope for the best?
      bp_l = 1;
    }
    else
      bp_l /= good_samples;
  }

  for (int j = i, time_idx=0; j < ntime*NCHAN; j+= NCHAN,time_idx++)
  {
    // detect
    //float w = kur_weights_dev[time_idx]*kur_weights_dev[time_idx];
    // NB that this formulation works because the weights are 0 or 1; if
    // we write out the expectation for the Fourier transform of the voltage
    // squared, the weights go in squared, so we normalize by the sum over
    // the weights squared, which is the same as the sum of the weights
    // (here kur_weights_dev) since they are 0 or 1
    float w = kur_weights_dev[time_idx];

    if (0.==w) {
      // if no samples are available, replace with mean bandpass
      fft_out[j].x = 0;
    }

    else {

      float pow = (fft_out[j].x*fft_out[j].x + fft_out[j].y*fft_out[j].y)/w;

      // apply a rough filter; values in excess of 11xmean shouldn't happen
      // more often than every 1.5 s, so we can clip values above this
      // without substantial distortion and possibly prevent bandpass
      // saturation; when we do, don't update the bandpass

      if (pow > bp_l*11)
        fft_out[j].x = 10;

      else {

        // update bandpass
        bp_l = scale*pow + (1-scale)*bp_l;

        // scale to bandpass and mean-subtract; this assumes the powers are 
        // chi^2_2 distributed, var(pow)=4, mean(pow)=std(pow)=2.  Therefore
        // dividing by mean will give standard deviation of 1 centred at 1.
        fft_out[j].x = pow/bp_l-1;
      }
    }

  }
  // write out current bandpass
  bp[i] = bp_l;
}

// sum polarizations in place
__global__ void pscrunch (cufftComplex *fft_out, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    //fft_out[i].x += fft_out[i+n].x;
    fft_out[i].x = M_SQRT1_2*(fft_out[i].x + fft_out[i+n].x);
  }
}

// sum polarizations in place
__global__ void pscrunch_weights (cufftComplex *fft_out, cufftReal* kur_weights_dev, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    // this formulation excludes samples with more than 80% RFI
    float w1_f = kur_weights_dev[i/NCHAN];
    float w2_f = kur_weights_dev[(i+n)/NCHAN];
    int w1 = w1_f >= MIN_WEIGHT;
    int w2 = w2_f >= MIN_WEIGHT;
    switch (w1+w2)
    {
      case 2:
        // both samples OK, account for variance with sqrt(2)
        fft_out[i].x = M_SQRT1_2*(fft_out[i].x + fft_out[i+n].x);
        kur_weights_dev[i/NCHAN] = 0.5*(w1_f + w2_f);
        break;
      case 1:
        // only one sample OK, variance = 1
        fft_out[i].x = w1*fft_out[i].x + w2*fft_out[i+n].x;
        //kur_weights_dev[i/NCHAN] = 0.5*(w1_f*w1 + w2_f*w2);
        kur_weights_dev[i/NCHAN] = w1_f*w1 + w2_f*w2;
        break;
      case 0:
        // no good samples, average bandpass (NB isn't this just 0?)
        //fft_out[i].x = 0.5*(fft_out[i].x + fft_out[i+n].x);
        fft_out[i].x = 0.;
        kur_weights_dev[i/NCHAN] = 0;
        break;
    }
  }
}

// average time samples
// TODO -- review normalization and make sure it's correct with polarization
__global__ void tscrunch (cufftComplex *fft_out, cufftReal* fft_ave,size_t n)
{
  // loop over the output indices; calculate corresponding input index,
  // then add up the subsequent NSCRUNCH samples
  float scale = sqrt (1./NSCRUNCH);
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    // explicit calculation of indices for future reference
    //int out_time_idx = i/NCHAN;
    //int out_chan_idx = i - out_time_idx*NCHAN;
    //int src_idx = out_time_idx * NSCRUNCH* NCHAN + out_chan_idx;
    //int src_idx = i+(NSCRUNCH-1)*out_time_idx*NCHAN;
    int src_idx = i+(NSCRUNCH-1)*(i/NCHAN)*NCHAN;
    fft_ave[i] = 0.;
    for (int j=0; j < NSCRUNCH; ++j, src_idx += NCHAN)
    {
      fft_ave[i] += fft_out[src_idx].x;
    }
    fft_ave[i] *= scale;
  }
}

__global__ void tscrunch_weights (cufftComplex *fft_out, cufftReal* fft_ave, cufftReal* kur_weights_dev, size_t n)
{
  // loop over the output indices; calculate corresponding input index,
  // then add up the subsequent NSCRUNCH samples
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {

    // explicit calculation of indices for future reference
    // int out_time_idx = i/NCHAN;
    // int out_chan_idx = i - out_time_idx*NCHAN;
    // int src_idx = out_time_idx * NSCRUNCH* NCHAN + out_chan_idx;
    // int src_idx = i+(NSCRUNCH-1)*out_time_idx*NCHAN;

    // TODO -- we might not want an additional cut on MIN_WEIGHT here
    int src_idx = i+(NSCRUNCH-1)*(i/NCHAN)*NCHAN;
    fft_ave[i] = 0.;
    int wt_sum = 0;
    float wt_sumf = 0;
    for (int j=0; j < NSCRUNCH; ++j, src_idx += NCHAN)
    {
      float wt = kur_weights_dev[src_idx/NCHAN];
      if (wt < MIN_WEIGHT) continue;
      wt_sum++;
      wt_sumf += wt;
      fft_ave[i] += wt*fft_out[src_idx].x;
    }
    if (wt_sumf/NSCRUNCH >= MIN_WEIGHT)
      fft_ave[i] /= sqrt(float(wt_sum));
    else
      // this just copies the bandpass in; NB the average is needed, I'm not
      // entirely sure why
      //fft_ave[i] = fft_out[src_idx].x/NSCRUNCH;
      fft_ave[i] = 0;
  }
}

// select and digitize
__global__ void sel_and_dig_2b (
    cufftReal *fft_ave, unsigned char* fft_trim_u, size_t n, int npol)
{
  int NCHANOUT = CHANMAX-CHANMIN+1;
  //int NTIME = (n*4) / (NCHANOUT*npol); // total time samples
  int NTIME = (VLITE_RATE/SEG_PER_SEC/NFFT)/NSCRUNCH;
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    // compute index into input array
    // correct for packing of 4 (because i indexes a *byte*, not a sample)
    //int time_idx = (i*4)/(NCHANOUT);
    //int chan_idx = i*4 - time_idx*NCHANOUT;
    int time_idx = (i*4)/(NCHANOUT*npol);
    int pol_idx = (i*4 - time_idx*NCHANOUT*npol)/NCHANOUT;
    int chan_idx = i*4 - time_idx*npol*NCHANOUT - pol_idx*NCHANOUT;
    fft_trim_u[i] = 0;
    for (int j = 0; j < 4; ++j)
    {
      // from Table 3 of Jenet & Anderson 1998
      //float tmp = fft_ave[time_idx*NCHAN+chan_idx+CHANMIN+j]/0.9674 + 1.5;
      float tmp = fft_ave[pol_idx*NTIME*NCHAN + time_idx*NCHAN+chan_idx+CHANMIN+j]/0.9675 + 1.5;
      if (tmp < 1) // do nothing, bit already correctly set
        continue;
      if (tmp < 2)
        fft_trim_u[i] += 1 << 2*j;
      else if (tmp < 3)
        fft_trim_u[i] += 2 << 2*j;
      else
        fft_trim_u[i] += 3 << 2*j;
    }
  }
}

// select and digitize
__global__ void sel_and_dig_4b (
    cufftReal *fft_ave, unsigned char* fft_trim_u, size_t n, int npol)
{
  int NCHANOUT = CHANMAX-CHANMIN+1;
  //int NTIME = n / (NCHANOUT*npol); // total time samples
  int NTIME = (VLITE_RATE/SEG_PER_SEC/NFFT)/NSCRUNCH;
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    // compute index into input array
    // correct for packing of 2 (because i indexes a *byte*, not a sample)
    //int time_idx = (i*2)/(NCHANOUT);
    //int chan_idx = i*2 - time_idx*NCHANOUT;
    int time_idx = (i*2)/(NCHANOUT*npol);
    int pol_idx = (i*2 - time_idx*NCHANOUT*npol)/NCHANOUT;
    int chan_idx = i*2 - time_idx*npol*NCHANOUT - pol_idx*NCHANOUT;
    // from Table 3 of Jenet & Anderson 1998
    //float tmp = fft_ave[time_idx*NCHAN+chan_idx+CHANMIN]/0.3188 + 7.5;
    float tmp = fft_ave[pol_idx*NTIME*NCHAN + time_idx*NCHAN+chan_idx+CHANMIN]/0.3188 + 7.5;
    if (tmp <= 0)
      fft_trim_u[i] = 0;
    else if (tmp >= 15)
      fft_trim_u[i] = 15;
    else
      fft_trim_u[i] = (unsigned char)(tmp);
    //tmp = fft_ave[time_idx*NCHAN+chan_idx+CHANMIN+1]/0.3188 + 7.5;
    tmp = fft_ave[pol_idx*NTIME*NCHAN + time_idx*NCHAN+chan_idx+CHANMIN+1]/0.3188 + 7.5;
    if (tmp <= 0)
      ;
    else if (tmp >= 15)
      fft_trim_u[i] += 15 << 4;
    else
      fft_trim_u[i] += (unsigned char)(tmp) << 4;
  }
}

// select and digitize
__global__ void sel_and_dig_8b (
    cufftReal *fft_ave, unsigned char* fft_trim_u, size_t n, int npol)
{
  int NCHANOUT = CHANMAX-CHANMIN+1;
  //int NTIME = n / (NCHANOUT*npol); // total time samples
  int NTIME = (VLITE_RATE/SEG_PER_SEC/NFFT)/NSCRUNCH;
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    // compute index into input array
    int time_idx = i/(NCHANOUT*npol);
    int pol_idx = (i - time_idx*NCHANOUT*npol)/NCHANOUT;
    int chan_idx = i - time_idx*npol*NCHANOUT - pol_idx*NCHANOUT;
    // from Table 3 of Jenet & Anderson 1998
    float tmp = fft_ave[pol_idx*NTIME*NCHAN + time_idx*NCHAN+chan_idx+CHANMIN]/0.02957 + 127.5;
    if (tmp <= 0)
      fft_trim_u[i] = 0;
    else if (tmp >= 255)
      fft_trim_u[i] = 255;
    else
      fft_trim_u[i] = (unsigned char) tmp;
  }
}

/*
// convert floating point to integer
__global__ void digitizearray(cufftReal *fft_ave, unsigned char* fft_ave_u, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    float tmp = fft_ave[i]/0.02957 + 127.5;
    if (tmp <= 0)
      fft_ave_u[i] = 0;
    else if (tmp >= 255)
      fft_ave_u[i] = 255;
    else
      fft_ave_u[i] = (unsigned char) tmp;
  }
}

// remove extraneous channels; in practice, this means 
__global__ void selectchannels(unsigned char* fft_ave_u, unsigned char* fft_trim_u, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    int nchan = CHANMAX-CHANMIN+1;
    int time_idx = i/nchan;
    int chan_idx = i - time_idx*nchan;
    fft_trim_u[i] = fft_ave_u[time_idx*NCHAN+chan_idx+CHANMIN];
  }
}

// detect total power and normalize the polarizations in place
// use a monolithic kernel pattern here since the total number of threads
// should be relatively small
// TODO -- this will probably be more efficient if done using lots of 
// threads and syncing; however, it doesn't seem to be a bottleneck
__global__ void detect_and_normalize (cufftComplex *fft_out, size_t ntime)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x; 
  if (i >= NCHAN*2) return;
  if (i >= NCHAN) // advance pointer to next polarization
  {
    fft_out += ntime*NCHAN;
    i -= NCHAN;
  }
  float sum1 = 0;
  float sum2 = 0;
  for (int j = i; j < ntime*NCHAN; j+= NCHAN)
  {
    float pow = fft_out[j].x*fft_out[j].x + fft_out[j].y*fft_out[j].y;
    fft_out[j].x = pow;
    sum1 += pow;
    sum2 += pow*pow;
  }
  sum1 *= 1./ntime;
  sum2 = sqrt(1./(sum2/ntime-sum1*sum1));
  for (int j = i; j < ntime*NCHAN; j+= NCHAN)
  {
    fft_out[j].x = (fft_out[j].x-sum1)*sum2;
  }
}

*/
