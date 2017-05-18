#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fcntl.h>

#include "vdifio.h"

#include <cufft.h>

#include "dada_def.h"
#include "dada_hdu.h"
#include "ipcio.h"
#include "ascii_header.h"
#include "multilog.h"

#include "util.h"

// from Julia's code
extern "C" {
#include "utils.h"
#include "def.h"
#include "executor.h"
}

// GPU options
#define DEVICE 0
#define NTHREAD 512

#define VD_FRM 5032
#define VD_DAT 5000
#define VLITE_RATE 128000000
#define VLITE_FRAME_RATE 25600
#define NFFT 12500 // number of samples to filterbank
#define NCHAN 6251 //output filterbank channels, including DC
#define NSCRUNCH 8 // time scrunch factor
#define CHANMIN 2411 // minimum output channel (counting DC), from low
#define CHANMAX 6250 // maximum output channel (counting DC), from low

//#define PROFILE 1

static volatile int NBIT = 2;

__global__ void convertarray (cufftReal *, unsigned char *, size_t);
__global__ void normalize (cufftComplex *, size_t);
__global__ void pscrunch (cufftComplex *, size_t n);
__global__ void tscrunch (cufftComplex *, cufftReal *, size_t);
__global__ void sel_and_dig_2b (cufftReal*, unsigned char*, size_t n);
__global__ void sel_and_dig_4b (cufftReal*, unsigned char*, size_t n);
__global__ void sel_and_dig_8b (cufftReal*, unsigned char*, size_t n);

void usage ()
{
  fprintf(stdout,"Usage: process [options]\n"
	  "-k hexadecimal shared memory key for input (default: 40)\n"
	  "-K hexadecimal shared memory key for output (default: 0=disabled)\n"
	  "-p listening port number (default: %zu)\n"
	  "-o print logging messages to stderr (as well as logfile)\n"
	  "-w output filterbank data (0=no sources, 1=listed sources, 2=all sources [def])\n"
	  "-b reduce output to b bits (2 [def], 4, and 8 are supported))\n"
	  ,(uint64_t)READER_SERVICE_PORT);
}


void print_cuda_properties () {
  cudaDeviceProp prop;
  int device_idx=DEVICE;
  cudaGetDeviceProperties (&prop,device_idx);
  printf("%s:\n",prop.name);
  printf("maxThreadsPerBlock= %d\n",prop.maxThreadsPerBlock);
  printf("maxThreadsDim[0]= %d\n",prop.maxThreadsDim[0]);
  printf("maxThreadsDim[1]= %d\n",prop.maxThreadsDim[1]);
  printf("maxThreadsDim[2]= %d\n",prop.maxThreadsDim[2]);
  printf("maxGridSize[0]= %d\n",prop.maxGridSize[0]);
  printf("maxGridSize[1]= %d\n",prop.maxGridSize[1]);
  printf("maxGridSize[2]= %d\n",prop.maxGridSize[2]);
}

int write_psrdada_header (dada_hdu_t* hdu, char* inchdr)
{

  // handle values from incoming header; set defaults if not available
  short int station_id = 0;
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
  int nchan = CHANMAX-CHANMIN + 1;
  double tsamp = double(NFFT)/VLITE_RATE*NSCRUNCH*1e6; // NB in mus
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
  dadacheck (ascii_header_set (ascii_hdr, "NPOL", "%d", 1) );
  dadacheck (ascii_header_set (ascii_hdr, "NBIT", "%d", NBIT) );
  dadacheck (ascii_header_set (ascii_hdr, "NPOL", "%d", NBIT) );
  dadacheck (ascii_header_set (ascii_hdr, "TSAMP", "%lf", tsamp) );
  dadacheck (ascii_header_set (ascii_hdr, "UTC_START", "%s", dada_utc) );
  multilog (hdu->log, LOG_INFO, "%s",ascii_hdr);
  ipcbuf_mark_filled (hdu->header_block, 4096);
  return 0;
}

int write_sigproc_header (FILE* output_fp, char* inchdr, vdif_header* vdhdr)
{
  // handle values from incoming header; set defaults if not available
  short int station_id = 0;
  ascii_header_get (inchdr, "STATIONID", "%d", &station_id);
  double ra = 0;
  ascii_header_get (inchdr, "RA", "%lf", &ra);
  double dec = 0;
  ascii_header_get (inchdr, "DEC", "%lf", &dec);
  char name[OBSERVATION_NAME_SIZE];
  ascii_header_get (inchdr, "NAME", "%s", name);
  name[OBSERVATION_NAME_SIZE-1] = '\0';

  double chbw = -64./NCHAN;
  int nchan = CHANMAX-CHANMIN + 1;
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
  send_double ("src_decj",sigproc_dec,output_fp);
  send_int ("data_type",1,output_fp);
  send_double ("fch1",384+(CHANMIN-0.5)*chbw,output_fp);
  send_double ("foff",chbw,output_fp);//negative foff, fch1 is highest freq
  send_int ("nchans",nchan,output_fp);
  send_int ("nbits",NBIT,output_fp);
  send_double ("tstart",mjd,output_fp);
  send_double ("tsamp",tsamp,output_fp);//[sec]
  send_int ("nifs",1,output_fp);
  send_string ("HEADER_END",output_fp);
  return 0;
}


void get_fbfile (char* fbfile, char* inchdr, vdif_header* vdhdr)
{
  // Open up filterbank file using timestamp and antenna
  short int station_id = 0;
  ascii_header_get (inchdr, "STATIONID", "%d", &station_id);
  char currt_string[128];
  struct tm tm_epoch = {0};
  int vdif_epoch = getVDIFEpoch (vdhdr);
  tm_epoch.tm_year = 100 + vdif_epoch/2;
  tm_epoch.tm_mon = 6*(vdif_epoch%2);
  time_t epoch_seconds = mktime (&tm_epoch) + getVDIFFrameEpochSecOffset (vdhdr);
  struct tm* utc_time = gmtime (&epoch_seconds);
  strftime(currt_string,sizeof(currt_string), "%Y%m%d_%H%M%S", utc_time);
  *(currt_string+15) = 0;
  snprintf (fbfile,128,"%s/%s_ea%02d.fil",DATADIR,currt_string,station_id);
}

int main (int argc, char *argv[])
{
  int exit_status = EXIT_SUCCESS;
  key_t key_in = 0x40;
  key_t key_out = 0x0;
  uint64_t port = READER_SERVICE_PORT;
  int stderr_output = 0;
  int write_fb = 2;

  int arg = 0;
  while ((arg = getopt(argc, argv, "hk:K:p:ow:b:")) != -1) {

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
      stderr_output = 1;
      break;

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

    }
  }

  cudacheck (cudaSetDevice (DEVICE));

  cudaDeviceProp prop;
  int device_idx=DEVICE;
  cudaGetDeviceProperties (&prop,device_idx);
  //print_cuda_properties();
  int nsms;
  cudaDeviceGetAttribute (&nsms,cudaDevAttrMultiProcessorCount,DEVICE);

  struct timespec ts_1ms;
  ts_1ms.tv_sec = 0;
  ts_1ms.tv_nsec = 1000000;

  //create and start timing
#ifdef PROFILE
  cudaEvent_t start,stop,start_total,stop_total;
  cudaEventCreate (&start);
  cudaEventCreate (&stop);
  cudaEventCreate (&start_total);
  cudaEventCreate (&stop_total);
  cudaEventRecord (start_total,0);
  float alloc_time=0, hdr_time=0, read_time=0, todev_time=0, convert_time=0, fft_time=0, normalize_time=0, tscrunch_time=0, pscrunch_time=0, digitize_time=0, write_time=0, flush_time=0, misc_time=0,elapsed=0;
#endif

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
      "%s/%s_%s_process_%06d.log",LOGDIR,hostname,currt_string,pid);
  FILE *logfile_fp = fopen (logfile, "w");
  multilog_add (log, logfile_fp);
  if (stderr_output)
    multilog_add (log, stderr);

  // connect to input buffer
  dada_hdu_t* hdu_in = dada_hdu_create (log);
  dada_hdu_set_key (hdu_in,key_in);
  if (dada_hdu_connect (hdu_in) != 0) {
    multilog (log, LOG_ERR, "Unable to connect to incoming PSRDADA buffer!\n");
    exit (EXIT_FAILURE);
  }
  char incoming_hdr[4096];

  // connect to output buffer (optional)
  dada_hdu_t* hdu_out = NULL;
  if (key_out) {
    // connect to the output buffer
    hdu_out = dada_hdu_create (log);
    dada_hdu_set_key (hdu_out,key_out);
    if (dada_hdu_connect (hdu_out) != 0) {
      multilog (log, LOG_ERR, "Unable to connect to outgoing PSRDADA buffer!\n");
      exit (EXIT_FAILURE);
    }
  }

#ifdef PROFILE
  cudaEventRecord(start,0);
#endif
  
  // allocate memory for 1 second; in the future, consider keeping
  // multiple buffers to allow for out-of-order data; this implementation
  // will handle dropped data
  unsigned char* udat; cudacheck (
  cudaMallocHost ((void**)&udat, 2*VLITE_RATE) );
  // device memory for 8-bit samples
  unsigned char* udat_dev; cudacheck (
  cudaMalloc ((void**)&udat_dev,2*VLITE_RATE/10) );

  // set up memory for FFTs
  // device memory for 32-bit floats
  int numfft = 2*VLITE_RATE/10/NFFT;
  cufftHandle plan;
  cufftcheck (cufftPlan1d (&plan,NFFT,CUFFT_R2C,numfft));

  cufftReal* fft_in; cudacheck (
  cudaMalloc ((void**)&fft_in,sizeof(cufftReal)*numfft*NFFT) );

  cufftComplex* fft_out; cudacheck (
  cudaMalloc ((void**)&fft_out,sizeof(cufftComplex)*numfft*NCHAN) );

  // NB, reduce by a further factor of 2 since pols summed in this buffer
  int scrunch = (numfft*NCHAN)/(2*NSCRUNCH);
  cufftReal* fft_ave; cudacheck (
  cudaMalloc ((void**)&fft_ave,sizeof(cufftReal)*scrunch) );

  // error check that NBIT is commensurate with trimmed array size
  int trim = (numfft*(CHANMAX-CHANMIN+1))/(2*NSCRUNCH);
  if (trim % (8/NBIT) != 0)
  {
    multilog (log, LOG_ERR,
      "Selected channel and bit scheme is not commensurate!.\n");
    exit(EXIT_FAILURE);
  }
  // reduce array size by packing of samples into byte
  trim /= (8/NBIT);
  unsigned char* fft_trim_u; cudacheck (
  cudaMalloc ((void**)&fft_trim_u,trim) );

  unsigned char* fft_trim_u_host; cudacheck (
  cudaMallocHost ((void**)&fft_trim_u_host,trim) );

#ifdef PROFILE
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&alloc_time,start,stop);
#endif

  // allocate a frame buffer on the stack
  char frame_buff[5032];
  vdif_header* vdhdr = (vdif_header*) frame_buff;
  uint8_t* dat = (uint8_t*)(frame_buff + 32);

  // connect to control socket
  Connection c;
  c.sockoptval = 1; //release port immediately after closing connection
  if (serve (port, &c) < 0) {
    multilog (log, LOG_ERR,
      "Failed to create control socket on port %d.\n", port);
    exit (EXIT_FAILURE);
  }
  fcntl (c.rqst, F_SETFL, O_NONBLOCK); // set up for polling
  char cmd_buff[32];

  int quit = 0;

  // This point start to loop over observations.
  while(true)
  {

  if (quit)
    break;

  multilog (log, LOG_INFO, "Waiting for DADA header.\n");
  dadacheck (dada_hdu_lock_read (hdu_in) );

  // this will block until a header has been created by writer, signalling the start of
  // an observation
  uint64_t hdr_size = 0;
  char* ascii_hdr = ipcbuf_get_next_read (hdu_in->header_block,&hdr_size);

  if (ascii_hdr == NULL) {
    // this could be an error condition, or it could be a result of shutting down the
    // data acquisition system; check to see if a CMD_QUIT has been issued in order to
    // log the appropriate outcome
    ssize_t npoll_bytes = read (c.rqst, cmd_buff, 32);
    if (npoll_bytes >  0) {
      for (int ib = 0; ib < npoll_bytes; ib++) {
        printf("Read command character %c.\n",cmd_buff[ib]);
        if (cmd_buff[ib] == CMD_QUIT) {
          multilog (log, LOG_INFO,
              "Received CMD_QUIT, indicating data taking is ceasing.  Exiting.\n");
          break;
        }
      }
    }
    // did not receive CMD_QUIT, so this is an abnormal state
    multilog (log, LOG_ERR, "PSRDADA read failed unexpectedly.  Terminating.\n");
    exit_status = EXIT_FAILURE;
    break;
  }

  multilog (log, LOG_INFO, "Beginning new observation.\n\n");
  memcpy (incoming_hdr,ascii_hdr,hdr_size);
  dadacheck (ipcbuf_mark_cleared (hdu_in->header_block));

  int bad_packets = 0;
  int skip_observation = 1;
  // consume from the buffer until we reach 1-s boundary
  for (int i=0; i < VLITE_FRAME_RATE*2; ++i)
  {

    ssize_t nread = ipcio_read(hdu_in->data_block,frame_buff,VD_FRM);

    if (nread != VD_FRM)
    {
      bad_packets++;
      // this will skip the short data segment and move to the next "observation"
      // TODO -- not sure why this sleep is in here, do we need it?
      nanosleep (&ts_1ms,NULL);
      continue;
    }
    if (i==0)
    {
      int sample = getVDIFFrameNumber(vdhdr);
      int current_sec = getVDIFFrameSecond(vdhdr);
      multilog (log, LOG_INFO, 
        "Starting trim at frame %d and second %d.\n",sample,current_sec);
    }
    if (i==(2*VLITE_FRAME_RATE-1))
    {
      int sample = getVDIFFrameNumber(vdhdr);
      int current_sec = getVDIFFrameSecond(vdhdr);
      multilog (log, LOG_INFO, 
        "Ending unsuccessfully at frame %d and second %d.\n",
        sample,current_sec);
    }
    if (nread < 0)
    {
      multilog (log, LOG_ERR, "Problem with nread! Skipping observation.\n");
      skip_observation = 1;
      break;
    }
    int sample = getVDIFFrameNumber(vdhdr);
    if (sample==0)
    {
      skip_observation = 0;
      break;
    }
  }

  if (bad_packets > 0)
    multilog (log, LOG_INFO,
      "Exiting buffer consumption with %d bad packets.\n",bad_packets);

  if (skip_observation)
  {
    multilog (log, LOG_INFO, "Skipping observation.\n");
    dadacheck (dada_hdu_unlock_read (hdu_in));
    continue;
  }

  // set up "previous" samples and copy first frame into 1s-buffer
  int current_sec = getVDIFFrameSecond(vdhdr);
  int current_sample[2] = {-1,-1};
  multilog (log, LOG_INFO, "Starting sec=%d\n",current_sec);

  // Open up filterbank file at appropriate time
  char fbfile[128];
  get_fbfile (fbfile, incoming_hdr, vdhdr);
  if (write_fb == 0) {
    multilog (log, LOG_INFO,
        "Filterbank output disabled.  Would have written to %s.\n",fbfile);
    snprintf (fbfile,128,"/dev/null"); 
  }
  if (write_fb == 1) {
    char name[OBSERVATION_NAME_SIZE];
    ascii_header_get (incoming_hdr, "NAME", "%s", name);
    name[OBSERVATION_NAME_SIZE-1] = '\0';
    if (check_name (name)) {
      multilog (log, LOG_INFO,
          "Source %s matches target list, recording filterbank data.\n", name);
    }
    else {
      multilog (log, LOG_INFO,
          "Source %s not on target list, disabling filterbank data.\n", name);
      multilog (log, LOG_INFO,
          "Filterbank output disabled.  Would have written to %s.\n",fbfile);
      snprintf (fbfile,128,"/dev/null"); 
    }
  }

  FILE *fb_fp = myopen (fbfile, "wb", true);
  multilog (log, LOG_INFO, "Writing data to %s.\n",fbfile);
  uint64_t fb_bytes_written = 0;

#ifdef PROFILE
  cudaEventRecord(start,0);
#endif

  if (key_out) {
    //Write psrdada output header values
    write_psrdada_header (hdu_out, incoming_hdr);
  }

  // write out a sigproc header
  write_sigproc_header (fb_fp, incoming_hdr, vdhdr);

#ifdef PROFILE
  cudaEventRecord (stop,0);
  cudaEventSynchronize (stop);
  cudaEventElapsedTime (&hdr_time, start,stop);
#endif

  while(true) // loop over data packets
  {
    int thread = getVDIFThreadID (vdhdr) != 0;
    int sample = getVDIFFrameNumber (vdhdr);
    int second = getVDIFFrameSecond (vdhdr);

    if (second == current_sec)
    {

#ifdef PROFILE
      cudaEventRecord (start,0);
#endif 

      size_t idx = VLITE_RATE*thread + sample*VD_DAT;
      memcpy (udat + idx,dat,VD_DAT);

      // check for skipped data
      if (sample!=current_sample[thread] + 1)
      {
        multilog (log, LOG_INFO, "Found skipped data in current second!\n");
        // set the intervening memory to zeros
        idx = VLITE_RATE*thread + (current_sample[thread]+1)*VD_DAT;
        memset (udat+idx,0,VD_DAT*(sample-current_sample[thread]-1));
      }
      current_sample[thread] = sample;

      // read the next frame into the buffer
      ssize_t nread = ipcio_read (hdu_in->data_block,frame_buff,VD_FRM);

#ifdef PROFILE
      cudaEventRecord (stop,0);
      cudaEventSynchronize (stop);
      cudaEventElapsedTime (&elapsed,start,stop);
      read_time += elapsed;
#endif

      if (nread < 0 ) {
        multilog (log, LOG_ERR, "Error on nread.\n");
        return (EXIT_FAILURE);
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
        "Major data skip!  (%d vs. %d) Aborting this observation.\n",
         second,current_sec);
         break;
    }

    // every 1s, check for a QUIT command
    int npoll_bytes = read (c.rqst, cmd_buff, 32);
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

    // check that both buffers are full
    current_sec = second;

    // need to dispatch the current buffer; zero fill any discontinuity
    // at the end
    for (int i = 0; i < 2; ++i)
    {
      int target = VLITE_FRAME_RATE - 1;
      if (current_sample[i] !=  target)
      {
        multilog (log, LOG_ERR,
          "Found skipped data in dispatch!  Aborting this observation.\n");
        //size_t idx = VLITE_RATE*i + current_sample[i]*VD_DAT;
        //memset(udat+idx,0,(target-current_sample[i])*VD_DAT);
        break;
      }
      current_sample[i] = -1;
    }

    // do dispatch -- break into 0.1s chunks to fit in GPU memory
    size_t nchunk = numfft*NFFT; // samples in a chunk, 2*VLITE_RATE/10
    for (int i = 0; i < 10; ++i)
    {
#ifdef PROFILE
      cudaEventRecord(start,0);
#endif 
      // need to copy pols separately
      for (int ipol=0; ipol < 2; ++ipol)
      {
        unsigned char* dst = udat_dev + ipol*nchunk/2;
        unsigned char* src = udat + i*nchunk/2;
        cudacheck (cudaMemcpy (dst,src,nchunk/2,cudaMemcpyHostToDevice) );
      }
#ifdef PROFILE
      cudaEventRecord(stop,0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&elapsed,start,stop);
      todev_time += elapsed;
#endif

#ifdef PROFILE
      cudaEventRecord(start,0);
#endif 
      convertarray<<<nsms*32,NTHREAD>>> (fft_in,udat_dev,nchunk);
      cudacheck (cudaGetLastError() );
#ifdef PROFILE
      cudaEventRecord(stop,0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&elapsed,start,stop);
      convert_time += elapsed;
#endif

#ifdef PROFILE
      cudaEventRecord(start,0);
#endif 
      cufftcheck (cufftExecR2C (plan,fft_in,fft_out) );
#ifdef PROFILE
      cudaEventRecord(stop,0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&elapsed,start,stop);
      fft_time += elapsed;
#endif

#ifdef PROFILE
      cudaEventRecord(start,0);
#endif 
      normalize<<<(NCHAN*2)/NTHREAD+1,NTHREAD>>> (fft_out,numfft/2);
      cudacheck ( cudaGetLastError() );
#ifdef PROFILE
      cudaEventRecord(stop,0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&elapsed,start,stop);
      normalize_time += elapsed;
#endif

#ifdef PROFILE
      cudaEventRecord(start,0);
#endif 
      size_t maxn = (numfft*NCHAN)/2; // reduce by 2 for polarization
      pscrunch<<<nsms*32,NTHREAD>>> (fft_out,maxn);
      cudacheck ( cudaGetLastError() );
#ifdef PROFILE
      cudaEventRecord(stop,0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&elapsed,start,stop);
      pscrunch_time += elapsed;
#endif

#ifdef PROFILE
      cudaEventRecord(start,0);
#endif 
      maxn /= NSCRUNCH;
      tscrunch<<<nsms*32,NTHREAD>>> (fft_out,fft_ave,maxn);
      cudacheck ( cudaGetLastError() );
#ifdef PROFILE
      cudaEventRecord(stop,0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&elapsed,start,stop);
      tscrunch_time += elapsed;
#endif

#ifdef PROFILE
      cudaEventRecord(start,0);
#endif 
      maxn = (CHANMAX-CHANMIN+1)*(maxn/NCHAN)/(8/NBIT);
      switch (NBIT)
      {
        case 2:
          sel_and_dig_2b<<<nsms*32,NTHREAD>>> (fft_ave,fft_trim_u,maxn);
          break;
        case 4:
          sel_and_dig_4b<<<nsms*32,NTHREAD>>> (fft_ave,fft_trim_u,maxn);
          break;
        case 8:
          sel_and_dig_8b<<<nsms*32,NTHREAD>>> (fft_ave,fft_trim_u,maxn);
          break;
        default:
          sel_and_dig_2b<<<nsms*32,NTHREAD>>> (fft_ave,fft_trim_u,maxn);
          break;
      }
      cudacheck ( cudaGetLastError() );
#ifdef PROFILE
      cudaEventRecord(stop,0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&elapsed,start,stop);
      digitize_time += elapsed;
#endif

#ifdef PROFILE
      cudaEventRecord(start,0);
#endif 

      // copy filterbanked data back to host
      cudacheck (cudaMemcpy (
            fft_trim_u_host,fft_trim_u,maxn,cudaMemcpyDeviceToHost) );

#ifdef PROFILE
      cudaEventRecord(stop,0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&elapsed,start,stop);
      misc_time += elapsed;
#endif

      // finally, push the filterbanked time samples onto psrdada buffer
      // and/or write out to sigproc

#ifdef PROFILE
      cudaEventRecord(start,0);
#endif 
      if (key_out)
        ipcio_write (hdu_out->data_block,(char *)fft_trim_u_host,maxn);
      fwrite(fft_trim_u_host,1,maxn,fb_fp);
      fb_bytes_written += maxn;
#ifdef PROFILE
      cudaEventRecord(stop,0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&elapsed,start,stop);
      write_time += elapsed;
#endif

    }
  }

#ifdef PROFILE
  cudaEventRecord(start,0);
#endif 

  dadacheck (dada_hdu_unlock_read (hdu_in));
  if (key_out)
    dadacheck (dada_hdu_unlock_write (hdu_out));

  // close files
  fclose(fb_fp);
  multilog (log, LOG_INFO, "Wrote %.2f MB to %s\n",fb_bytes_written*1e-6,fbfile);

#ifdef PROFILE
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&flush_time,start,stop);

  cudaEventRecord(stop_total,0);
  cudaEventSynchronize(stop_total);
  float total_time;
  cudaEventElapsedTime(&total_time,start_total,stop_total);
  float sub_time = alloc_time + hdr_time + read_time + todev_time + convert_time + fft_time + normalize_time + pscrunch_time + tscrunch_time + digitize_time + write_time + flush_time + misc_time;
  multilog (log, LOG_INFO, "Alloc Time..%.3f\n", alloc_time*1e-3);
  multilog (log, LOG_INFO, "Read Time...%.3f\n", read_time*1e-3);
  multilog (log, LOG_INFO, "Copy To Dev.%.3f\n", todev_time*1e-3);
  multilog (log, LOG_INFO, "Convert.....%.3f\n", convert_time*1e-3);
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

  } // end loop over observations

  // close sockets and files
  shutdown (c.rqst, 2);
  fclose (logfile_fp);

  // clean up psrdada data structures
  dada_hdu_disconnect (hdu_in);
  dada_hdu_destroy (hdu_in);
  if (key_out) {
    dada_hdu_disconnect (hdu_out);
    dada_hdu_destroy (hdu_out);
  }

  //free memory
  cudaFreeHost(udat);
  cudaFree(udat_dev);
  cudaFree(fft_in);
  cudaFree(fft_out);
  cudaFree(fft_ave);
  cudaFree(fft_trim_u);
  cudaFreeHost(fft_trim_u_host);

  return exit_status;
    
}


//convert unsigned char time array to float
__global__ void convertarray(cufftReal *time, unsigned char *utime, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    if (utime[i] == 0) 
      time[i] = 0;
    else
      time[i] = (cufftReal)(utime[i])/128-1;
  }
}

// detect total power and normalize the polarizations in place
// use a monolithic kernel pattern here since the total number of threads
// should be relatively small
__global__ void normalize(cufftComplex *fft_out, size_t ntime)
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

// detect and sum polarizations; do in place
__global__ void pscrunch(cufftComplex *fft_out, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    //cufftComplex pol1 = fft_out[i];
    //cufftComplex pol2 = fft_out[i+n];
    //fft_out[i].x = pol1.x*pol1.x + pol1.y*pol1.y + 
    //               pol2.x*pol2.x + pol2.y*pol2.y;
    fft_out[i].x += fft_out[i+n].x;
  }
}

// average time samples
__global__ void tscrunch(cufftComplex *fft_out, cufftReal* fft_ave,size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    // explicit calculation of indices for future references
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
    // additional factor of 2 here corrects for polarization sum
    // TODO precompute this or is compiler smart enough?
    // and added an extra factor of 2 to get output std. correct
    fft_ave[i] *= sqrt(0.25/NSCRUNCH);
  }
}

// select and digitize
__global__ void sel_and_dig_2b(
    cufftReal *fft_ave, unsigned char* fft_trim_u, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    // compute index into input array
    int nchan = CHANMAX-CHANMIN+1;
    // correct for packing of 4
    int time_idx = (i*4)/(nchan);
    int chan_idx = i*4 - time_idx*nchan;
    fft_trim_u[i] = 0;
    for (int j = 0; j < 4; ++j)
    {
      // from Table 3 of Jenet & Anderson 1998
      float tmp = fft_ave[time_idx*NCHAN+chan_idx+CHANMIN+j]/0.9674 + 1.5;
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
__global__ void sel_and_dig_4b(
    cufftReal *fft_ave, unsigned char* fft_trim_u, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    // compute index into input array
    int nchan = CHANMAX-CHANMIN+1;
    // correct for packing of 2
    int time_idx = (i*2)/(nchan);
    int chan_idx = i*2 - time_idx*nchan;
    // from Table 3 of Jenet & Anderson 1998
    float tmp = fft_ave[time_idx*NCHAN+chan_idx+CHANMIN]/0.3188 + 7.5;
    if (tmp <= 0)
      fft_trim_u[i] = 0;
    else if (tmp >= 15)
      fft_trim_u[i] = 15;
    else
      fft_trim_u[i] = (unsigned char)(tmp);
    tmp = fft_ave[time_idx*NCHAN+chan_idx+CHANMIN+1]/0.3188 + 7.5;
    if (tmp <= 0)
      ;
    else if (tmp >= 15)
      fft_trim_u[i] += 15 << 4;
    else
      fft_trim_u[i] += (unsigned char)(tmp) << 4;
  }
}

// select and digitize
__global__ void sel_and_dig_8b(
    cufftReal *fft_ave, unsigned char* fft_trim_u, size_t n)
{
  for (
      int i = threadIdx.x + blockIdx.x*blockDim.x; 
      i < n; 
      i += blockDim.x*gridDim.x)
  {
    // compute index into input array
    int nchan = CHANMAX-CHANMIN+1;
    int time_idx = i/(nchan);
    int chan_idx = i - time_idx*nchan;
    // from Table 3 of Jenet & Anderson 1998
    float tmp = fft_ave[time_idx*NCHAN+chan_idx+CHANMIN]/0.02957 + 127.5;
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
*/
