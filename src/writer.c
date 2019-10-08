#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <pthread.h>
#include <math.h>

#include "dada_def.h"
#include "dada_hdu.h"
#include "ipcio.h"
#include "ascii_header.h"
#include "multilog.h"

#include "def.h"
#include "utils.h"
#include "executor.h"
#include "vdifio.h"
#include "multicast.h"

#define MSGMAXSIZE 8192
#define OD_CACHE_SIZE 15

static FILE* logfile_fp = NULL;
static struct timespec ts_1ms;
static const char cmdstop[] = {CMD_STOP};
static int mc_control_sock = 0;
static int mc_info_sock = 0;
static multilog_t* mlog = NULL;
static dada_hdu_t* hdu = NULL;

void usage ()
{
  fprintf(stdout,"Usage: writer [options]\n"
	  "-k hexadecimal shared memory key  (default: 40)\n"
	  "-e ethernet device id (default: eth0)\n"
	  "-o print logging messages to stderr (as well as logfile)\n"
    );
}

void cleanup (void)
{
  fprintf (stderr,"called cleanup! [WRITER]\n");
  fflush (stderr);
  if (hdu != NULL)
  {
    dada_hdu_disconnect (hdu);
    dada_hdu_destroy (hdu);
  }
  fprintf (stderr,"h1w\n");
  fflush (stderr);
  if (mc_control_sock > 0)
    shutdown (mc_control_sock, 2);
  if (mc_info_sock > 0)
    shutdown (mc_info_sock, 2);
  fprintf (stderr,"h2w\n");
  fflush (stderr);
  multilog (mlog, LOG_INFO, "Completed shutdown [WRITER].\n");
  multilog_close (mlog);
  if (logfile_fp) fclose (logfile_fp);
  fprintf (stderr,"h3w\n");
  fflush (stderr);
}

void exit_handler (void)
{
  fprintf (stderr, "called exit_handler\n");
  fflush (stderr);
  cleanup ();
}

void sigint_handler (int dummy) {
  fclose (logfile_fp);
}

void dadacheck (int rcode)
{
  if (rcode < 0)
  {
    printf("dada problem %d",rcode);
    exit (EXIT_FAILURE);
  }
}

int write_psrdada_header(dada_hdu_t* hdu, vdif_header* hdr, ObservationDocument* od)
{
  char* ascii_hdr = ipcbuf_get_next_write (hdu->header_block);
  int station_id = getVDIFStationID (hdr);
  dadacheck (ascii_header_set (ascii_hdr, "STATIONID", "%d",station_id) );
  dadacheck (ascii_header_set (ascii_hdr, "NCHAN", "%d", 1) );
  dadacheck (ascii_header_set (ascii_hdr, "BANDWIDTH", "%lf", -64.) );
  dadacheck (ascii_header_set (ascii_hdr, "CFREQ", "%lf", 352.) );
  dadacheck (ascii_header_set (ascii_hdr, "NPOL", "%d", 2) );
  dadacheck (ascii_header_set (ascii_hdr, "NBIT", "%d", 8) );
  dadacheck (ascii_header_set (ascii_hdr, "TSAMP", "%lf", 0.0078125) );
  dadacheck (ascii_header_set (ascii_hdr, "RA", "%lf", od->ra) );
  dadacheck (ascii_header_set (ascii_hdr, "DEC", "%lf", od->dec) );
  dadacheck (ascii_header_set (ascii_hdr, "NAME", "%s", od->name) );
  dadacheck (ascii_header_set (ascii_hdr, "SCANSTART", "%lf", od->startTime) );
  dadacheck (ascii_header_set (ascii_hdr, "DATAID", "%s", od->datasetId) );
  // get time from VDIF frame
  struct tm tm_epoch = {0};
  int vdif_epoch = getVDIFEpoch (hdr);
  tm_epoch.tm_year = 100 + vdif_epoch/2;
  tm_epoch.tm_mon = 6*(vdif_epoch%2);
  time_t epoch_seconds = mktime (&tm_epoch) + getVDIFFrameEpochSecOffset (hdr);
  struct tm* utc_time = gmtime (&epoch_seconds);
  char dada_utc[64];
  strftime (dada_utc, 64, DADA_TIMESTR, utc_time);
  dadacheck (ascii_header_set (ascii_hdr, "UTC_START", "%s", dada_utc) );
  dadacheck (ascii_header_set (ascii_hdr, "UNIX_TIMET", "%ld", epoch_seconds) );
  multilog (hdu->log, LOG_INFO, "psrdada header:\n%s",ascii_hdr);
  ipcbuf_mark_filled (hdu->header_block, 4096);
  return 0;
}

void print_observation_document(char* result, ObservationDocument* od)
{
  int lens = 128;
  char s[lens];
  snprintf(result,2048,"    datasetId = %s\n", od->datasetId);
  snprintf(s,lens,"    configId = %s\n", od->configId);
  strcat(result,s);
  snprintf(s,lens,"    startTime = %10.8f\n", od->startTime);
  strcat(result,s);
  snprintf(s,lens,"    name = %s\n", od->name);
  strcat(result,s);
  snprintf(s,lens,"    ra = %10.8f\n", od->ra);
  strcat(result,s);
  snprintf(s,lens,"    dec = %10.8f\n", od->dec);
  strcat(result,s);
  snprintf(s,lens,"    dra = %10.8f\n", od->dra);
  strcat(result,s);
  snprintf(s,lens,"    ddec = %10.8f\n", od->ddec);
  strcat(result,s);
  snprintf(s,lens,"    azoffs = %10.8f\n", od->azoffs);
  strcat(result,s);
  snprintf(s,lens,"    eloffs = %10.8f\n", od->eloffs);
  strcat(result,s);
  snprintf(s,lens,"    startLST = %10.8f\n", od->startLST);
  strcat(result,s);
  snprintf(s,lens,"    scanNo = %d\n", od->scanNo);
  strcat(result,s);
  snprintf(s,lens,"    subscanNo = %d\n", od->subscanNo);
  strcat(result,s);
  snprintf(s,lens,"    primaryBand = %s\n", od->primaryBand);
  strcat(result,s);
  snprintf(s,lens,"    usesPband = %d\n", od->usesPband);
  strcat(result,s);
}

void fprint_observation_document(FILE* fd, ObservationDocument* od, int brief)
{
  // convert to hhmmss
  double fracday = od->startTime-(int)od->startTime;
  fracday *= 24;
  int hh = (int)(fracday);
  fracday -= hh;
  fracday *= 60;
  int mm = (int)(fracday);
  fracday -= mm;
  fracday *= 60;
  int ss = (int)(fracday+0.5);
  fprintf(fd,"    datasetId = %s\n", od->datasetId);
  fprintf(fd,"    configId = %s\n", od->configId);
  fprintf(fd,"    startTime = %10.8f [%02d:%02d:%02d]\n", od->startTime,
      hh, mm, ss);
  fprintf(fd,"    name = %s\n", od->name);
  fprintf(fd,"    ra = %10.8f, dec = %10.8f\n", od->dec, od->ra);
  fprintf(fd,"    startLST = %10.8f\n", od->startLST);
  fprintf(fd,"    scanNo = %d; subscanNo = %d\n", od->scanNo, od->subscanNo);
  //fprintf(fd,"    subscanNo = %d\n", od->subscanNo);
  fprintf(fd,"    primaryBand = %s\n", od->primaryBand);
  if (brief)
    return;
  fprintf(fd,"    dra = %10.8f\n", od->dra);
  fprintf(fd,"    ddec = %10.8f\n", od->ddec);
  fprintf(fd,"    azoffs = %10.8f\n", od->azoffs);
  fprintf(fd,"    eloffs = %10.8f\n", od->eloffs);
  fprintf(fd,"    usesPband = %d\n", od->usesPband);
}

/*
 * Fill the ObservationDocument pointed to with fake values.  If the scan
 * time is set correctly, this will initiate or terminate data taking.
 * To terminate, specify finish==True.
 */
void fake_observation_document (ObservationDocument* od, double scan_start,
  int finish, int change_ra)
{
  static double last_ra = 0;
  od->startTime = scan_start+1./86400; // in decimal/double MJD(UTC)
  if (finish)
		snprintf(od->name,OBSERVATION_NAME_SIZE,"FINISH");
  else
		snprintf(od->name,OBSERVATION_NAME_SIZE,"FAKE");
  if (change_ra)
  {
    od->ra = last_ra;
    last_ra += 0.1;
    if (last_ra > 3.14*2)
      last_ra = 0;
  }
  else
    od->ra = 1;
  od->dec = 1;
}

// TODO -- currently, we are only signalling stop by EOD.  I think
// this is good, but we need to ensure consistency with PB, such that
// the voltage buffer doesn't fill up if PB for some reason aborts
// early.  Potentially could be handle by a multicast message.  Or just
// quit and restart...

void stop_observation (dada_hdu_t* hdu, multilog_t* log, 
    uint64_t* packets_written, int* state)
{
  if (*state != STATE_STARTED)
    return;

  *state = STATE_STOPPED;
  if (dada_hdu_unlock_write (hdu) < 0) {
    multilog (log, LOG_ERR,
        "[STATE_STARTED->STOP]: unable to unlock psrdada HDU, exiting\n");
    exit (EXIT_FAILURE);
  }
  multilog (log, LOG_INFO,
      "[STATE_STARTED->STOP] Wrote %d packets to psrdada buffer.\n",*packets_written);
  *packets_written = 0;

  nanosleep (&ts_1ms, NULL);
}

/* NO longer needed
// signal other writers to stop observation, if started
int mc_stop_observation (dada_hdu_t* hdu, multilog_t* log, int* state)
{
  if (*state != STATE_STARTED)
    return 1;
  if (MulticastSend (mc_vlitegrp, MC_WRITER_PORT, cmdstop, 1) < 0)
  {
    multilog (log, LOG_ERR, "send: %s\n", strerror (errno));
    return 0;
  }
  return 1;
}
*/

/*
 * Return all buffers that have not yet been written out and that overlap
 * any of the trigger window.  By fiat, all buffers are now aligned to a
 * UTC second.
 */
int get_buffer_trigger_overlap (char** bufs, int nbufs,
    trigger_t** tqueue, time_t* dump_times,
    char** bufs_to_write, time_t* times_to_write)
{
  int nbuf_to_write = 0;
  double min_time = 1e100;
  double max_time = -1e100;
  for (int ibuf=0; ibuf < nbufs; ++ibuf)
  {
    // get time for frame at start of buffer
    // TODO -- we can actually cache the seconds when we do the 1-s
    // boundary, since we create alignment by design.  Leave for now.
    vdif_header* vdhdr = (vdif_header*) bufs[ibuf];
    time_t utc_seconds;
    double buf_t0 = vdif_to_dunixepoch (vdhdr, &utc_seconds); // beginning of buffer
    if (buf_t0 < min_time) min_time = buf_t0;
    if (buf_t0 > max_time) max_time = buf_t0;

    /* TMP -- leave for debugging purposes
    char dada_utc[32];
    struct tm utc_time;
    gmtime_r (&utc_seconds, &utc_time);
    strftime (dada_utc, 64, DADA_TIMESTR, &utc_time);
    if (ibuf==0)
      fprintVDIFHeader(stderr, vdhdr, VDIFHeaderPrintLevelColumns);
    multilog (mlog, LOG_INFO, "buffer UTC %s\n", dada_utc);
    printf ("utc_seconds = %ld\n",utc_seconds);
    fprintVDIFHeader (stderr, vdhdr, VDIFHeaderPrintLevelShort);
    */

    int trigger_overlap = 0;
    for (int itrig=0; itrig < MAX_TRIGGERS; ++itrig)
    {
      if (tqueue[itrig] == NULL)
        break;
      // condition for overlap is 
      // !((trigger_start > buffer_end) || (trigger_end < buffer_start)
      // --> (trigger_start < buffer_end && trigger_end > buffer_start)
      if ((tqueue[itrig]->t0 <= (buf_t0+1)) &&
          (tqueue[itrig]->t1 >= buf_t0))
      {
        trigger_overlap = 1;
        break;
      }
    } // end loop over trigger times

    // continue to next buffer if no trigger overlaps
    if (!trigger_overlap) continue;

    // check to make sure we haven't already dumped it
    int already_dumped = 0;
    for (int idump=0; idump < MAX_THREADIOS; ++idump)
    {
      if (utc_seconds == dump_times[idump])
      {
        already_dumped = 1;
        break;
      }
    } // end loop over dump times

    if (already_dumped) continue;

    times_to_write[nbuf_to_write] = utc_seconds;
    bufs_to_write[nbuf_to_write++] = bufs[ibuf];
  } // end loop over buffers
  fprintf (stderr,"earliest buffer t0 = %.3f\n",min_time);
  fprintf (stderr,"latest   buffer t0 = %.3f\n",max_time);
  return nbuf_to_write;
}

/*
 * Return an ObservationDocument from the cache if its startTime matches
 * the current second according to the VDIF stream.
 */
ObservationDocument* search_od_cache (ObservationDocument* first_od, int od_cache_size, double vdif_mjd)
{
  // round VDIF timestamp to nearest second in the day
  int vdif_sec = (int)(0.5+86400*(vdif_mjd-(int)vdif_mjd));
  for (int i = 0; i < od_cache_size; ++i, ++first_od) 
  {
    double od_mjd = first_od->startTime;
    int od_sec = (int)(0.5+86400*(od_mjd-(int)od_mjd));
    if (od_sec == vdif_sec)
      return first_od;
  }
  return NULL;
}

/*
 * Check to see if the pointing position is close enough between scans
 * to consider them contiguous.  This small amount allows for things
 * like VLASS.  For speed, we are also ignoring projection effects!
 */
int check_od_consistency (ObservationDocument* od1, 
    ObservationDocument* od2, multilog_t* log, int obs_length)
{

  // check for edge case of FINISH scan
  if (strcasecmp (od2->name,"FINISH") == 0)
    return 0;

  double rtol = 0.00873; // 0.5 deg
  int max_integ = 480; // 8 minutes
  if ((fabs(od1->ra-od2->ra)< rtol) && fabs(od1->dec-od2->dec) < rtol)
  {
     if (obs_length < max_integ)
     {
       multilog (log, LOG_INFO,
           "Pointing unchanged, will continue to integrate [t=%d].\n",obs_length);
       return 1;
     }
     multilog (log, LOG_INFO,
         "Integration exceeded %ds, starting a new one.\n", max_integ);
  }
  return 0;
}

/*
 * Return the difference, in frames, between two packets represented by
 * the provided VDIF headers.  It is assumed hdr2 will follow in time.
 * This will include the total number of frames, that is, for both threads/
 * polarizations in the data stream.  It will be 1 for contiguous data.
 */
int vdif_frame_difference (vdif_header* hdr1, vdif_header* hdr2)
{
  // recall our convention for thread, thread = threadid != 0
  // viz. threadid == 0 = 0, anything else mapped to 0
  return ((int)hdr2->seconds-hdr1->seconds)*(2*FRAMESPERSEC) +
         ((int)hdr2->frame-hdr1->frame)*2 + 
         ((int)(hdr2->threadid!=0)-(hdr1->threadid!=0));
}

/*
 * Change the fields of the supplied VDIF header to emulate the arrival of the
 * subequent frame/packet.
 */
void increment_vdif_header (vdif_header* hdr)
{
  int old_thread = hdr->threadid != 0; // our internal thread id
  if (old_thread == 1)
  { 
    // if we are on thread 1, we increment frame and switch to thread 0
    hdr->frame += 1;
    hdr->threadid = 0;
    
    // check if we wrapped to the next second
    if (FRAMESPERSEC==hdr->frame)
    {
      hdr->seconds += 1;
      hdr->frame = 0;
    }
  }
  else
    // otherwise, simply switch from thread 0 to thread 1
    hdr->threadid = 1;
}

int main(int argc, char** argv)
{
  // register SIGINT handling
  signal (SIGINT, sigint_handler);

  // register exit function
  atexit (exit_handler);

  int exit_status = EXIT_SUCCESS;
  key_t key = 0x40;
  mlog = multilog_open ("writer",0);
  hdu = dada_hdu_create (mlog);

  // get this many bytes at a time from data socket
  int FRAME_BUFSIZE = UDP_HDR_SIZE + VDIF_FRAME_SIZE; 
  char frame_buff[FRAME_BUFSIZE];
  char* vdif_buff = frame_buff + UDP_HDR_SIZE;
  vdif_header* vdhdr = (vdif_header*) (vdif_buff);
  char fill_vdif_buff[VDIF_FRAME_SIZE] = {127};
  vdif_header* fill_vdhdr = (vdif_header*) fill_vdif_buff;
  char dev[16] = ETHDEV;
  char hostname[MAXHOSTNAME];
  gethostname(hostname,MAXHOSTNAME);
  //char *starthost = strstr(hostname, "difx");

  // keep track of received frames for checking for packet loss, and also
  // of latest time based on VDIF packets (updates every 1s)
  int last_frame[2] = {-2,-2};
  int time_since_timet = 10;
  time_t packet_monitor_timet = 0;
  int obs_length = 0;
  double vdif_dmjd = -1;
  
  int state = STATE_STOPPED;
  //uint64_t port = WRITER_SERVICE_PORT;
  //uint64_t iport = WRITER_INFO_PORT;
  char cmd = CMD_NONE;
  int arg, maxsock = 0; 
  int stderr_output = 0;
  
  while ((arg = getopt(argc, argv, "hk:e:o")) != -1) {
    switch(arg) {

    case 'h':
      usage ();
      return 0;
      
    case 'k':
      if (sscanf (optarg, "%x", &key) != 1) {
        fprintf (stderr, "writer: could not parse key from %s\n", optarg);
        return -1;
      }
      break;

    case 'e':
      if (sscanf (optarg, "%s", dev) != 1) {
        fprintf (stderr, 
            "writer: could not parse ethernet device from %s\n", optarg);
        return -1;
      }
      break;

    case 'o':
      stderr_output = 1;
      break;
    }
  }

  fd_set readfds;
  struct timeval tv_500mus; //timeout for select()
  tv_500mus.tv_sec = 0;
  tv_500mus.tv_usec = 500; // set this to block long enough to receive a VLITE packet

  // empirically a reasonable timeout for raw socket when data are flowing
  struct timeval tv_10ms;
  tv_10ms.tv_sec = 0;
  tv_10ms.tv_usec = 10000;

  // for nanosleep
  ts_1ms = get_ms_ts(1);

  // observation metadata
  ScanInfoDocument D; //multicast message struct
  int od_cache_idx = 0;
  ObservationDocument* current_od = NULL;
  ObservationDocument od_cache[OD_CACHE_SIZE];

  // logging
  time_t currt;
  struct tm tmpt;
  currt = time (NULL);
  localtime_r (&currt,&tmpt);
  char currt_string[128];
  strftime (currt_string,sizeof(currt_string), "%Y%m%d_%H%M%S", &tmpt);
  currt_string[15] = 0;
  char logfile[128];
  pid_t pid = getpid ();
  snprintf (logfile,128,
      "%s/%s_%s_writer_%06d.log",LOGDIR,currt_string,hostname,pid);
  //FILE *logfile_fp = fopen (logfile, "w");
  logfile_fp = fopen (logfile, "w");
  multilog_add (mlog, logfile_fp);
  if (stderr_output)
    multilog_add (mlog, stderr);
  printf("writing log to %s\n",logfile);

  // print out command used to invoke
  char* incoming_hdr = &frame_buff[0]; // borrow this buffer
  strncpy (incoming_hdr, argv[0], FRAME_BUFSIZE-1);
  for (int i = 1; i < argc; ++i) {
    strcat (incoming_hdr, " ");
    strncat (incoming_hdr, argv[i], 4095-strlen(incoming_hdr));
  }
  multilog (mlog, LOG_INFO, "[WRITER] invoked with: \n%s\n", incoming_hdr);

  // have a 1-1 mapping between threadios and recent dump times.  Given
  // size of the buffer, about 100 seems appropriate.
  // Since the buffers are 1s, use UTC time stamps.
  threadio_t* threadios[MAX_THREADIOS] = {NULL};
  time_t dump_times[MAX_THREADIOS] = {0};
  int dump_idx = 0; // index of next dump / threadio
  // this really needs to be the size of the number of buffers!
  // but 100 here should be plenty
  char* bufs_to_write[100] = {NULL};
  time_t times_to_write[100] = {0};

  //Try to connect to existing data block
  dada_hdu_set_key (hdu,key);
  if (dada_hdu_connect(hdu) < 0)
  {
    multilog (mlog, LOG_ERR, "Unable to connect to psrdada HDU.\n");
    exit (EXIT_FAILURE);
  }
  multilog (mlog, LOG_INFO, "Connected to psrdada HDU.\n");

  
  // set up multi-cast sockets
  int mc_control_sock = open_mc_socket (mc_vlitegrp, MC_WRITER_PORT,
      "Control Socket [Writer]", &maxsock, mlog);
  char mc_from[24]; //ip address of multicast sender
  char mc_control_buf[32];
  
  int mc_trigger_sock = open_mc_socket (mc_vlitegrp, MC_TRIGGER_PORT,
      "Trigger Socket [Writer]", &maxsock, mlog);
  char mc_trigger_buf[sizeof(trigger_t)*MAX_TRIGGERS];
  trigger_t* trigger_queue[MAX_TRIGGERS] = {NULL};
  
  char mc_info_buf[MSGMAXSIZE];
  int mc_info_sock = open_mc_socket (mc_obsinfogrp, MULTI_OBSINFO_PORT,
      "Observation Info [Writer]", &maxsock, mlog);

  // Open raw data stream socket
  Connection raw;
  raw.alen = sizeof(raw.rem_addr);
  // size of data socket buffer; enough for 131 ms less overheads
  raw.sockoptval = 32*1024*1024;

  raw.svc = openRawSocket (dev,0);
  if (raw.svc < 0) {
    multilog (mlog, LOG_ERR, 
       "Cannot open raw socket on %s. Error code %d\n", dev, raw.svc);
    exit (EXIT_FAILURE);
  }
  setsockopt (raw.svc, SOL_SOCKET, SO_RCVBUF, (char *)&(raw.sockoptval), 
      sizeof(raw.sockoptval));
  if (raw.svc > maxsock)
    maxsock = raw.svc;
  // Set a short timeout on raw socket.  Native rate should be large.
  setsockopt (raw.svc, SOL_SOCKET, SO_RCVTIMEO,
      (const char *)&(tv_10ms), sizeof(tv_10ms));

  //int write_voltages = 0;
  uint64_t packets_written = 0;
  ssize_t raw_bytes_read = 0;
  ssize_t ipcio_bytes_written = 0;

  int break_loop = 0;

  while (1) {

    // new, simplified control loop handling stateless triggering.  (In
    // general, there is no reason to expect we might not have finished an
    // observation by the time we realize there is a transient in the data!)
    //
    // Thus, on each iteration, we check the control socket, the trigger
    // socket, and the raw data socket.  Regardless of state, we read any
    // available data.  However, if we are not in STATE_STARTED, we simply
    // discard it.  As before, when we read data we try to read 20 frames
    // simply to avoid calling select 20 times.  (We know there will be data.)

    FD_ZERO (&readfds);
    FD_SET (mc_control_sock, &readfds);
    FD_SET (mc_trigger_sock, &readfds);
    FD_SET (mc_info_sock, &readfds);
    FD_SET (raw.svc, &readfds);
    if (select (maxsock+1,&readfds,NULL,NULL,&tv_500mus) < 0)
    {
      multilog (mlog, LOG_ERR, "[Main Control Loop] Error calling select.");
      fflush (logfile_fp);
      exit (EXIT_FAILURE);
    }
    
    // if input is waiting on control socket, read it
    // in principle only commands we accept now are trigger and quit
    if (FD_ISSET (mc_control_sock, &readfds))
    {
      fprintf (stderr, "receiving command\n");
      fflush (stderr);
      cmd = check_for_cmd (mc_control_sock, mc_control_buf, 32, logfile_fp);
    }
    else
      cmd = CMD_NONE;

    // CMD_START can now only be issued by comparing packets to scan starts
    if ((cmd == CMD_START) && (state == STATE_STOPPED))
    {
      multilog (mlog, LOG_ERR, "[Main Control Loop] rxed CMD_START.");
      exit (EXIT_FAILURE);
    }

    // must accept CMD_STOP to handle skips on other nodes
    if (cmd == CMD_STOP)
      stop_observation (hdu, mlog, &packets_written, &state);
    
    //
    //if ((cmd == CMD_FAKE_START) && (state == STATE_STOPPED))
    if ((cmd == CMD_FAKE_START))
    {
      od_cache_idx += 1;
        if (od_cache_idx == OD_CACHE_SIZE)
          od_cache_idx = 0;
      fake_observation_document (&od_cache[od_cache_idx], vdif_dmjd, 0, 1);
      multilog (mlog, LOG_INFO,
          "Inserting a fake START observation document.\n");
    }

    //if ((cmd == CMD_FAKE_STOP) && (state == STATE_STARTED))
    if ((cmd == CMD_FAKE_STOP))
    {
      od_cache_idx += 1;
        if (od_cache_idx == OD_CACHE_SIZE)
          od_cache_idx = 0;
      fake_observation_document (&od_cache[od_cache_idx], vdif_dmjd, 1, 0);
      multilog (mlog, LOG_INFO,
          "Inserting a fake FINISH observation document.\n");
    }

    if (cmd == CMD_QUIT)
    {
      if (state == STATE_STARTED)
      {
        multilog (mlog, LOG_INFO,"[STATE_STARTED->QUIT], exiting.\n");
        multilog (mlog, LOG_INFO,"dada_hdu_unlock_write result = %d.\n",
            dada_hdu_unlock_write (hdu));
      }
      else if (state == STATE_STOPPED)
        multilog (mlog, LOG_INFO, "[STATE_STOPPED->QUIT], exiting.\n");
      break;
    }

    // check for a new scan message
    if (FD_ISSET (mc_info_sock, &readfds))
    {
      int nbytes = MultiCastReceive (
          mc_info_sock, mc_info_buf, MSGMAXSIZE, NULL);
      if (nbytes <= 0) {
        multilog (mlog, LOG_ERR,
            "Error on Obsinfo socket, return value = %d\n",nbytes);
        exit (EXIT_FAILURE);
      }
      parseScanInfoDocument (&D, mc_info_buf);
      if (D.type == SCANINFO_OBSERVATION)
      {
        multilog (mlog, LOG_INFO, "Incoming ObservationDocument:\n");
        //printScanInfoDocument (&D);
        fprint_observation_document(stderr, &D.data.observation, 1);
        od_cache_idx += 1;
        if (od_cache_idx == OD_CACHE_SIZE)
          od_cache_idx = 0;
        ObservationDocument* od = &od_cache[od_cache_idx];
        // check to see if we have run out of room in the OD cache
        if (od == current_od)
        {
          multilog (mlog, LOG_ERR,
              "Attempting to overwrite current ObservationDocument!");
          exit( EXIT_FAILURE);
        }
        memcpy (od, &D.data.observation, sizeof(ObservationDocument));
        // write to log or file?
      }
    }

    // check for triggers
    if (FD_ISSET (mc_trigger_sock, &readfds))
    {
      int nbytes = MultiCastReceive (mc_trigger_sock, mc_trigger_buf,
          sizeof(trigger_t)*MAX_TRIGGERS, mc_from);
      int ntrigger = nbytes / sizeof(trigger_t);
      multilog (mlog, LOG_INFO, "Received %d triggers.\n", ntrigger);
      if (ntrigger > MAX_TRIGGERS)
      {
        multilog (mlog, LOG_INFO, "Truncating trigger list to fit.\n");
        ntrigger = MAX_TRIGGERS;
      }
      if ((nbytes == 0)  || ((nbytes-ntrigger*sizeof(trigger_t))!=0))
      {
        multilog (mlog, LOG_ERR,
            "Error on Trigger socket, return value = %d\n",nbytes);  
        continue;
      }
      for (int i=0; i < ntrigger; ++i)
      {
        fprintf (stderr, "working on trigger %d\n",i);
        trigger_t* trig = (trigger_t*) (mc_trigger_buf) + i;
        trigger_queue[i] = trig;
        multilog (mlog, LOG_INFO,
            "(writer) received a trigger with t0=%.3f and t1=%.3f and meta %s.\n",
            trig->t0,trig->t1,trig->meta);
      }
    }

    // process triggers
    if (trigger_queue[0] != NULL)
    {
      fprintf (stderr, "here5\n");
      fflush (stderr);
      // determine which buffers to dump
      ipcbuf_t* buf = (ipcbuf_t *) hdu->data_block;; 
      int nbufs = ipcbuf_get_nbufs(buf);
      int bufsz = ipcbuf_get_bufsz(buf);
      int nbufs_to_write = get_buffer_trigger_overlap (
          buf->buffer, nbufs, trigger_queue,
          dump_times, bufs_to_write, times_to_write);
      multilog (mlog, LOG_INFO,
          "Dumping out %d buffers with the following UTC timestamps:\n",nbufs_to_write);
      for (int ibuf=0; ibuf < nbufs_to_write; ++ibuf)
      {
        threadio_t* tio = threadios[dump_idx];
        if (tio != NULL)
        {
          // check status of previous I/O and clean up memory
          if (tio->status != 0)
          {
            if (tio->status < 0)
              multilog (mlog, LOG_ERR,
                  "Previous I/O still ongoing, likely an error condition.");
            else
              multilog (mlog, LOG_ERR,
                  "Previous I/O encountered error.");
          }
          free (tio);
        }

        tio = (threadio_t*)(malloc (sizeof (threadio_t)));
        if (tio==NULL)
        {
          multilog (mlog, LOG_ERR, "Unable to allocate ThreadIO mem.");
          exit( EXIT_FAILURE);
        }
        tio->status = -1;
        tio->buf = bufs_to_write[ibuf];
        tio->bufsz = bufsz;
        tio->ms_delay = ibuf > 1?ibuf*500:0;
        int sid = getVDIFStationID ((vdif_header*)bufs_to_write[ibuf]);
        snprintf (tio->fname, 255,"%s/%s_ea%02d_%li.vdif",
            EVENTDIR,currt_string,sid,times_to_write[ibuf]);
        multilog (mlog, LOG_INFO, ".... UTC %li writing to %s\n",
            times_to_write[ibuf], tio->fname);
        // TODO -- could actually do a better job of this since we KNOW how far
        // a given buffer is from dropping off the end of the ring buffer, and
        // rather than launching the threads all at once, we can either launch them
        // staggered over time (on the 1-s boundaries) or at least in between 20-packet
        // reads so that we don't overflow the raw socket buffer
        // TODO also -- we can make the dump threads less problematic if we use
        // zero copy by memory mapping the output file to the system memory and
        // faulting it
        if (pthread_create (&tio->tid, NULL, &buffer_dump, (void*)tio) != 0)
            multilog (mlog, LOG_ERR, "pthread_create failed\n");
        if (pthread_detach (tio->tid) != 0)
            multilog (mlog, LOG_ERR, "pthread_detach failed\n");
        threadios[dump_idx] = tio;
        dump_times[dump_idx++] = times_to_write[ibuf];
        if (MAX_THREADIOS == dump_idx)
          dump_idx = 0;
      }
      // remove triggers from queue
      for (int i=0; i < MAX_TRIGGERS; ++i)
        trigger_queue[i] = NULL;
    }

    // read packets from the raw socket
    if (FD_ISSET (raw.svc, &readfds))
    {
      // To reduce CPU utilization, read multiple packets.
      for (int ipacket = 0; ipacket < 20; ++ipacket) {

        raw_bytes_read = recvfrom (raw.svc, frame_buff, FRAME_BUFSIZE, 0, 
            (struct sockaddr *)&(raw.rem_addr), &(raw.alen));

        if (raw_bytes_read != FRAME_BUFSIZE)
        {
          if (raw_bytes_read <= 0)
          {
            multilog (mlog, LOG_ERR,
                "Raw socket read failed: %d\n.", raw_bytes_read);
            //exit (EXIT_FAILURE);
            // change 9/5/2019 -- break out of multi-packet loop instead of failing
            break;

          }
          else
          {
            //multilog (mlog, LOG_ERR, "Received anomalous packet size: %d, aborting obs.\n", raw_bytes_read);
            multilog (mlog, LOG_ERR,
                "Received anomalous packet size: %d.  Continuing.\n", raw_bytes_read);
            //stop_observation (hdu, mlog, &packets_written, &state);
            //abnormal_stop = 1;
            continue;
          }
        }

        int read_thread = getVDIFThreadID (vdhdr) != 0;
        int read_frame = getVDIFFrameNumber (vdhdr);
        vdif_header* vdhdr_to_use;
        char* vdbuff_to_use;
        int frame_diff;

        // handle special case of first frame
        if (last_frame[read_thread] == -2)
          frame_diff = 1;
        else 
          frame_diff = vdif_frame_difference (fill_vdhdr, vdhdr);

        // TMP -- sanity check
        if (frame_diff < 1)
        {
          fprintf (stderr, "frame difference is %d, should not be possible!\n",frame_diff);
          fprintf (stderr, "hdr1: sec, frame, threadid=%d %d %d\n",fill_vdhdr->seconds,fill_vdhdr->frame,fill_vdhdr->threadid);
          fprintf (stderr, "hdr2: sec, frame, threadid=%d %d %d\n",vdhdr->seconds,vdhdr->frame,vdhdr->threadid);
        }

        if (frame_diff > 1)
          multilog (mlog, LOG_ERR, "Found %d skipped frames.\nCurrent  frame=%d, thread=%d.\nPrevious frame=%d, thread=%d.\n",frame_diff-1,read_frame,read_thread,fill_vdhdr->frame,fill_vdhdr->threadid!=0);
        
        for (int iframe = frame_diff; iframe > 0; --iframe)
        {
          if (iframe>1)
          {
            // increment the last recorded header and use it for all processing
            increment_vdif_header (fill_vdhdr);
            vdhdr_to_use = fill_vdhdr;
            vdbuff_to_use = fill_vdif_buff;
          }
          else
          {
            vdhdr_to_use = vdhdr;
            vdbuff_to_use = vdif_buff;
          }

          int thread = getVDIFThreadID (vdhdr_to_use) != 0;
          int current_frame = getVDIFFrameNumber (vdhdr_to_use);


          // check if we are on a one-second boundary with thread==0
          if ((thread==0) && (MAXFRAMENUM == last_frame[0]))
          {
            vdif_dmjd = getVDIFFrameDMJD (vdhdr_to_use, FRAMESPERSEC);
            
            // reset sample counter
            last_frame[0] = -1;
            last_frame[1] = -1;

            // check synchronization between system time and packets every 10s
            if (time_since_timet==10)
              packet_monitor_timet = time (NULL);
            time_since_timet -= 1;
            if (time_since_timet==0)
            {
              time_t old_time = packet_monitor_timet;
              packet_monitor_timet = time (NULL);
              if ((packet_monitor_timet-old_time) > 11)
              {
                fprintf (stderr, "packet_monitor_time = %ld; old_time = %ld\n",packet_monitor_timet,old_time);
                multilog (mlog, LOG_ERR, "Packet times and system time are out of synch by more than 1s!\n");
                exit (EXIT_FAILURE);
              }
              time_since_timet = 10;
            }

            // monitor output buffer -- if full, abort
            ipcio_t* ipc = hdu->data_block;
            uint64_t m_nbufs = ipcbuf_get_nbufs ((ipcbuf_t *)ipc);
            uint64_t m_full_bufs = ipcbuf_get_nfull((ipcbuf_t*) ipc);
            if (m_full_bufs == (m_nbufs - 2))
            {
              multilog (mlog, LOG_ERR, "[WRITER] Only one free buffer left!  Aborting process.");
              exit (EXIT_FAILURE);
            }

            if (state == STATE_STARTED)
              obs_length += 1;

            // See if there is a new scan, defined such that the VDIF second
            // and the startTIME ceilinged to the next second match
            ObservationDocument* new_scan = search_od_cache (
                &od_cache[0], OD_CACHE_SIZE, vdif_dmjd);

            if (new_scan != NULL)
            {
              multilog (mlog, LOG_INFO, "Found a new scan at %10.8f.\n",
                  new_scan->startTime);

              // If recording, check to see if new scan is consistent.
              // with the old one.  If not, stop the observation.
              int new_obs = 1;
              if (state==STATE_STARTED)
              {
                if (!check_od_consistency (current_od, new_scan, mlog, obs_length))
                {
                  // stop current observation if they are inconsistent
                  stop_observation (hdu, mlog, &packets_written, &state);
                  // Check if new scan is FINISH
                  if (strcasecmp (new_scan->name, "FINISH") == 0)
                    new_obs = 0;
                }
                else
                  new_obs = 0;
              }

              current_od = new_scan;

              // Start new observation
              if (new_obs)
              {
                multilog (mlog, LOG_INFO,"[STATE_STOPPED->START]\n");
                state = STATE_STARTED;
                packets_written = 0;
                obs_length = 0;

                if (dada_hdu_lock_write (hdu) < 0) {
                  multilog (mlog, LOG_ERR, 
                      "Unable to lock psrdada HDU, exiting\n");
                  exit (EXIT_FAILURE);
                }
              
                multilog (mlog, LOG_INFO, "Writing psrdada header.\n");
                write_psrdada_header (hdu, vdhdr_to_use, current_od);
              }
            }

          }
          
          if (iframe==1)
            // copy over VDIF header for next loop
            *fill_vdhdr = *vdhdr;

          // update sample counter
          last_frame[thread] = current_frame;

          if (state != STATE_STARTED)
            continue;

          // TMP sanity check
          //if (vdbuff_to_use == fill_vdif_buff)
          //  fprintf (stderr, "Writing out of fill buffer.\n");

          ipcio_bytes_written = ipcio_write (
              hdu->data_block, vdbuff_to_use, VDIF_FRAME_SIZE);

          if (ipcio_bytes_written != VDIF_FRAME_SIZE) {
            multilog (mlog, LOG_ERR,
                "Failed to write VDIF packet to psrdada buffer.\n");
            exit (EXIT_FAILURE);
          }
          else
            packets_written++;

        } // end loop over received (and possibly) skipped) frame(s)
      } // end multiple packet read
    } // end raw socket read logic

    if (break_loop)
    {
      multilog (mlog, LOG_INFO, "Breaking main loop by request.\n");
      stop_observation (hdu, mlog, &packets_written, &state);
    }
  } // end main loop over state/packets
  
  for (int iio=0; iio < MAX_THREADIOS; ++iio)
    if (threadios[iio] != NULL)
      free (threadios[iio]);

  return exit_status;
} // end main

