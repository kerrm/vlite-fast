#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <pthread.h>

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

static FILE* logfile_fp = NULL;

void usage ()
{
  fprintf(stdout,"Usage: writer [options]\n"
	  "-k hexadecimal shared memory key  (default: 40)\n"
	  "-e ethernet device id (default: eth0)\n"
	  "-o print logging messages to stderr (as well as logfile)\n"
    );
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

void fprint_observation_document(FILE* fd, ObservationDocument* od)
{
  fprintf(fd,"    datasetId = %s\n", od->datasetId);
  fprintf(fd,"    configId = %s\n", od->configId);
  fprintf(fd,"    startTime = %10.8f\n", od->startTime);
  fprintf(fd,"    name = %s\n", od->name);
  fprintf(fd,"    ra = %10.8f\n", od->ra);
  fprintf(fd,"    dec = %10.8f\n", od->dec);
  fprintf(fd,"    dra = %10.8f\n", od->dra);
  fprintf(fd,"    ddec = %10.8f\n", od->ddec);
  fprintf(fd,"    azoffs = %10.8f\n", od->azoffs);
  fprintf(fd,"    eloffs = %10.8f\n", od->eloffs);
  fprintf(fd,"    startLST = %10.8f\n", od->startLST);
  fprintf(fd,"    scanNo = %d\n", od->scanNo);
  fprintf(fd,"    subscanNo = %d\n", od->subscanNo);
  fprintf(fd,"    primaryBand = %s\n", od->primaryBand);
  fprintf(fd,"    usesPband = %d\n", od->usesPband);
}

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
    multilog (log, LOG_INFO, "buffer UTC %s\n", dada_utc);
    printf ("utc_seconds = %ld\n",utc_seconds);
    fprintVDIFHeader (stderr, vdhdr, VDIFHeaderPrintLevelShort);
    */

    int trigger_overlap = 0;
    for (int itrig=0; itrig < MAX_TRIGGERS; ++itrig)
    {
      if (tqueue[itrig] == NULL)
        break;
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

int main(int argc, char** argv)
{
  // register SIGINT handling
  signal (SIGINT, sigint_handler);

  int exit_status = EXIT_SUCCESS;
  key_t key = 0x40;
  multilog_t* log = multilog_open ("writer",0);
  dada_hdu_t* hdu = dada_hdu_create (log);

  //get this many bytes at a time from data socket
  int BUFSIZE = UDP_HDR_SIZE + VDIF_PKT_SIZE; 
  char buf[BUFSIZE];
  char dev[16] = ETHDEV;
  char hostname[MAXHOSTNAME];
  gethostname(hostname,MAXHOSTNAME);
  //char *starthost = strstr(hostname, "difx");
  
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

  // observation metadata
  ObservationDocument od;

  // logging
  time_t currt;
  struct tm tmpt;
  currt = time(NULL);
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
  multilog_add (log, logfile_fp);
  if (stderr_output)
    multilog_add (log, stderr);
  printf("writing log to %s\n",logfile);

  // have a 1-1 mapping between threadios and recent dump times.  Given
  // size of the buffer, about 100 seems appropriate.
  // Since the buffers are 1s, use UTC time stamps.
  threadio_t* threadios[MAX_THREADIOS] = {NULL};
  time_t dump_times[MAX_THREADIOS] = {0};
  int dump_idx = 0; // index of next dump / threadio
  char* bufs_to_write[MAX_TRIGGERS] = {NULL};
  time_t times_to_write[MAX_TRIGGERS] = {0};

  //Try to connect to existing data block
  dada_hdu_set_key (hdu,key);
  if (dada_hdu_connect(hdu) < 0)
  {
    multilog (log, LOG_ERR, "Unable to connect to psrdada HDU.\n");
    exit (EXIT_FAILURE);
  }
  multilog (log, LOG_INFO, "Connected to psrdada HDU.\n");

  
  // set up multi-cast sockets
  int mc_control_sock = open_mc_socket (mc_vlitegrp, MC_WRITER_PORT,
      "Control Socket [Writer]", &maxsock, log);
  char mc_from[24]; //ip address of multicast sender
  char mc_control_buf[32];
  
  int mc_trigger_sock = open_mc_socket (mc_vlitegrp, MC_TRIGGER_PORT,
      "Trigger Socket [Writer]", &maxsock, log);
  char mc_trigger_buf[sizeof(trigger_t)*MAX_TRIGGERS];
  trigger_t* trigger_queue[MAX_TRIGGERS] = {NULL};
  
  char mc_info_buf[10000]; // TODO check this size
  int mc_info_sock = open_mc_socket (mc_vlitegrp, MC_INFO_PORT,
      "Observation Info [Writer]", NULL, log);

  // Open raw data stream socket
  Connection raw;
  raw.alen = sizeof(raw.rem_addr);
  raw.sockoptval = 16*1024*1024; //size of data socket internal buffer

  raw.svc = openRawSocket (dev,0);
  if (raw.svc < 0) {
    multilog (log, LOG_ERR, 
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

  int write_header = 0;
  //int write_voltages = 0;
  uint64_t packets_written = 0;
  //uint64_t voltage_packets_written = 0;
  //uint64_t voltage_packet_counter = 0;
  ssize_t info_bytes_read = 0;
  ssize_t raw_bytes_read = 0;
  ssize_t ipcio_bytes_written = 0;

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
    FD_SET (raw.svc, &readfds);
    if (select (maxsock+1,&readfds,NULL,NULL,&tv_500mus) < 0)
    {
      multilog (log, LOG_ERR, "[Main Control Loop] Error calling select.");
      fflush (logfile_fp);
      exit (EXIT_FAILURE);
    }
    
    // if input is waiting on control socket, read it
    if (FD_ISSET (mc_control_sock, &readfds))
    {
      fprintf (stderr, "receiving command\n");
      fflush (stderr);
      cmd = check_for_cmd (mc_control_sock, mc_control_buf, 32, logfile_fp);
    }
    else
      cmd = CMD_NONE;

    // handle commands; NB that commands can be set in this loop as
    // well as come in on the control socket
    if ((cmd == CMD_START) && (state == STATE_STOPPED))
    {
      fprintf (stderr, "here1\n");
      fflush (stderr);
      multilog (log, LOG_INFO,"[STATE_STOPPED->START]\n");
      // start an observation
      state = STATE_STARTED;
      packets_written = 0;
      // TODO -- block here until we get the observation info
      // TODO -- need to handle case of multiple observation documents
      info_bytes_read = read (mc_info_sock, mc_info_buf, 1024);
      multilog (log, LOG_INFO, "Read %d bytes from fd %d, expected %d.\n", 
          info_bytes_read,mc_info_sock,sizeof(ObservationDocument));
      // TODO -- what could happen here if heimdall runs behind, or if 
      // we have a long transient dump, is that we receive multiple ODs?
      // so check for a multiple of sizeof(OD) and then use the last one
      if (info_bytes_read != sizeof(ObservationDocument)) {
        // TODO will have to decide how to handle this error
        multilog (log, LOG_ERR, "Not a valid ObservationDocument!\n");
      }
      else {
        memcpy (&od, mc_info_buf, sizeof(ObservationDocument));
        char result[2048];
        print_observation_document (result,&od);
        multilog (log, LOG_INFO, "ObservationDocument:\n%s", result);
      }
      if (dada_hdu_lock_write (hdu) < 0) {
        multilog (log, LOG_ERR, "Unable to lock psrdada HDU, exiting\n");
        exit (EXIT_FAILURE);
      }
      write_header = 1;
      multilog (log, LOG_INFO, "After dada_hdu_lock_write\n");
      fflush (logfile_fp);

      // TODO -- recreate dumping here, but using trigger socket instead
      /*
      if (write_voltages)
      {
        dump_time = 0;
        multilog (log, LOG_INFO, "Recording voltages for %s.\n",od.name);
      }

      // otherwise, check for a set source to trigger on
      int do_dump = (write_voltages==0) & dump_check_name (od.name,od.datasetId);
      // dump a maximum of 3 times on known sources
      if (do_dump > 0 && do_dump < 3) {
        dump_time = time (NULL);
        multilog (log, LOG_INFO, "Setting a dump time %ld for %s.\n",
            dump_time, od.name);
      }
      */
    }

    if ((cmd == CMD_STOP) && (state == STATE_STARTED))
    {
      fprintf (stderr, "here2\n");
      fflush (stderr);
      //CMD_STOP --> change state to STOPPED, close data block
      state = STATE_STOPPED;
      if (dada_hdu_unlock_write (hdu) < 0) {
        multilog (log, LOG_ERR,
            "[STATE_STARTED->STOP]: unable to unlock psrdada HDU, exiting\n");
        exit_status = EXIT_FAILURE;
        break;
      }
      multilog (log, LOG_INFO, 
          "[STATE_STARTED->STOP] Wrote %d packets to psrdada buffer.\n",packets_written);
      packets_written = 0;
      fflush (logfile_fp);
      continue;
    }

    if (cmd == CMD_QUIT)
    {
      fprintf (stderr, "here3\n");
      fflush (stderr);
      if (state == STATE_STARTED)
      {
        multilog (log, LOG_INFO,"[STATE_STARTED->QUIT], exiting.\n");
        multilog (log, LOG_INFO,"dada_hdu_unlock_write result = %d.\n",
            dada_hdu_unlock_write (hdu));
      }
      else if (state == STATE_STOPPED)
        multilog (log, LOG_INFO, "[STATE_STOPPED->QUIT], exiting.\n");
      break;
    }

    // check for triggers
    if (FD_ISSET (mc_trigger_sock, &readfds))
    {
      fprintf (stderr, "here4\n");
      fflush (stderr);
      int nbytes = MultiCastReceive (
          mc_trigger_sock, mc_trigger_buf, sizeof(trigger_t)*MAX_TRIGGERS, mc_from);
      int ntrigger = nbytes / sizeof(trigger_t);
      multilog (log, LOG_INFO, "Received %d triggers.\n", ntrigger);
      fprintf (stderr, "Received %d triggers.\n", ntrigger);
      fflush (stderr);
      if (ntrigger > MAX_TRIGGERS)
      {
        multilog (log, LOG_INFO, "Truncating trigger list to fit.\n");
        ntrigger = MAX_TRIGGERS;
      }
      if ((nbytes == 0)  || ((nbytes-ntrigger*sizeof(trigger_t))!=0))
      {
        multilog (log, LOG_ERR,
            "Error on Trigger socket, return value = %d\n",nbytes);  
        continue;
      }
      for (int i=0; i < ntrigger; ++i)
      {
        fprintf (stderr, "working on trigger %d\n",i);
        trigger_t* trig = (trigger_t*) (mc_trigger_buf) + i;
        trigger_queue[i] = trig;
        multilog (log, LOG_INFO,
            "(writer) received a trigger with t0=%.3f and t1=%3.f and meta %s.\n",
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
      multilog (log, LOG_INFO,
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
              multilog (log, LOG_ERR,
                  "Previous I/O still ongoing, likely an error condition.");
            else
              multilog (log, LOG_ERR,
                  "Previous I/O encountered error.");
          }
          free (tio);
        }

        tio = (threadio_t*)(malloc (sizeof (threadio_t)));
        tio->status = -1;
        tio->buf = bufs_to_write[ibuf];
        tio->bufsz = bufsz;
        tio->ms_delay = ibuf > 1?ibuf*500:0;
        int sid = getVDIFStationID ((vdif_header*)bufs_to_write[ibuf]);
        snprintf (tio->fname, 255,"%s/%s_ea%02d_%li.vdif",
            EVENTDIR,currt_string,sid,times_to_write[ibuf]);
        multilog (log, LOG_INFO, ".... UTC %li writing to %s\n",
            times_to_write[ibuf], tio->fname);
        if (pthread_create (&tio->tid, NULL, &buffer_dump, (void*)tio) != 0)
            multilog (log, LOG_ERR, "pthread_create failed\n");
        if (pthread_detach (tio->tid) != 0)
            multilog (log, LOG_ERR, "pthread_detach failed\n");
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

        raw_bytes_read = recvfrom (raw.svc, buf, BUFSIZE, 0, 
            (struct sockaddr *)&(raw.rem_addr), &(raw.alen));

        if (raw_bytes_read == BUFSIZE)
        {
          if (state != STATE_STARTED)
            continue; // move to next packet

          if (write_header) {
            multilog (log, LOG_INFO, "Writing psrdada header.\n");
            write_psrdada_header (
                hdu, (vdif_header *)(buf + UDP_HDR_SIZE), &od);
            write_header = 0;
          }

          ipcio_bytes_written = ipcio_write (
              hdu->data_block,buf+UDP_HDR_SIZE, VDIF_PKT_SIZE);

          if (ipcio_bytes_written != VDIF_PKT_SIZE) {
            multilog (log, LOG_ERR,
                "Failed to write VDIF packet to psrdada buffer.\n");
            exit_status = EXIT_FAILURE;
            break;
          }
          else
            packets_written++;
        }
        else if (raw_bytes_read <= 0) {
          multilog (log, LOG_ERR,
              "Raw socket read failed: %d\n.", raw_bytes_read);
          fflush (logfile_fp);
          cmd = CMD_STOP;
          break; // break from multi-packet loop
        }
        else {
          multilog (log, LOG_ERR,
              "Received packet size: %d, ignoring.\n", raw_bytes_read);
          fflush (logfile_fp);
          cmd = CMD_STOP;
          break; // break from multi-packet loop
        }
      } // end multiple packet read
    } // end raw socket read logic

  } // end main loop over state/packets
  
  shutdown (mc_control_sock,2);
  shutdown (mc_info_sock,2);

  multilog (log, LOG_INFO,
      "dada_hdu_disconnect result = %d.\n", dada_hdu_disconnect (hdu));
  fclose (logfile_fp);
  //change_file_owner (logfile_fp, "vliteops");

  for (int iio=0; iio < MAX_THREADIOS; ++iio)
    if (threadios[iio] != NULL)
      free (threadios[iio]);

  return exit_status;
} // end main

