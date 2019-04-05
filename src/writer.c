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

#include "utils.h"
#include "executor.h"
#include "def.h"
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

static threadio_t** threadios = NULL;
static int nthreadios = 0;

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
      
   /*
    case 'p':
      if (sscanf (optarg, "%"PRIu64"", &port) != 1) {
        fprintf (stderr, "writer: could not parse port from %s\n", optarg);
        return -1;
      }
      break;

    case 'i':
      if (sscanf (optarg, "%"PRIu64"", &iport) != 1) {
        fprintf (stderr, "writer: could not parse port from %s\n", optarg);
        return -1;
      }
      break;
    */
    
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

  struct timeval tv_1s; //timeout for raw socket
  tv_500mus.tv_sec = 1;
  tv_500mus.tv_usec = 0;

  Connection raw;
  raw.alen = sizeof(raw.rem_addr);
  raw.sockoptval = 16*1024*1024; //size of data socket internal buffer

  // observation metadata
  ObservationDocument od;

  // logging
  time_t currt;
  struct tm tmpt;
  currt = time(NULL);
  localtime_r(&currt,&tmpt);
  char currt_string[128];
  strftime(currt_string,sizeof(currt_string), "%Y%m%d_%H%M%S", &tmpt);
  currt_string[15] = 0;
  char logfile[128];
  pid_t pid = getpid();
  snprintf (logfile,128,
      "%s/%s_%s_writer_%06d.log",LOGDIR,currt_string,hostname,pid);
  //FILE *logfile_fp = fopen (logfile, "w");
  logfile_fp = fopen (logfile, "w");
  multilog_add (log, logfile_fp);
  if (stderr_output)
    multilog_add (log, stderr);
  printf("writing log to %s\n",logfile);

  // keep track of last N (~200) dump times.  The idea is both to avoid
  // dumping out the same data twice, and to potentially place a limit on
  // the number of dumps per unit time.  In current implementation, each
  // buffer is 1s, so keeping track of UTC time stamps is perfect.
  size_t size_dump_times=300, num_dump_times=0 ,idx_dump_times=0;
  time_t dump_times[size_dump_times];
  memset (dump_times, 0, sizeof(time_t)*size_dump_times);

  // set up some limits on dumps; a reasonable guess might be no more than
  // 30s of dumped data every five minutes
  size_t dump_window = 120;
  size_t max_dump_window = 60;

  dada_hdu_set_key (hdu,key);

  //Try to connect to existing data block
  if (dada_hdu_connect(hdu) < 0)
  {
    multilog (log, LOG_ERR, "Unable to connect to psrdada HDU.\n");
    exit (EXIT_FAILURE);
  }
  multilog (log, LOG_INFO, "Connected to psrdada HDU.\n");

  // set up multi-cast control socket
  int mc_control_sock = 0;
  char mc_control_buf[32];
  mc_control_sock = openMultiCastSocket ("224.3.29.71", 20001);
  if (mc_control_sock < 0) {
    multilog (log, LOG_ERR, "Failed to open Observation multicast; openMultiCastSocket = %d\n",mc_control_sock);
    exit (EXIT_FAILURE);
  }
  else
    multilog (log, LOG_INFO, "Control socket: %d\n",mc_control_sock);
  if (mc_control_sock > maxsock)
    maxsock = mc_control_sock;

  //Open raw data stream socket
  raw.svc = openRawSocket (dev,0);
  if(raw.svc < 0) {
    multilog (log, LOG_ERR, 
       "Cannot open raw socket on %s. Error code %d\n", dev, raw.svc);
    exit (EXIT_FAILURE);
  }
  setsockopt (raw.svc, SOL_SOCKET, SO_RCVBUF, (char *)&(raw.sockoptval), 
      sizeof(raw.sockoptval));
  if (raw.svc > maxsock)
    maxsock = raw.svc;
  // TODO -- see if we need this -- set timeout on raw socket read
  // (I doubt we need this unless things are really broken.  Could use it
  // to check for a network failure and shut down, I suppose.)
  setsockopt (raw.svc, SOL_SOCKET, SO_RCVTIMEO, (const char *)&(tv_1s), sizeof(tv_1s));

  // set up multi-cast observation information socket
  int mc_info_sock = 0;
  char mc_info_buf[1024];
  mc_info_sock = openMultiCastSocket ("224.3.29.71", 20002);
  if (mc_info_sock < 0) {
    multilog (log, LOG_ERR, "Failed to open Observation multicast; openMultiCastSocket = %d\n",mc_info_sock);
    exit (EXIT_FAILURE);
  }
  else
    multilog (log, LOG_INFO, "Obsinfo socket: %d\n",mc_info_sock);

  // skip_frames is a klugey construct to handle the continual start/stop
  // of reading from the raw socket.  Between observations, the buffer will
  // fill and packets will drop, so the first packets will be old data.
  // Current solution is to simply skip enough frames that the buffer is
  // exhausted.  A better way might be to set the buffer to 0 size when
  // we aren't clocking samples, or to continually clock them in and use
  // valid flags for start/stop of observation.  However, tests with this
  // suggest the OS doesn't actually change the buffer size, so we're stuck
  // with the kluge.
  int skip_frames = 0;
  int write_header = 0;
  int write_voltages = 0;
  uint64_t packets_written = 0;
  uint64_t voltage_packets_written = 0;
  uint64_t voltage_packet_counter = 0;
  ssize_t info_bytes_read = 0;
  ssize_t raw_bytes_read = 0;
  ssize_t ipcio_bytes_written = 0;
  time_t dump_time = 0;

  while (1) {

    if (state == STATE_STOPPED) {
      if (cmd == CMD_NONE)
        cmd = check_for_cmd (mc_control_sock, mc_control_buf, 32, logfile_fp);

      if (cmd == CMD_START) 
      {
        state = STATE_STARTED;
        skip_frames = 0;
        packets_written = 0;
        write_voltages = 0;
        voltage_packets_written = 0;
        voltage_packet_counter = 0;
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
        cmd = CMD_NONE;
        //
        // check for a set source to write voltages on
        write_voltages = voltage_check_name(od.name,od.datasetId);
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

      }
      else if (cmd == CMD_EVENT) {
        multilog (log, LOG_INFO, "[STATE_STOPPED] ignored CMD_EVENT.\n");
        cmd = CMD_NONE;
      }
      else if (cmd == CMD_STOP) {
        multilog (log, LOG_INFO, "[STATE_STOPPED] ignored CMD_STOP.\n");
        cmd = CMD_NONE;
      }
      else if (cmd == CMD_QUIT) {
        multilog (log, LOG_INFO, "[STATE_STOPPED] received CMD_QUIT, exiting.\n");
        break;
      }
    } // end STATE_STOPPED logic

    //if state is STARTED, poll the command listening socket and the VDIF raw data socket
    if (state == STATE_STARTED) {

      // this construct adds the two sockets to the FDs and then blocks
      // until input is available on one or the other; multiplexing
      FD_ZERO (&readfds);
      FD_SET (mc_control_sock, &readfds);
      FD_SET (raw.svc, &readfds);
      if (select (maxsock+1,&readfds,NULL,NULL,&tv_500mus) < 0)
      {
        multilog (log, LOG_ERR, "[STATE_STARTED] Error calling select.");
        fflush (logfile_fp);
        exit (EXIT_FAILURE);
      }
      
      //if input is waiting on listening socket, read it
      if (FD_ISSET (mc_control_sock,&readfds)) {
        cmd = check_for_cmd (mc_control_sock, mc_control_buf, 32, logfile_fp);
      }

      // note to future self -- dump_time is not actually used below, I think
      // it's mostly just a shorthand for a triggered dump, so refactor to
      // make it more apparent!

      // check if we are enough packets removed from an automatic dump
      // request and, if so, execute it; wait 10s
      int do_dump = (dump_time && (packets_written > 25600*20));
      if (do_dump)
        dump_time = 0;

      // only trigger a "write voltage" dump every 10 seconds for overhead
      if (write_voltages && (packets_written > 25600*16) && ((packets_written-voltage_packet_counter)>(25600*4)))
      {
        do_dump = 1;
        voltage_packet_counter = packets_written;
        multilog (log, LOG_INFO, "Issuing a new 10-s dump for %s.\n",od.name);
      }


      if (cmd == CMD_EVENT)
      {
        do_dump += 1; // will be 1 for CMD_EVENT or 2 for triggered dump
        cmd = CMD_NONE;
        // CMD_EVENT overrides dump_time
        dump_time = 0;
      }

      if (do_dump)
      //if (0) // TMP -- disable dumps
      {
        
        // dump voltages to file
        
        // check to see if there is a previous dump ongoing
        int ready_to_write = 1;

        for (int iio=0; iio < nthreadios; ++iio)
        {
          if (threadios[iio]->status != 0)
          {
            // previous dump still happening
            multilog (log, LOG_INFO,
                "Previous I/O still ongoing, skipping this CMD_EVENT.");
            ready_to_write = 0;
            break;
          }
          if (threadios[iio]->status > 0)
          {
            // some problem with previous dump; skip new dumps for now
            multilog (log, LOG_ERR,
                "Previous I/O errored, skipping this CMD_EVENT.");
            ready_to_write = 0;
            break;
          }
        }

        if (!ready_to_write)
          continue;

        multilog (log, LOG_INFO, "Launching threads for buffer dump.\n");

        // free previous structs
        for (int iio=0; iio < nthreadios; ++iio)
          if (threadios[iio])
            free (threadios[iio]);
        if (threadios)
          free (threadios); threadios = NULL;
        nthreadios = 0;

        // determine which buffers to write
        ipcbuf_t* buf = (ipcbuf_t *) hdu->data_block;; 
        char** buffer = buf->buffer;
        int nbufs = ipcbuf_get_nbufs(buf);
        int bufsz = ipcbuf_get_bufsz(buf);

        // for testing, pretend trigger time is current time -14s
        // with a 20s transient
        time_t trigger_time = time (NULL)-1;
        //double trigger_len = 20;
        time_t tmin = trigger_time - 1;
        //time_t tmax = trigger_time + trigger_len + 8;
        time_t tmax = trigger_time;
        printf ("tmin = %ld tmax = %ld\n",tmin,tmax);
        char* bufs_to_write[32] = {NULL};

        // check to see whether we have exceeded the dump limit; if we have
        // a triggered dump, ignore the limit
        if (do_dump==1)
        {
          int dumps_in_window = 0;
          for (int idump=0; idump < num_dump_times; ++idump)
          {
            dumps_in_window += abs (trigger_time-dump_times[idump]) < dump_window;
          } 
          multilog (log, LOG_INFO, "Found %d seconds of dumps in the window.\n",
              dumps_in_window);
          if (dumps_in_window > max_dump_window)
          {
            multilog (log, LOG_INFO, "Too many dumps, skipping this trigger.\n");
            continue;
          }
        }

        for (int ibuf=0; ibuf < nbufs; ++ibuf)
        {
          // get time for frame at start of buffer
          vdif_header* vdhdr = (vdif_header*) buffer[ibuf];
          time_t utc_seconds = vdif_to_unixepoch (vdhdr);

          // TMP
          //
          /*
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
          //
          // end TMP

          // first condition -- within trigger time (TODO -- if we want
          // a voltage recording mode, could override this)
          int cond1 = (utc_seconds >= tmin) && (utc_seconds <= tmax);

          // second condition -- see if we have already dumped it
          int cond2 = 1;
          for (int idump=0; idump < num_dump_times; ++idump)
          {
            if (utc_seconds == dump_times[idump])
            {
              cond2 = 0;
              break;
            }
          }
          if (cond1 && cond2)
          {
            // dump this buffer
            bufs_to_write[nthreadios++] = buffer[ibuf];
            multilog (log, LOG_INFO, "dumping buffer %02d\n", ibuf);
            num_dump_times += num_dump_times < size_dump_times;
            if (idx_dump_times == size_dump_times)
              idx_dump_times = 0;
            dump_times[idx_dump_times++] = utc_seconds;
            voltage_packets_written += 2*25600;
          }  
        } // end loop over buffers

        // TODO -- replace this file name with appropriate format
        currt = time (NULL);
        gmtime_r (&currt,&tmpt);
        strftime (
            currt_string,sizeof(currt_string), "%Y%m%d_%H%M%S", &tmpt);
        *(currt_string+15) = 0;
        multilog (log, LOG_INFO, "GM time is %s\n", currt_string);

        // make new threadio_ts and pthreads
        threadios = malloc (sizeof (threadio_t*) * nthreadios);
        for (int iio=0; iio < nthreadios; ++iio)
        {
          threadios[iio] = malloc (sizeof (threadio_t));
          threadios[iio]->status = -1;
          threadios[iio]->buf = bufs_to_write[iio];
          threadios[iio]->bufsz = bufsz;
          int sid = getVDIFStationID ((vdif_header*)bufs_to_write[iio]);
          snprintf (threadios[iio]->fname, 255,
              "%s/%s_ea%02d_buff%02d.vdif",
              EVENTDIR,currt_string,sid,iio);
          if (pthread_create (&threadios[iio]->tid, NULL, &buffer_dump, 
              (void*)threadios[iio]) != 0)
              multilog (log, LOG_ERR, "pthread_create failed\n");
          if (pthread_detach (threadios[iio]->tid) != 0)
              multilog (log, LOG_ERR, "pthread_detach failed\n");
        }

        // write out the observation document for this dump
        char dump_od_fname[256];
        snprintf (dump_od_fname, 255,
            "/home/vlite-master/mtk/events/%s.od", currt_string);
        FILE *dump_od_fp = fopen (dump_od_fname, "w");
        fprint_observation_document(dump_od_fp, &od);
        fclose (dump_od_fp);

        // TODO -- want to transfer antenna properties and write them out

        // and we're done
        multilog (log, LOG_INFO,
            "Launched %d threads for buffer dump.\n", nthreadios);

      } // end do_dump logic
      
      //CMD_STOP --> change state to STOPPED, close data block
      else if (cmd == CMD_STOP) {
        state = STATE_STOPPED;
        skip_frames = 0;
        if (dada_hdu_unlock_write (hdu) < 0) {
          multilog (log, LOG_ERR,
              "[STATE_STARTED->STOP]: unable to unlock psrdada HDU, exiting\n");
          exit_status = EXIT_FAILURE;
          break;
        }
        multilog (log, LOG_INFO, 
            "[STATE_STARTED->STOP] Wrote %d packets to psrdada buffer.\n",packets_written);
        packets_written = 0;
        write_voltages = 0;
        cmd = CMD_NONE;
        fflush (logfile_fp);
        continue;
      }

      //CMD_QUIT --> close data block, shutdown listening socket, return
      else if (cmd == CMD_QUIT) {
        multilog (log, LOG_INFO,
            "[STATE_STARTED->QUIT], exiting.\n");
        multilog (log, LOG_INFO,
            "dada_hdu_unlock_write result = %d.\n",
            dada_hdu_unlock_write (hdu));
        break;
      }

      else if (cmd == CMD_START) {
        multilog (log, LOG_INFO,
            "[STATE_STARTED->START] ignored CMD_START.\n");
        cmd = CMD_NONE;
      }

      //read a packet from the data socket and write it to the ring buffer
      if (FD_ISSET (raw.svc,&readfds)) {
        
        // This is a little kluge to reduce CPU utilization.
        // If we are in START_STARTED, read multiple packets before 
        // looping back to multiplex and check for commands
        for (int ipacket = 0; ipacket < 20; ++ipacket) {

          raw_bytes_read = recvfrom (raw.svc, buf, BUFSIZE, 0, 
              (struct sockaddr *)&(raw.rem_addr), &(raw.alen));

          if (raw_bytes_read == BUFSIZE) {

            // this makes sure we read 0.5s of data from the buffer; based
            // on the maximum requested raw socket size (2x16MB) this is
            // more than enough to clear out any old data from the buffer
            if (skip_frames <  25600) {
              skip_frames ++;
              continue;
            }

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
    } // end STATE_STARTED logic

  } // end main loop over state/packets
  
  shutdown (mc_control_sock,2);
  shutdown (mc_info_sock,2);

  multilog (log, LOG_INFO,
      "dada_hdu_disconnect result = %d.\n", dada_hdu_disconnect (hdu));
  fclose (logfile_fp);
  //change_file_owner (logfile_fp, "vliteops");

  for (int iio=0; iio < nthreadios; ++iio)
    free (threadios[iio]);
  free (threadios);

  return exit_status;
} // end main

