#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <time.h>

#include "dada_def.h"
#include "dada_hdu.h"
#include "ipcio.h"
#include "ascii_header.h"
#include "multilog.h"

#include "utils.h"
#include "executor.h"
#include "def.h"
#include "vdifio.h"

static FILE* logfile_fp = NULL;

void usage ()
{
  fprintf(stdout,"Usage: writer [options]\n"
	  "-k hexadecimal shared memory key  (default: 40)\n"
	  "-p listening port number (default: %"PRIu64")\n"
	  "-i information port number (default: %"PRIu64")\n"
	  "-e ethernet device id (default: eth0)\n"
	  "-o print logging messages to stderr (as well as logfile)\n"
	  ,(uint64_t)WRITER_SERVICE_PORT,(uint64_t)WRITER_INFO_PORT);
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
  // big TODO is to get the information from the messenger process
  char* ascii_hdr = ipcbuf_get_next_write (hdu->header_block);
  short int station_id = getVDIFStationID (hdr);
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
  char *starthost = strstr(hostname, "difx");
  
  int state = STATE_STOPPED;
  uint64_t port = WRITER_SERVICE_PORT;
  uint64_t iport = WRITER_INFO_PORT;
  char cmd = CMD_NONE;
  int arg, maxsock = 0; 
  int stderr_output = 0;
  
  while ((arg = getopt(argc, argv, "hk:p:i:e:o")) != -1) {
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

  Connection c;
  c.sockoptval = 1; //release port immediately after closing connection

  Connection raw;
  raw.alen = sizeof(raw.rem_addr);
  raw.sockoptval = 32*1024*1024; //size of data socket internal buffer

  // TODO -- make this for passing source information
  Connection ci;
  ci.sockoptval = 1;

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


  dada_hdu_set_key (hdu,key);

  //Try to connect to existing data block
  if (dada_hdu_connect(hdu) < 0)
  {
    multilog (log, LOG_ERR, "Unable to connect to psrdada HDU.\n");
    exit (EXIT_FAILURE);
  }
  multilog (log, LOG_INFO, "Connected to psrdada HDU.\n");

  //create listening socket 
  if (serve (port, &c) < 0) {
    multilog (log, LOG_ERR,
      "Failed to create control socket on port %d.\n", port);
    exit (EXIT_FAILURE);
  }
  maxsock = c.rqst;

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
  setsockopt(raw.svc, SOL_SOCKET, SO_RCVTIMEO, (const char *)&(tv_1s), sizeof(tv_1s));

  // Open observation information socket
  if (serve (iport, &ci) < 0) {
    multilog (log, LOG_ERR,
      "Failed to create information socket on port %d.\n", iport);
    exit (EXIT_FAILURE);
  }

  int skip_frames = 0;
  int write_header = 0;
  uint64_t packets_written = 0;
  ssize_t info_bytes_read = 0;
  ssize_t raw_bytes_read = 0;
  ssize_t ipcio_bytes_written = 0;

  while(1) {

    if (state == STATE_STOPPED) {
      if (cmd == CMD_NONE)
        cmd = wait_for_cmd (&c, logfile_fp);

      if (cmd == CMD_START) 
      {
        state = STATE_STARTED;
        skip_frames = 0;
        packets_written = 0;
        // TODO -- block here until we get the observation info
        // TODO -- need to handle case of multiple observation documents
        info_bytes_read = read (ci.rqst,ci.buf,1024);
        multilog (log, LOG_INFO, "Read %d bytes from fd %d, expected %d.\n", 
            info_bytes_read,ci.rqst,sizeof(ObservationDocument));
        if (info_bytes_read != sizeof(ObservationDocument)) {
          // TODO will have to decide how to handle this error
          multilog (log, LOG_ERR, "Not a valid ObservationDocument!\n");
        }
        else {
          memcpy (&od, ci.buf, sizeof(ObservationDocument));
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
      FD_SET (c.rqst, &readfds);
      FD_SET (raw.svc, &readfds);
      if (select (maxsock+1,&readfds,NULL,NULL,&tv_500mus) < 0)
      {
        multilog (log, LOG_ERR, "[STATE_STARTED] Error calling select.");
        fflush (logfile_fp);
        exit (EXIT_FAILURE);
      }
      
      //if input is waiting on listening socket, read it
      if (FD_ISSET (c.rqst,&readfds)) {
        cmd = wait_for_cmd (&c, logfile_fp);
      }

      if (cmd == CMD_EVENT) {
        // NB -- am ignoring this logic branch for now, hence below printf
        printf("in CMD_EVENT, should not be here\n");
        //close data block
        if (dada_hdu_unlock_write(hdu) < 0) {
          fprintf(stderr,"unable to unlock psrdada HDU, exiting\n");
          exit(EXIT_FAILURE);
        }
        //fprintf(stderr,"after dada_hdu_unlock_write\n");
        
        //dump DB to file
        currt = time(NULL);
        localtime_r(&currt,&tmpt);
        strftime(currt_string,sizeof(currt_string), "%Y%m%d_%H%M%S", &tmpt);
        *(currt_string+15) = 0;
        char eventfile[128];
        snprintf(eventfile, 128,
            "%s/%s%s_%s_ev.out",EVENTDIR,starthost,dev,currt_string);
        
        FILE* evfd;
        if((evfd = fopen(eventfile,"w")) == NULL) {
          fprintf(stderr,
              "Writer: Could not open file %s for writing.\n",eventfile);
        }
        else {
          event_to_file(hdu->data_block,evfd);
          fclose(evfd);
        }

        //Get pending commands, resume writing to ring buffer if there are none
        FD_ZERO (&readfds);
        FD_SET (c.rqst, &readfds);
        FD_SET (raw.svc, &readfds);
        if (select (maxsock+1,&readfds,NULL,NULL,&tv_500mus) < 0)
        {
          fprintf(stderr,"Error calling select.");
          exit(1);
        }

        state = STATE_STOPPED;
        if (FD_ISSET (c.rqst,&readfds)) {
          cmd = wait_for_cmd (&c, logfile_fp);
          fprintf(stderr,
              "Writer: flushed out command socket after event_to_file.\n");
        }
        else
          cmd = CMD_START; 
        
        continue; 
      }
      
      //CMD_STOP --> change state to STOPPED, close data block
      else if (cmd == CMD_STOP) {
        state = STATE_STOPPED;
        skip_frames = 0;
        if (dada_hdu_unlock_write (hdu) < 0) {
          multilog (log, LOG_ERR,
              "Writer: unable to unlock psrdada HDU, exiting\n");
          exit_status = EXIT_FAILURE;
          break;
        }
        multilog (log, LOG_INFO, 
            "Wrote %d packets to psrdada buffer.\n",packets_written);
        packets_written = 0;
        cmd = CMD_NONE;
        fflush (logfile_fp);
        continue;
      }

      //CMD_QUIT --> close data block, shutdown listening socket, return
      else if (cmd == CMD_QUIT) {
        multilog (log, LOG_INFO,
            "[STATE_STARTED] received CMD_QUIT, exiting.\n");
        multilog (log, LOG_INFO,
            "dada_hdu_unlock_write result = %d.\n",
            dada_hdu_unlock_write (hdu));
        break;
      }

      else if (cmd == CMD_START) {
        multilog (log, LOG_INFO,
            "[STATE_STARTED] ignored CMD_START.\n");
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

            if (skip_frames <  25600*2) {
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
  
  shutdown(c.rqst,2);
  multilog (log, LOG_INFO,
      "dada_hdu_disconnect result = %d.\n", dada_hdu_disconnect (hdu));
  fclose (logfile_fp);
  //change_file_owner (logfile_fp, "vliteops");

  return exit_status;
} // end main

