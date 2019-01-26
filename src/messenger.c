#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <signal.h>
#include <time.h>
#include <errno.h>
#include <math.h>

#include "utils.h"
#include "multicast.h"
#include "def.h"
#include "executor.h"
#include "vlite_xml.h"
#include "alert.h"

#include "multilog.h"

#define MSGMAXSIZE 8192

//Multicast ports
#define MULTI_OBSINFO_PORT 53001
#define MULTI_ANTPROP_PORT 53000
#define MULTI_ALERT_PORT 20011
#define MULTI_TEST_PORT 53901

//Multicast group IPs
char testgrp[] = "239.199.3.2";
char obsinfogrp[] = "239.192.3.2";
char antpropgrp[] = "239.192.3.1";
char alertgrp[] = "239.192.2.3";

int mc_reader_port = 20000;
int mc_writer_port = 20001;
int mc_info_port = 20002;

// for signal handling and graceful termination
static volatile int keep_running = 1;
static volatile int signal_received = 0;

void usage ()
{
  fprintf(stdout,"Usage: messenger configuration_file [options]\n"
	  "-m enable manual control on specified port (e.g. via telnet)\n"
	  "-o print logging messages to stdout (as well as logfile)\n"
  );
}

void sigint_handler (int dummy) {
  signal_received = dummy;
  keep_running = 0;
}

void fill_dummy_obs_doc (ObservationDocument* od)
{
  snprintf(od->name,OBSERVATION_NAME_SIZE,"ManualTest");
  od->ra = 1;
  od->dec = 1;
  od->startTime = 55555;
  /*
  datasetId = X_osro.57889.63195775463
  configId = X_osro.57889.63195775463.1
  startTime = 57889.63464120
  name = 3C84
  ra = 0.87180360
  dec = 0.72451575
  dra = 0.00000000
  ddec = 0.00000000
  azoffs = 0.00000000
  eloffs = 0.00000000
  startLST = 0.98718036
  scanNo = 12
  subscanNo = 1
  primaryBand = 10GHz
  usesPband = 0
  */
}

// just a version with built in error logging
void mc_send(const char *group, int port, const char *message,
    int length, multilog_t* log)
{
  if (MulticastSend (group, port, message, length) < 0)
    multilog (log, LOG_ERR, "send: %s\n", strerror (errno));
}
    

int main(int argc, char** argv)
{
  // register SIGINT handling
  //signal (SIGINT, sigint_handler);
  // TEMPORARY -- handle all signals with quit
  for (int i=0; i < 15; ++i)
    signal (i+1, sigint_handler);

  uint64_t manual_control_port = 0;
  int stdout_output = 1;
  int write_alerts = 0;
  int write_obsinfos = 1;

  int arg = 0;
  while ((arg = getopt(argc, argv, "hom:")) != -1) {
    switch(arg) {

    case 'h':
      usage ();
      return 0;
      
    case 'm':
      if (sscanf (optarg, "%"PRIu64"", &manual_control_port) != 1) {
        fprintf (stderr, "writer: could not parse port from %s\n", optarg);
        return -1;
      }
      break;

    case 'o':
      stdout_output = 1;
      break;

    } // end switch statement
  } // end argument parsing

  // check for mandatory configuration file
  if (optind==argc)
    fprintf (stderr, "Must specify configuration file!.\n");


  multilog_t* log = multilog_open ("messenger",0);
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
      "%s/%s_%s_messenger_%06d.log",LOGDIR,currt_string,hostname,pid);
  FILE *logfile_fp = fopen (logfile, "w");
  multilog_add (log, logfile_fp);
  if (stdout_output)
    multilog_add (log, stdout);

  int nvconf = 0;
  VFASTConfig** vconf = parse_vfast_config (argv[optind], &nvconf);
  if (NULL == vconf) {
    multilog (log, LOG_ERR, "Could not parse configuration file: %s!\n",
        argv[optind]);
    exit (EXIT_FAILURE);
  }
  multilog (log, LOG_INFO, "Loaded configuration from %s.\n",argv[optind]);
  for (int iconfig = 0; iconfig < nvconf; ++iconfig)
    multilog (log, LOG_INFO, "Configuration %d:\n%s", iconfig, 
        print_vfast_config (vconf[iconfig], NULL));

  const char cmdstop[] = {CMD_STOP};
  const char cmdevent[] = {CMD_EVENT};
  const char cmdquit[] = {CMD_QUIT};
  const char cmdstart[] = {CMD_START};

  fd_set readfds;
  int obsinfosock, antpropsock, alertsock, maxnsock = 0;
  //int ant_on_src = 0;
  int recording = 0;

  ScanInfoDocument D; //multicast message struct
  const ObservationDocument *od;
  AlertDocument A; 
  AntPropDocument last_antprop; // keep a copy of antenna properties

  struct timespec ts_500ms = get_ms_ts (500);
  struct timespec ts_100ms = get_ms_ts (100);
  struct timespec ts_10ms = get_ms_ts (10);

  char scaninfofile[128];
  //char eventlogfile[] = "eventlog.txt";
  char msg[MSGMAXSIZE];
  char from[24]; //ip address of multicast sender

  // keep track of previous pointing and combine observation if unchanged
  double last_ra = 100;
  double last_dec = 100;
  double max_integ = 480;
  time_t start_integ, stop_integ;

  FILE *sfd;//, *efd;

  // create manual control socket
  Connection mc;
  if (manual_control_port > 0) {
    mc.port = manual_control_port;
    if (serve (manual_control_port, &mc) < 0) {
      multilog (log, LOG_ERR,
        "Failed to create manual control socket on port %"PRIu64".\n", manual_control_port);
      exit (EXIT_FAILURE); // TODO -- cleaner shutdown?
    }
    if (mc.rqst > maxnsock)
      maxnsock = mc.rqst;
  }
  //

  //connect to VLA obsinfo multicast
  obsinfosock = openMultiCastSocket (obsinfogrp, MULTI_OBSINFO_PORT);
  if (obsinfosock < 0) {
    multilog (log, LOG_ERR, "Failed to open Observation multicast; openMultiCastSocket = %d\n",obsinfosock);
    exit (EXIT_FAILURE);
  }
  else {
    multilog (log, LOG_INFO, "Obsinfo socket: %d\n",obsinfosock);
    if(obsinfosock > maxnsock)
      maxnsock = obsinfosock;
  }

  //connect to VLA antprop multicast
  antpropsock = openMultiCastSocket (antpropgrp, MULTI_ANTPROP_PORT);
  if (antpropsock < 0) {
    multilog (log, LOG_ERR, "Failed to open Antprop multicast; openMultiCastSocket = %d\n",antpropsock);
    exit (EXIT_FAILURE);
  }
  else {
    multilog (log, LOG_INFO, "Antprop socket: %d\n",antpropsock);
    if (antpropsock > maxnsock)
      maxnsock = antpropsock;
  }

  //connect to VLA alert multicast
  alertsock = openMultiCastSocket (alertgrp, MULTI_ALERT_PORT);
  if (alertsock < 0) {
    multilog (log, LOG_ERR,"Failed to open Alert multicast; openMultiCastSocket = %d\n",alertsock);
    exit (EXIT_FAILURE);
  }
  else {
    multilog (log, LOG_INFO, "Alert socket: %d\n",alertsock);
    if (alertsock > maxnsock)
      maxnsock = alertsock;
  }

  /*
  //Open event log
  if ((efd = fopen (eventlogfile, "a")) == NULL) {
    multilog (log, LOG_ERR, "Messenger: Could not open event log file %s for writing.\n",eventlogfile);
    exit (EXIT_FAILURE);
  }
  */

  while (keep_running) {
    /* From Walter: The antprop message comes whenever a new scheduling 
     * block starts and at each UT midnight (where EOP values might change) 
     * Observation document happens before each new scan and a FINISH 
     * document at the end of a scheduling block. At the beginning of the 
     * scheduling block the antprop document always precedes the observation
     * document.
     */

    //Blocking select() on multicast sockets and Heimdall event trigger socket
    //Have to call FS_SET on each socket before calling select()
    FD_ZERO (&readfds);
    if (manual_control_port > 0)
    {
      FD_SET (mc.rqst, &readfds);
      FD_SET (obsinfosock, &readfds);
    }
    else {
      FD_SET (obsinfosock, &readfds);
      FD_SET (antpropsock, &readfds);
      FD_SET (alertsock, &readfds);
    }

    select (maxnsock+1,&readfds,NULL,NULL,NULL);
    fflush (stdout);
    fflush (stderr);
    
    //Obsinfo socket
    if (FD_ISSET (obsinfosock,&readfds)) {
      int nbytes = MultiCastReceive (obsinfosock, msg, MSGMAXSIZE, from);
      if (nbytes <= 0) {
        multilog (log, LOG_ERR,
            "Error on Obsinfo socket, return value = %d\n",nbytes);
        continue;
      }
      
      parseScanInfoDocument (&D, msg);
      printScanInfoDocument (&D);
      multilog (log, LOG_INFO,
          "Message type: %d = %s\n",D.type,ScanInfoTypeString[D.type]);

      // only do this part if not under manual control
      if (0==manual_control_port)
      {
        if (D.type == SCANINFO_OBSERVATION) {
          od = &(D.data.observation);
    
          if (write_obsinfos) {
            sprintf (scaninfofile,"%s/%s.%s.obsinfo.%04d.%04d.txt",
              OBSINFODIR,od->datasetId,od->name,od->scanNo,od->subscanNo);

            if ((sfd = fopen(scaninfofile,"w")) == NULL) {
              multilog (log, LOG_ERR, 
                "Messenger: Could not open file %s for writing.\n",scaninfofile);
            }
            else {
              fprintScanInfoDocument (&D,sfd);
              fclose (sfd);
            }
          }
        
          // TO ADD: Check that Writers are still connected before sending; change their isonnected elements if they are not
          if (strcasecmp (od->name,"FINISH") == 0)
          {
            multilog (log, LOG_INFO, "sending STOP command\n");
            mc_send ("224.3.29.71",mc_reader_port,cmdstop,1,log);
            mc_send ("224.3.29.71",mc_writer_port,cmdstop,1,log);
            recording = 0;
            // TEMP -- see if this sleep is helpful
            nanosleep (&ts_500ms, NULL);
          }
          else if (od->scanNo == 1)
            multilog (log, LOG_INFO, "scanNo==1: ignoring DUMMY scan\n");
          else
          {

            // TODO: Check that Readers/Writers are still connected
            // before sending; change their isonnected elements if they 
            // are not

            if (recording) {
              // check to see if the pointing position has changed; tweak
              // this to allow for small changes to keep integrating during
              // VLASS
              if ((fabs(od->ra-last_ra)< 0.5) && fabs(od->dec-last_dec) < 0.5)
              {
                multilog (log, LOG_INFO,
                    "Pointing unchanged, will continue to integrate.\n");
                last_ra = od->ra;
                last_dec = od->dec;
                time (&stop_integ);
                double dt = difftime (stop_integ, start_integ);
                if (dt < max_integ)
                  continue;
                else
                  multilog (log, LOG_INFO,
                      "Integration exceeded %lf s, starting a new one.\n",
                      max_integ);
              }

              multilog (log, LOG_INFO, "sending STOP command\n");
              mc_send ("224.3.29.71",mc_writer_port,cmdstop,1,log);
              // this is empirical, but seems to be long enough to allow
              // messages to propagate before sending the next command;
              // can also avoid this by altering wait_for_cmd to not discard
              // the older commands
              nanosleep (&ts_500ms, NULL);
            }

            // update pointing position and reset integration timer
            last_ra = od->ra;
            last_dec = od->dec;
            time (&start_integ);
            multilog (log, LOG_INFO,"sending START command\n");

            mc_send ("224.3.29.71",mc_reader_port,cmdstart,1,log);
            mc_send ("224.3.29.71",mc_writer_port,cmdstart,1,log);
            mc_send ("224.3.29.71",mc_info_port,(char*)(od),sizeof(ObservationDocument),log);
            multilog (log, LOG_INFO,"sent START command\n");
            fflush (stderr); fflush (stdout);
            nanosleep (&ts_10ms, NULL);
            recording = 1;
          }
        }
      }
    }

    //Antprop socket
    if (FD_ISSET (antpropsock,&readfds)) {
      int nbytes = MultiCastReceive (antpropsock, msg, MSGMAXSIZE, from);
      if (nbytes <= 0) {
        multilog (log, LOG_ERR,
            "Error on Antprop socket, return value = %d\n",nbytes);  
        continue;
      }

      parseScanInfoDocument (&D, msg);
      printScanInfoDocument (&D);
      multilog (log, LOG_INFO,
          "Message type: %d = %s\n",D.type,ScanInfoTypeString[D.type]);
      if (D.type == SCANINFO_ANTPROP) {

        const AntPropDocument* ap = &(D.data.antProp);
        // copy
        last_antprop = *ap;

        sprintf (scaninfofile,"%s/%s.antprop.txt",OBSINFODIR,ap->datasetId);
        sfd = fopen (scaninfofile,"w");

        if ((sfd = fopen (scaninfofile,"w")) == NULL) {
          multilog (log, LOG_ERR,
            "Messenger: Could not open file %s for writing.\n",scaninfofile);
        }
        else {
          fprintScanInfoDocument (&D,sfd);
          fclose (sfd);
        }
      }
    } // end Antprop socket block

    //Alert socket
    if (FD_ISSET (alertsock,&readfds)) {
      int nbytes = MultiCastReceive (alertsock, msg, MSGMAXSIZE, from);
      if (nbytes <= 0) {
        multilog (log, LOG_ERR,
            "Error on Alert socket, return value = %d\n",nbytes);
        continue;
      }

      parseAlertDocument (&A, msg);
      if (write_alerts && (strcmp (A.monitorName,"ELPosError") == 0 || strcmp (A.monitorName,"AZPosError") == 0)) {
        //printAlertDocument(&A);
        //printf("alertState = %d\n", A.alertState);
        //sprintf(scaninfofile,"%f.alert.txt",A.timeStamp);
        sprintf (scaninfofile,"%s/%f.alert.txt",OBSINFODIR,A.timeStamp);
        sfd = fopen (scaninfofile,"w");
        fprintAlertDocument (&A,sfd);
        fclose (sfd);
      }
    } // end Alert socket block

    // Manual control socket
    if ((manual_control_port > 0) && FD_ISSET (mc.rqst,&readfds)) {
    //if (FD_ISSET (mc.rqst,&readfds)) {
      char cmd = wait_for_cmd (&mc, NULL);
      if (cmd == CMD_START) {
        ObservationDocument dummy;
        fill_dummy_obs_doc (&dummy);
        multilog (log, LOG_INFO, "sending manual START command\n");
        mc_send ("224.3.29.71",mc_reader_port,cmdstart,1,log);
        mc_send ("224.3.29.71",mc_writer_port,cmdstart,1,log);
        mc_send ("224.3.29.71",mc_info_port,(char *)(&dummy),sizeof(ObservationDocument),log);
      }
      else if (cmd == CMD_EVENT) {
        multilog (log, LOG_INFO, "sending manual EVENT command\n");
        mc_send ("224.3.29.71",mc_reader_port,cmdstop,1,log);
        mc_send ("224.3.29.71",mc_writer_port,cmdevent,1,log);
      }
      else if (cmd == CMD_STOP) {
        multilog (log, LOG_INFO, "sending STOP command\n");
        mc_send ("224.3.29.71",mc_reader_port,cmdstop,1,log);
        mc_send ("224.3.29.71",mc_writer_port,cmdstop,1,log);
      }
      else if (cmd == CMD_QUIT) {
        multilog (log, LOG_INFO, "breaking from main loop\n");
        break;
      }
    } // end manual control block

    // TODO -- could have messages from transient detection here
    
  } //end main program loop

  multilog (log, LOG_INFO, "Shutting down...\n");
  if (signal_received != 0)
    multilog (log, LOG_ERR, "received signal %d\n", signal_received+1);
  fclose (logfile_fp);

  // shut down readers
  mc_send ("224.3.29.71",mc_reader_port,cmdquit,1,log);
  // run writers for a little longer to avoid readers hanging
  sleep (2);
  mc_send ("224.3.29.71",mc_writer_port,cmdquit,1,log);

  // clean up configuration memory
  for (int ii=0; ii < nvconf; ++ii)
    free (vconf[ii]);
  //fclose(efd);
}
