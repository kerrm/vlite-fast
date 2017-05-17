#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <signal.h>
#include <time.h>

#include "utils.h"
#include "multicast.h"
#include "def.h"
#include "executor.h"
#include "vlite_xml.h"
#include "alert.h"

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

// for signal handling and graceful termination
static volatile int keep_running = 1;

//How many Reader/Writer pairs are there
#define NRWPAIRS 1
#define HOST0 0 //read data from difxN hosts where N >= HOST0

void usage ()
{
  //XXX: add options to connect to writers/readers at certain ports?
  fprintf(stdout,"Usage: new_messenger [options]\n"
	  "-m enable manual control on specified port (e.g. via telnet)\n"
  );
}

void sigint_handler (int dummy) {
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

int main(int argc, char** argv)
{
  // register SIGINT handling
  signal (SIGINT, sigint_handler);

  uint64_t manual_control_port = 0;

  int arg = 0;
  while ((arg = getopt(argc, argv, "hm:")) != -1) {
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
    }  
  }

  const char cmdstop[] = {CMD_STOP};
  const char cmdevent[] = {CMD_EVENT};
  const char cmdquit[] = {CMD_QUIT};
  const char cmdstart[] = {CMD_START};

  //XXX: in final version, there will be an array of readers and an array of writers, so initial connection attempts as well as each command send will loop over those. (How to get port numbers for all? Specify on command line?)  
  Connection cr[NRWPAIRS];
  Connection cw[NRWPAIRS];
  Connection ci[NRWPAIRS];
  char hostname[MAXHOSTNAME]; // = "vlite-difx1.evla.nrao.edu";
  fd_set readfds;
  int obsinfosock, antpropsock, alertsock, maxnsock = 0, nbytes, ii, nhost = HOST0-1;
  int ant_on_src = 0;
  int recording = 0;

  ScanInfoDocument D; //multicast message struct
  const ObservationDocument *od;
  const AntPropDocument *ap;
  AlertDocument A; 

  struct timespec ts_100ms;
  ts_100ms.tv_sec = 0;
  ts_100ms.tv_nsec = 100000000;

  struct timespec ts_10ms;
  ts_10ms.tv_sec = 0;
  ts_10ms.tv_nsec = 10000000;

  struct timespec ts_1ms;
  ts_1ms.tv_sec = 0;
  ts_1ms.tv_nsec = 1000000;

  char scaninfofile[128];
  char eventlogfile[] = "eventlog.txt";
  char msg[MSGMAXSIZE];
  char src[SRCMAXSIZE];
  char from[24]; //ip address of multicast sender

  FILE *sfd, *efd;

  //Initialize Connections and connect to Readers/Writers: two pairs per difx host
  for(ii=0; ii<NRWPAIRS; ii++) {
    if(ii%2 == 0) {
      nhost++;
      cr[ii].port = READER_SERVICE_PORT;
      cw[ii].port = WRITER_SERVICE_PORT;
      ci[ii].port = WRITER_SERVICE_PORT + 21;
    }
    else {
      cr[ii].port = READER_SERVICE_PORT + 1;
      cw[ii].port = WRITER_SERVICE_PORT + 1;
      ci[ii].port = WRITER_SERVICE_PORT + 21 + 1;
    }
    
    sprintf(cr[ii].hostname,"vlite-difx%d.evla.nrao.edu",nhost+1);
    sprintf(cw[ii].hostname,"vlite-difx%d.evla.nrao.edu",nhost+1);
    sprintf(ci[ii].hostname,"vlite-difx%d.evla.nrao.edu",nhost+1);
    cr[ii].isconnected = 0;
    cw[ii].isconnected = 0;
    ci[ii].isconnected = 0;

    printf("Reader on %s, port %d\n", cr[ii].hostname, cr[ii].port);
    printf("Writer on %s, port %d\n", cw[ii].hostname, cw[ii].port);
    printf("Writer info on %s, port %d\n", ci[ii].hostname, ci[ii].port);

    if(conn(&cw[ii]) < 0) {
      fprintf(stderr,"Messenger: could not connect to Writer on %s port %d\n", cw[ii].hostname, cw[ii].port);
      exit(1);
    }
    else
      cw[ii].isconnected = 1;

    // wait 10 ms for connection to happen and writer to begin listening on next port
    nanosleep(&ts_100ms,NULL);

    if(conn(&ci[ii]) < 0) {
      fprintf(stderr,"Messenger: could not connect to Writer info on %s port %d\n", ci[ii].hostname, ci[ii].port);
      exit(1);
    }
    else
      ci[ii].isconnected = 1;
    
    if(conn(&cr[ii]) < 0 ) {
      fprintf(stderr,"Messenger: could not connect to Reader on %s port %d\n", cr[ii].hostname, cr[ii].port);
      exit(1);
    }
    else
      cr[ii].isconnected = 1;
  }

  // create manual control socket
  Connection mc;
  if (manual_control_port > 0) {
    mc.port = manual_control_port;
    if (serve (manual_control_port, &mc) < 0) {
      //multilog (log, LOG_ERR,
      fprintf (stderr,
        "Failed to create manual control socket on port %d.\n", manual_control_port);
      exit (EXIT_FAILURE);
    }
    if (mc.rqst > maxnsock)
      maxnsock = mc.rqst;
  }

  //connect to VLA obsinfo multicast
  obsinfosock = openMultiCastSocket(obsinfogrp, MULTI_OBSINFO_PORT);
  if (obsinfosock < 0) {
    fprintf(stderr,"Failed to open Observation multicast; openMultiCastSocket = %d\n",obsinfosock);
    exit (EXIT_FAILURE);
  }
  else {
    printf ("Obsinfo socket: %d\n",obsinfosock);
    if(obsinfosock > maxnsock)
      maxnsock = obsinfosock;
  }

  //connect to VLA antprop multicast
  antpropsock = openMultiCastSocket(antpropgrp, MULTI_ANTPROP_PORT);
  if (antpropsock < 0) {
    fprintf (stderr,"Failed to open Antprop multicast; openMultiCastSocket = %d\n",antpropsock);
    exit (EXIT_FAILURE);
  }
  else {
    printf ("Antprop socket: %d\n",antpropsock);
    if (antpropsock > maxnsock)
      maxnsock = antpropsock;
  }

  //connect to VLA alert multicast
  alertsock = openMultiCastSocket (alertgrp, MULTI_ALERT_PORT);
  if (alertsock < 0) {
    fprintf(stderr,"Failed to open Alert multicast; openMultiCastSocket = %d\n",alertsock);
    exit(EXIT_FAILURE);
  }
  else {
    printf ("Alert socket: %d\n",alertsock);
    if (alertsock > maxnsock)
      maxnsock = alertsock;
  }

  //Open event log
  if ((efd = fopen(eventlogfile, "a")) == NULL) {
    fprintf(stderr,"Messenger: Could not open event log file %s for writing.\n",eventlogfile);
    exit(EXIT_FAILURE);
  }

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
      FD_SET (mc.rqst, &readfds);
    else {
      FD_SET (obsinfosock, &readfds);
      FD_SET (antpropsock, &readfds);
      FD_SET (alertsock, &readfds);
    }

    //printf("Waiting for multicast messages...\n");
    select (maxnsock+1,&readfds,NULL,NULL,NULL);
    
    //Obsinfo socket
    if (FD_ISSET (obsinfosock,&readfds)) {
      nbytes = MultiCastReceive(obsinfosock, msg, MSGMAXSIZE, from);
      //fprintf(stderr,"Received %d bytes from Obsinfo multicast.\n",nbytes);
      if(nbytes <= 0) {
        fprintf(stderr,
            "Error on Obsinfo socket, return value = %d\n",nbytes);
        continue;
      }
      
      parseScanInfoDocument(&D, msg);
      printScanInfoDocument(&D);
      printf("Message type: %d = %s\n",D.type,ScanInfoTypeString[D.type]);
      if(D.type == SCANINFO_OBSERVATION) {
        od = &(D.data.observation);
        strcpy(src,od->name);
	
        sprintf(scaninfofile,"%s/%s.%s.obsinfo.%04d.%04d.txt",OBSINFODIR,od->datasetId,od->name,od->scanNo,od->subscanNo);

        if((sfd = fopen(scaninfofile,"w")) == NULL) {
          fprintf(stderr, "Messenger: Could not open file %s for writing.\n",scaninfofile);
        }
        else {
          fprintScanInfoDocument(&D,sfd);
          fclose(sfd);
        }
      
        // TO ADD: Check that Writers are still connected before sending; change their isonnected elements if they are not
        if (strcasecmp (src,"FINISH") == 0) {
          for(ii=0; ii<NRWPAIRS; ii++) {
            printf("sending STOP command\n");
            if (send(cr[ii].svc, cmdstop, 1, 0) == -1)
              perror("send");
            if (send(cw[ii].svc, cmdstop, 1, 0) == -1)
              perror("send");
          }
          recording = 0;
        }
        else {
          // TO ADD: Check that Readers/Writers are still connected before sending; change their isonnected elements if they are not
          if (recording) {
            for(ii=0; ii<NRWPAIRS; ii++) {
              printf("sending STOP command\n");
              if (send(cw[ii].svc, cmdstop, 1, 0) == -1)
                perror("send");
            }
            // this is empirical, but seems to be long enough to allow
            // messages to propagate
            nanosleep(&ts_100ms, NULL);
          }
          for(ii=0; ii<NRWPAIRS; ii++) {
            printf("sending START command\n");
            if (send(cr[ii].svc, cmdstart, 1, 0) == -1)
              perror("send");
            if (send(cw[ii].svc, cmdstart, 1, 0) == -1)
              perror("send");
            if (send(ci[ii].svc, od, sizeof(ObservationDocument), 0) == -1)
              perror("send");
          }
          nanosleep(&ts_10ms, NULL);
          recording = 1;
        }
      }
    }

    //Antprop socket
    if (FD_ISSET (antpropsock,&readfds)) {
      nbytes = MultiCastReceive (antpropsock, msg, MSGMAXSIZE, from);
      //fprintf(stderr,"Received %d bytes from Antprop multicast.\n",nbytes);
      if(nbytes <= 0) {
        fprintf(stderr,
            "Error on Antprop socket, return value = %d\n",nbytes);  
        continue;
      }

      parseScanInfoDocument(&D, msg);
      printScanInfoDocument(&D);
      printf("Message type: %d = %s\n",D.type,ScanInfoTypeString[D.type]);
      if (D.type == SCANINFO_ANTPROP) {
        ap = &(D.data.antProp);

        sprintf (scaninfofile,"%s/%s.antprop.txt",OBSINFODIR,ap->datasetId);
        sfd = fopen (scaninfofile,"w");

        if ((sfd = fopen (scaninfofile,"w")) == NULL) {
          fprintf (stderr,"Messenger: Could not open file %s for writing.\n",scaninfofile);
        }
        else {
          fprintScanInfoDocument (&D,sfd);
          fclose (sfd);
        }
      }
    }

    //Alert socket
    if (FD_ISSET (alertsock,&readfds)) {
      nbytes = MultiCastReceive (alertsock, msg, MSGMAXSIZE, from);
      //fprintf(stderr,"Received %d bytes from Alert multicast.\n",nbytes);
      if (nbytes <= 0) {
        fprintf (stderr,"Error on Alert socket, return value = %d\n",nbytes);
        continue;
      }

      parseAlertDocument(&A, msg);
      //if((strcmp(A.monitorName,"ELPosError") == 0 || strcmp(A.monitorName,"AZPosError") == 0) && A.alertState == 1) {
      if(strcmp(A.monitorName,"ELPosError") == 0 || strcmp(A.monitorName,"AZPosError") == 0) {
      //if (1) {
        //printAlertDocument(&A);
        //printf("alertState = %d\n", A.alertState);
        //sprintf(scaninfofile,"%f.alert.txt",A.timeStamp);
        sprintf(scaninfofile,"%s/%f.alert.txt",OBSINFODIR,A.timeStamp);
        sfd = fopen(scaninfofile,"w");
        fprintAlertDocument(&A,sfd);
        fclose(sfd);
      }
    }

    // Manual control socket
    if (FD_ISSET (mc.rqst,&readfds)) {
      char cmd = wait_for_cmd (&mc);
      if (cmd == CMD_START) {
        ObservationDocument dummy;
        fill_dummy_obs_doc (&dummy);
        printf("sending manual START command\n");
        for(ii=0; ii<NRWPAIRS; ii++) {
          if (send(cr[ii].svc, cmdstart, 1, 0) == -1)
            perror("send");
          if (send(cw[ii].svc, cmdstart, 1, 0) == -1)
            perror("send");
          if (send(ci[ii].svc, &dummy, sizeof(ObservationDocument), 0) == -1)
            perror("send");
        }
      }
      else if (cmd == CMD_STOP) {
        printf("sending STOP command\n");
        for(ii=0; ii<NRWPAIRS; ii++) {
          if (send(cr[ii].svc, cmdstop, 1, 0) == -1)
            perror("send");
          if (send(cw[ii].svc, cmdstop, 1, 0) == -1)
            perror("send");
        }
      }

      else if (cmd == CMD_QUIT) {
        printf("breaking from main loop");
        break;
        /*
        printf("sending QUIT command\n");
        for(ii=0; ii<NRWPAIRS; ii++) {
          if (send(cr[ii].svc, cmdquit, 1, 0) == -1)
            perror("send");
          if (send(cw[ii].svc, cmdquit, 1, 0) == -1)
            perror("send");
        }
        */
      }
    } // end manual control logic

    /*
    if(FD_ISSET(heimdall.rqst,&readfds)) {
      //read() or recv()--form of message TBD
      //send EVENT to Writer
    }
    */
    
  } //end while

  printf("Shutting down...\n");

  for(ii=0; ii<NRWPAIRS; ii++) {
    if (send(cr[ii].svc, cmdquit, 1, 0) == -1)
      perror("send");
  }
  // run writers for a little longer to avoid them "running out" of packets before
  // processing signal
  sleep (2);
  for(ii=0; ii<NRWPAIRS; ii++) {
    if (send(cw[ii].svc, cmdquit, 1, 0) == -1)
      perror("send");
  }

  fclose(efd);
}
