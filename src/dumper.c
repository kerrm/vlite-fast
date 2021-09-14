// Writer is attached to the buffer we want to read from, so it is a
// natural source of the (mapped) memory address of the system memory.
//
// So we keep writer "in the loop" for the triggering process.  It
// catches the triggers and determines which addresses to write.  It
// then issues a vcopier_t struct with e.g. the filename, address, and
// other info for a 1-s dump and puts it on the socket.
//
// When copier receives it, it copies the data into its own memory
// space (if we're not zero copying or sharing, this can just be a
// malloc).  Perhaps at this point we can use mmap to fault it to
// disk?  Need to look into that.
//
// From some research, zero copy is still tough.  (In a newer kernel,
// copy_file_range is the way to go.)  Moreover, I'm worried about the
// SSD write speed itself.  So let's just make a standard copy for now
// and move along.  then issues an directly the memory addresses to
// dump over the socket.

#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>

#include "multilog.h"
#include "def.h"
#include "utils.h"
#include "multicast.h"

#define MAX_DUMP_REQUESTS 32
//#define DUMP_BUFSZ 64
#define DUMP_BUFSZ (2*VDIF_FRAME_SIZE*FRAMESPERSEC)

static multilog_t* mlog = NULL;
static FILE* logfile_fp = NULL;
static int mc_control_sock = 0;
static int mc_dumper_sock = 0;
char* mem_slots[MAX_DUMP_REQUESTS] = {NULL};
static dump_req_t reqs[MAX_DUMP_REQUESTS];
void cleanup (int status, int quit);

void send_quit () {
  const char cmdstop[] = {CMD_STOP};
  MulticastSend (mc_vlitegrp, MC_WRITER_PORT, cmdstop, 1);
}

void clear_req (dump_req_t* req) {
  req->buf_in=NULL;
  req->buf_local=NULL;
  req->fname[0] = '\0';
}

/* Move existing requests to beginning of queue, and return the index
 * of the first empty slot.
 */
int compress_reqs () {
  int first_empty = 0;
  // Move all existing dumps to the beginning of the queue
  // This is dumb double loop, but whatever
  for (int ireq1=0; ireq1 < MAX_DUMP_REQUESTS; ++ireq1) {
    // Is this slot empty?  If not, don't overwrite.
    if (reqs[ireq1].buf_in != NULL) {
      ++first_empty;
      continue;
    }
    for (int ireq2=ireq1; ireq2 < MAX_DUMP_REQUESTS; ++ireq2) {
      // Is this slot active?  If so, move it into the empty slot.
      if (reqs[ireq2].buf_in != NULL) {
        reqs[ireq1] = reqs[ireq2]; //full copy
        clear_req (reqs + ireq2);
        ++first_empty;
        break;
      }
    }
  }
  return first_empty;
}

/* Return the first unused local buffer.*/
char* get_local_buf () {
  for (int ireq1=0; ireq1 < MAX_DUMP_REQUESTS; ++ireq1) {
    int success = 1;
    for (int ireq2=0; ireq2 < MAX_DUMP_REQUESTS; ++ireq2) {
      if (reqs[ireq2].buf_local == mem_slots[ireq1]) {
        success = 0;
        break;
      }
    }
    if (success == 1)
      return mem_slots[ireq1];
  }
  return NULL;
}

void make_buffs() {
  // Make local memory for dumps
  for (int ireq=0; ireq < MAX_DUMP_REQUESTS; ++ireq) {
    mem_slots[ireq] = calloc (DUMP_BUFSZ, 1);
    if (mem_slots[ireq] == NULL) {
      // Failure.  Quit, and clean up memory
      multilog (mlog, LOG_ERR, "[DUMPER] Error allocating memory.\n");
      cleanup (-1, 1);
    }
  }
}

int addr_to_slot (char* addr) {
  for (int ireq=0; ireq < MAX_DUMP_REQUESTS; ++ireq) {
    if (addr==mem_slots[ireq])
      return ireq;
  }
  return -1;
}

void copy_voltages (dump_req_t* req) {
  if (req->bufsz != DUMP_BUFSZ) {
    multilog (mlog, LOG_ERR, "[DUMPER] Mismatched buffer size.\n");
    cleanup (-1, 1);
  }
  memcpy (req->buf_local, req->buf_in, req->bufsz);
}

void dump_voltages (dump_req_t* req)
{
  if (req->buf_local==NULL)
    return;
  printf ("Dumping from slot %d.\n",(addr_to_slot(req->buf_local)));
  // dump to disk
  // TODO -- check out O_DIRECT option
  int fd = open (req->fname, O_WRONLY | O_CREAT, 0664);
  if (fd < 0) {
    multilog (mlog, LOG_ERR, "[DUMPER] Could not open file %s for writing.\n",req->fname);
    cleanup (5,1);
  }

  uint64_t written = 0;
  while (written < req->bufsz)
  {
    ssize_t result = write (fd, req->buf_local+written, req->bufsz-written);
    if (result < 0)
    {
      if (EINTR != errno)
      {
        multilog (mlog, LOG_ERR, "[DUMPER] Could complete write to file %s, encountering errno %d.\n",req->fname,errno);
        cleanup (5,1);
      }
      // TODO -- need to sleep if interrupted?
    }
    else
      written += result;
  }
  fsync (fd);
  // needed since we don't reset umask, which is 0027
  fchmod (fd, 0664);
  close (fd);
  // hardcode for now vliteops=6362, vlite=5636
  chown (req->fname, 6362, 5636);
  req->buf_in = NULL;
  req->buf_local = NULL;
}

void open_log (char* hostname) {
  // logging
  gethostname (hostname,MAXHOSTNAME);
  mlog = multilog_open ("dumper",0);
  time_t currt = time (NULL);
  struct tm tmpt;
  localtime_r (&currt,&tmpt);
  char currt_string[128];
  strftime (currt_string,sizeof(currt_string), "%Y%m%d_%H%M%S", &tmpt);
  currt_string[15] = 0;
  char logfile[128];
  pid_t pid = getpid ();
  snprintf (logfile,128,
      "%s/%s_%s_dumper_%06d.log",LOGDIR,currt_string,hostname,pid);
  logfile_fp = fopen (logfile, "w");
  multilog_add (mlog, logfile_fp);
  //if (stderr_output)
  //multilog_add (mlog, stderr);
  //printf("writing log to %s\n",logfile);
}

void cleanup (int status, int quit) {
  if (quit!=0)
    send_quit();
  if (mc_control_sock > 0)
    shutdown (mc_control_sock, 2);
  if (mc_dumper_sock > 0)
    shutdown (mc_dumper_sock, 2);
  // Dump any remaining requests
  for (int ireq=0; ireq < MAX_DUMP_REQUESTS; ++ireq)
    dump_voltages (&reqs[ireq]);
  for (int ireq=0; ireq < MAX_DUMP_REQUESTS; ++ireq)
    free (mem_slots[ireq]);
  multilog (mlog, LOG_INFO, "Completed cleanup [DUMPER].\n");
  multilog_close (mlog);
  if (logfile_fp) {
    fflush (logfile_fp);
    fclose (logfile_fp);
  }
  exit (status);
}

void sigint_handler (int dummy) {
  cleanup (0, 0);
}



int main(int argc, char** argv) {

  // register SIGINT handling
  signal (SIGINT, sigint_handler);

  // register exit function
  //atexit (exit_handler);

  // TODO -- make this a utility function with "dumper" as arg etc.?
  char hostname[MAXHOSTNAME];
  open_log (hostname);
  
  // Multi-cast sockets: control and dump requests.
  // For simplicity, just subscribe to the "writer" group.  Why not just
  // use a single control group for all processes?
  int maxsock = 0;
  mc_control_sock = open_mc_socket (mc_vlitegrp, MC_WRITER_PORT,
      "Control Socket [Dumper]", &maxsock, mlog);
  char mc_from[24]; //ip address of multicast sender
  char mc_control_buf[32];
  
  // Listen on a specific group which should *only* have dump requests.
  mc_dumper_sock = open_mc_socket (mc_vlitegrp, MC_DUMPER_PORT,
      "Dump Socket [Dumper]", &maxsock, mlog);

  // Dump request queue
  for (int ireq=0; ireq < MAX_DUMP_REQUESTS; ++ireq)
    clear_req (reqs + ireq);

  // Set up local memory for voltage copies.
  make_buffs ();

  /*
  // TMP
  for (int ireq = 0; ireq < MAX_DUMP_REQUESTS; ++ireq) {
    sprintf (mem_slots[ireq], "%d test test test %d\n", ireq, ireq);
  }
  */

  fd_set readfds;
  struct timeval select_timeout; //timeout for select(), 10ms
  select_timeout.tv_sec = 0;
  select_timeout.tv_usec = 10000;

  /*
  // TMP
  int test_send = 0;
  dump_req_t dump = {.buf_in=mem_slots[0],.buf_local=NULL,.bufsz=DUMP_BUFSZ };
  char tmp[64];
  // end TMP
  */
  while (1) {

    /*
    // Send a dump request
    if (test_send < 64) {
      for (int idummy=0; idummy<4; ++idummy) {
        sprintf (tmp, "test_dump_%02d", test_send);
        dump.buf_in = mem_slots[16+test_send%16];
        strcpy (dump.fname, tmp);
        if (MulticastSend (mc_vlitegrp, MC_DUMPER_PORT, (const char *)(&dump), sizeof(dump_req_t)) < 0)
        {
          multilog (mlog, LOG_ERR, "send: %s\n", strerror (errno));
          cleanup (-1, 0);
        }
        test_send++;
      }
      sleep (0.01);
    }
    else {
      cleanup (0, 0);
    }
    */

    FD_ZERO (&readfds);
    FD_SET (mc_control_sock, &readfds);
    FD_SET (mc_dumper_sock, &readfds);
    if (select (maxsock+1,&readfds,NULL,NULL,&select_timeout) < 0)
    {
      multilog (mlog, LOG_ERR, "[DUMPER] Error calling select.\n");
      cleanup (4,1);
    }
    
    // Only command we are accepting now is quit!
    if (FD_ISSET (mc_control_sock, &readfds)) {
      char cmd = check_for_cmd (
          mc_control_sock, mc_control_buf, 32, logfile_fp);
      if (cmd == CMD_QUIT)
        cleanup (0,0);
    }

    // Check for dump requests.
    if (FD_ISSET (mc_dumper_sock, &readfds))
    {
      int first_empty = compress_reqs ();
      int slots = MAX_DUMP_REQUESTS - first_empty;
      if (slots==0) {
        multilog (mlog, LOG_INFO, "Dump request queue full!  Aborting.\n");
        cleanup (1,1);
      }
      //printf ("first_empty=%d\n",first_empty);
      //printf ("slots=%d\n",slots);

      // Read as many dump requests as will fit in the buffer.
      int new_dumps = 0;
      for (int ireq=0; ireq < 10*MAX_DUMP_REQUESTS; ++ireq) {
        if (first_empty==MAX_DUMP_REQUESTS) {
          multilog (mlog, LOG_INFO, "Out of slots! Aborting.\n");
          cleanup (-1,1);
        }
        int nbytes = MultiCastReceive (mc_dumper_sock, 
            (char*)(&reqs[first_empty]), sizeof(dump_req_t), mc_from);
        if (nbytes == 0)
          break;
        if (nbytes == -1)  {
          if ((EWOULDBLOCK == errno) || (EAGAIN == errno))
            break;
        }
        if (nbytes != sizeof(dump_req_t)) {
          // TODO this is an error
          multilog (mlog, LOG_INFO, "Received a corruped dump request.  Aborting.\n");
          cleanup (2,1);
        }
        if (strcmp (hostname, reqs[first_empty].hostname) != 0)
          continue;
        new_dumps++;
        first_empty++;
      } // end dump socket loop
      multilog (mlog, LOG_INFO, "Received %d dump requests.\n",new_dumps);
    }

    // Copy any new dump requests to local memory.  These are sorted such
    // that the earliest dumps are the oldest.
    int pending_writes = 0;
    for (int ireq=0; ireq < MAX_DUMP_REQUESTS; ++ireq) {

      // No request waiting
      if (reqs[ireq].buf_in == NULL)
        continue; // in principle could break because they are sorted

      // Already copied.
      if ((reqs[ireq].buf_in != NULL) && (reqs[ireq].buf_local != NULL)) {
        pending_writes++;
        continue;
      }
      
      // Need to copy.
      reqs[ireq].buf_local = get_local_buf();
      if (reqs[ireq].buf_local == NULL) {
        multilog (mlog, LOG_INFO, "No local memory. Aborting.\n");
        send_quit ();
        cleanup (3,1);
      }
      copy_voltages (reqs+ireq);
    }

    if (pending_writes == 0)
      continue;

    multilog (mlog, LOG_INFO, "Pending %d writes.\n", pending_writes);

    // Finally, write one (and only one) request to disk on each loop.
    for (int ireq=0; ireq < MAX_DUMP_REQUESTS; ++ireq) {
      if (reqs[ireq].buf_local != NULL) {
        dump_voltages (reqs+ireq);
        break;
      }
    }

  }

  cleanup (0,0);
  return 0; // will never get here
} // end main
