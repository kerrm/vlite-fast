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
#include <sys/mman.h>
#include <sys/shm.h>

#include "dada_def.h"
#include "dada_hdu.h"
#include "ipcio.h"
#include "ascii_header.h"
#include "multilog.h"
#include "multilog.h"
#include "def.h"
#include "utils.h"
#include "multicast.h"

#define MAX_DUMP_REQUESTS 40
#define DUMP_HISTORY 180
#define DUMP_BUFSZ (2*VDIF_FRAME_SIZE*FRAMESPERSEC)

static multilog_t* mlog = NULL;
static dada_hdu_t* hdu = NULL;
static FILE* logfile_fp = NULL;
static char currt_string[64];
static char dump_fname[256];
static int mc_control_sock = 0;
static int mc_trigger_sock = 0;
static char* mem_slots[MAX_DUMP_REQUESTS] = {NULL};
static dump_req_t reqs[MAX_DUMP_REQUESTS];
static trigger_t trigs[MAX_TRIGGERS];
static time_t working_times [MAX_DUMP_REQUESTS] = {0};
static char* working_bufs [MAX_DUMP_REQUESTS] = {0};
static time_t dumped_times[DUMP_HISTORY] = {0};

void cleanup (int status, int quit);
void clear_req (dump_req_t* req);
void dumpchk (int cond, const char* msg);
void dump_voltages (dump_req_t* req);
void send_quit ();

void usage ()
{
  fprintf(stdout,"Usage: dumper [options]\n"
	  "-k hexadecimal shared memory key  (default: 40)\n"
	  "-o print logging messages to stderr (as well as logfile)\n"
    );
}

void cleanup (int status, int quit) {
  multilog (mlog, LOG_INFO, "[DUMPER] Beginning cleanup.\n");
  if (quit!=0)
    send_quit();
  if (hdu != NULL) {
    multilog (mlog, LOG_INFO, "[DUMPER] Beginning HDU cleanup.\n");
    fflush (logfile_fp);
    fflush (stderr);
    dada_hdu_disconnect (hdu);
    dada_hdu_destroy (hdu);
  }
  if (mc_control_sock > 0)
    shutdown (mc_control_sock, 2);
  if (mc_trigger_sock > 0)
    shutdown (mc_trigger_sock, 2);
  // Dump any remaining requests
  for (int ireq=0; ireq < MAX_DUMP_REQUESTS; ++ireq) {
    dump_voltages (&reqs[ireq]);
    clear_req (&reqs[ireq]);
  }
  for (int ireq=0; ireq < MAX_DUMP_REQUESTS; ++ireq) {
    munlock (mem_slots[ireq],DUMP_BUFSZ);
    free (mem_slots[ireq]);
  }
  multilog (mlog, LOG_INFO, "Completed cleanup [DUMPER].\n");
  multilog_close (mlog);
  if (logfile_fp) {
    fflush (logfile_fp);
    fclose (logfile_fp);
  }
  exit (status);
}


void send_quit () {
  const char cmdstop[] = {CMD_QUIT};
  MulticastSend (mc_vlitegrp, MC_WRITER_PORT, cmdstop, 1);
}

void clear_req (dump_req_t* req) {
  req->buf_in = NULL;
  req->buf_local=NULL;
  req->bufsz=0;
  req->utc_epoch = 0;
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

int ttcmp (const void *p, const void *q) {
  time_t x = *(const time_t *)p;
  time_t y = *(const time_t *)q;
  if (x < y) return -1;
  if (x > y) return 1;
  return 0;
}

/* Examine the trigger buffer, determine which voltage buffers overlap
 * the event times, and populate the "working" buffers with the times and
 * memory addresses to write out.  Return the number of seconds to dump.
 */
int trigs_to_utcs (int ntrig) {

  if (ntrig==0)
    return 0;

  ipcbuf_t* buf = (ipcbuf_t *) hdu->data_block;
  int nbufs = ipcbuf_get_nbufs (buf);
  int nbuf_to_write = 0;
  char** bufs = buf->buffer;
  double min_time = 1e100;
  double max_time = -1e100;
  double min_trig = 1e100;
  double max_trig = -1e100;

  for (int itrig=0; itrig < ntrig; ++itrig) {
    if (trigs[itrig].t0 < min_trig)
      min_trig = trigs[itrig].t0;
    if (trigs[itrig].t1 > max_trig)
      max_trig = trigs[itrig].t1;
  }
  multilog (mlog, LOG_INFO, "earliest trigger t0 = %.3f\n",min_trig);
  multilog (mlog, LOG_INFO, "latest   trigger t0 = %.3f\n",max_trig);

  for (int ibuf=0; ibuf < nbufs; ++ibuf)
  {
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
    for (int itrig=0; itrig < ntrig; ++itrig) {
      if ((trigs[itrig].t0 <= (buf_t0+1)) && (trigs[itrig].t1 >= buf_t0)) {
        trigger_overlap = 1;
        break;
      }
    } // end loop over trigger times

    // continue to next buffer if no trigger overlaps
    if (!trigger_overlap) continue;

    // check to make sure we haven't already dumped it
    void* pItem = bsearch (
        &utc_seconds, dumped_times, DUMP_HISTORY, sizeof (time_t), ttcmp);
    if (pItem==NULL) {
      working_bufs[nbuf_to_write] = bufs[ibuf];
      working_times[nbuf_to_write++] = utc_seconds;
    }
  }
  multilog (mlog, LOG_INFO, "earliest buffer t0 = %.3f\n",min_time);
  multilog (mlog, LOG_INFO, "latest   buffer t0 = %.3f\n",max_time);

  // Insert the new dumps into the dump history and re-sort
  for (int ibuf=0; ibuf < nbuf_to_write; ++ibuf) {
    dumped_times[ibuf] = working_times[ibuf];
    multilog (mlog, LOG_INFO, "dumping UTC %d.\n",working_times[ibuf]);
  }
  qsort (dumped_times, DUMP_HISTORY, sizeof (time_t), ttcmp);

  multilog (mlog, LOG_INFO, "Will write %d buffers.\n",nbuf_to_write);

  return nbuf_to_write;
}

/* Examine the "working" buffers, and convert them to dump request structs,
 * populating the reqs queue.
 */
int utcs_to_reqs (int nbuf_to_write) {

  if (nbuf_to_write==0)
    return 0;

  int bufsz = ipcbuf_get_bufsz ((ipcbuf_t*)hdu->data_block);
  int first_empty = compress_reqs ();
  if ((nbuf_to_write + first_empty) > MAX_DUMP_REQUESTS) {
    multilog (mlog, LOG_ERR, "Dump request queue full!  Aborting.\n");
    cleanup (1,1);
  }

  // Write out the earliest times first.
  for (int i0=0; i0 < nbuf_to_write; ++i0) {
    int idx = -1;
    time_t mintime = 3000000000;
    for (int i1=0; i1 < nbuf_to_write; ++i1) {
      if (working_bufs[i1]==NULL)
        continue; // already processed
      if (idx==-1) { // first unprocessed automatically becomes min
        idx = i1;
        mintime = working_times[i1];
        continue;
      }
      if (working_times[i1] < mintime) {
        mintime = working_times[i1];
        idx = i1;
      }
    }
    dumpchk (idx, "IDX must not be <0!\n");
    clear_req (&reqs[first_empty]);
    reqs[first_empty].buf_in = working_bufs[idx];
    reqs[first_empty].bufsz = bufsz;
    reqs[first_empty].utc_epoch = working_times[idx];
    working_times[idx] = 0;
    working_bufs[idx] = NULL;
    first_empty++;
  }
  return nbuf_to_write;
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
    if (mlock (mem_slots[ireq], DUMP_BUFSZ) != 0) {
      multilog (mlog, LOG_ERR, "[DUMPER] Unable to lock memory.\n");
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

  // dump to disk -- TODO -- check out O_DIRECT option
  int sid = getVDIFStationID ((vdif_header*)req->buf_local);
  snprintf (dump_fname, 255,"%s/%s_ea%02d_%li.vdif",
      EVENTDIR,currt_string,sid,req->utc_epoch);
  multilog (mlog, LOG_INFO, "[DUMPER] writing voltages to %s.\n", dump_fname);
  int fd = open (dump_fname, O_WRONLY | O_CREAT, 0664);
  if (fd < 0) {
    multilog (mlog, LOG_ERR, "[DUMPER] Could not open file %s for writing.\n",dump_fname);
    cleanup (5,1);
  }

  uint64_t written = 0;
  while (written < req->bufsz)
  {
    ssize_t result = write (fd, req->buf_local+written, req->bufsz-written);
    if (result < 0)
    {
      int my_errno = errno;
      if (EINTR != my_errno)
      {
        multilog (mlog, LOG_ERR, "[DUMPER] Could not complete write to file %s, encountering errno %d=%s.\n",dump_fname,my_errno, strerror(my_errno));
        cleanup (5,1);
      }
      // TODO -- need to sleep if interrupted?
    }
    else
      written += result;
  }
  if (fsync (fd) != 0) {
    if (ENOSPC == errno)
      multilog (mlog, LOG_ERR, "[DUMPER] SSD is full!\n");
    else
      multilog (mlog, LOG_ERR, "[DUMPER] fsync failed encountering errno %d.\n",errno);
    cleanup (5,1);
  }
  // needed since we don't reset umask, which is 0027
  fchmod (fd, 0664);
  if (close (fd) != 0) {
    multilog (mlog, LOG_ERR, "[DUMPER] close failed encountering errno %d.\n",errno);
    cleanup (5,1);
  }

  // hardcode for now vliteops=6362, vlite=5636
  chown (dump_fname, 6362, 5636);
}

void open_log (char* hostname, int stderr_output) {
  // logging
  gethostname (hostname,MAXHOSTNAME);
  mlog = multilog_open ("dumper",0);
  time_t currt = time (NULL);
  struct tm tmpt;
  localtime_r (&currt,&tmpt);
  strftime (currt_string,sizeof(currt_string), "%Y%m%d_%H%M%S", &tmpt);
  currt_string[15] = '\0';
  char logfile[128];
  pid_t pid = getpid ();
  snprintf (logfile,128,
      "%s/%s_%s_dumper_%06d.log",LOGDIR,currt_string,hostname,pid);
  logfile_fp = fopen (logfile, "w");
  multilog_add (mlog, logfile_fp);
  if (stderr_output)
    multilog_add (mlog, stderr);
  //printf("writing log to %s\n",logfile);
}

void dumpchk (int cond, const char* msg) {
  if (cond < 0) {
    multilog (mlog, LOG_ERR, msg);
    cleanup (-1,1);
  }
}

void sigint_handler (int dummy) {
  cleanup (0, 0);
}



int main(int argc, char** argv) {

  int arg, stderr_output = 1;
  key_t key = 0x40;
  while ((arg = getopt(argc, argv, "hk:o")) != -1) {
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

    case 'o':
      stderr_output = 1;
      break;
    }
  }

  // register SIGINT handling
  signal (SIGINT, sigint_handler);

  // register exit function
  //atexit (exit_handler);

  // TODO -- make this a utility function with "dumper" as arg etc.?
  char hostname[MAXHOSTNAME];
  open_log (hostname, stderr_output);

  hdu = dada_hdu_create (mlog);
  dada_hdu_set_key (hdu,key);
  dumpchk (
      dada_hdu_connect (hdu), "Unable to connect to psrdada HDU.\n");
  multilog (mlog, LOG_INFO, "Connected to psrdada HDU.\n");
  //dumpchk (
      //dada_hdu_lock_read (hdu), "[DUMPER] Unable to lock_read.\n");
  // TODO -- need to make a loop over observations, a la writer/process etc.
  
  // Multi-cast sockets: control and dump requests.
  // For simplicity, just subscribe to the "writer" group.  Why not just
  // use a single control group for all processes?
  int maxsock = 0;
  mc_control_sock = open_mc_socket (mc_vlitegrp, MC_WRITER_PORT,
      "Control Socket [Dumper]", &maxsock, mlog);
  char mc_from[24]; //ip address of multicast sender
  char mc_control_buf[32];
  
  // Listen on a specific group which should *only* have dump requests.
  mc_trigger_sock = open_mc_socket (mc_vlitegrp, MC_TRIGGER_PORT,
      "Trigger Socket [Dumper]", &maxsock, mlog);

  // Dump request queue
  for (int ireq=0; ireq < MAX_DUMP_REQUESTS; ++ireq)
    clear_req (reqs + ireq);

  // Set up local memory for voltage copies.
  make_buffs ();

  fd_set readfds;
  struct timeval select_timeout; //timeout for select(), 100ms
  select_timeout.tv_sec = 0;
  select_timeout.tv_usec = 100000;

  // blocking call to begin observation
  // ASSUMPION: WE don't have to do anything else with this is a viewer; it
  // is now connected and we can view the shared memory buffers.  Therefore
  // we don't even worry about the loop over observations.
  // For that matter, do we even need to open this?  Perhaps connecting
  // is enough...
  //dumpchk (dada_hdu_open_view (hdu), "[DUMPER] Failed to open HDU.\n");

  while (1) { // begin loop over select / DADA reads
    FD_ZERO (&readfds);
    FD_SET (mc_control_sock, &readfds);
    FD_SET (mc_trigger_sock, &readfds);
    dumpchk (
      select (maxsock+1,&readfds,NULL,NULL,&select_timeout),
      "[DUMPER] Error calling select.\n");

    // TMP
    fflush (logfile_fp);
    fflush (stderr);
    
    // Only command we are accepting now is quit!
    if (FD_ISSET (mc_control_sock, &readfds)) {
      char cmd = check_for_cmd (
          mc_control_sock, mc_control_buf, 32, logfile_fp);
      if (cmd == CMD_QUIT) {
        multilog (mlog, LOG_INFO, "[DUMPER] rxd command QUIT, exiting.\n");
        cleanup (0,0);
      }
    }

    // Check for triggers.
    if (FD_ISSET (mc_trigger_sock, &readfds)) {
      // TODO -- handle partial rx?
      int ntrig = 0;
      for (; ntrig < MAX_TRIGGERS; ++ntrig) {
        int nbytes = MultiCastReceive (mc_trigger_sock, 
            (char*)(&trigs[ntrig]), sizeof(trigger_t), mc_from);
        if (nbytes == -1)  {
          if ((EWOULDBLOCK == errno) || (EAGAIN == errno))
            break;
        }
        if (nbytes == 0)
          break;
        if (nbytes != sizeof(trigger_t)) {
          multilog (mlog, LOG_ERR,
              "Received a corruped dump request.  Aborting.\n");
          cleanup (-2,1);
        }
        trigger_t* trig = trigs + ntrig;
        multilog (mlog, LOG_INFO,
            "[DUMPER] rxd trigger: t0=%.3f, t1=%.3f, meta=%s.\n",
            trig->t0,trig->t1,trig->meta);
      }
      multilog (mlog, LOG_INFO, "Received %d triggers.\n", ntrig);

      // Now, convert triggers to dump requests.
      int nbuf_to_dump = trigs_to_utcs (ntrig);
      utcs_to_reqs (nbuf_to_dump);
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
        clear_req (reqs+ireq);
        break;
      }
    }
  } // end loop over select / DADA reads

  cleanup (0,0);
  return 0; // will never get here
} // end main
