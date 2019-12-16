#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdio.h>
#include "Connection.h"
#include "ipcio.h"
#include "vdifio.h"
#include "multilog.h"
#include "def.h"

#ifdef _cplusplus
extern "C" {
#endif
typedef struct
{
  char hostname[MAXHOSTNAME]; // = "vlite-difx1.evla.nrao.edu";
  char iface[4];
  int gpu_dev_id;
  uint64_t reader_port;
  uint64_t writer_port;
  uint64_t info_port;
  key_t bb_dada_key;
  key_t fb_dada_key;
  int write_fb;
  int nbit;
} VFASTConfig;

typedef struct {
  
  char* buf;          // address to buffer to copy
  uint64_t bufsz;     // size of buffer in bytes
  char fname[256];    // name of output file
  int status;         // status: -1 == working, 0 == success, >0 == error
  int ms_delay;       // milliseconds to delay beginning of execution
  pthread_t tid;      // the thread context

} threadio_t;

typedef struct {
  
  double t0;          // start of trigger window [seconds since Unix epoch]
  double t1;          // end of trigger window [seconds since Unix epoch]
 // float  sn;          // s/n of trigger
 // float  dm;          // dm  of trigger
 // float  width;       // width of trigger
  float  peak_time;   // peaktime of trigger
  char meta[128];     // ancillary information like source name, ...

} trigger_t;

int serve(int port, Connection* c);
int wait_for_cmd(Connection* c, FILE* fp);
int check_for_cmd (int socket, char* buf, int maxlen, FILE* outstream);
int test_for_cmd (int testcmd, int socket, char* buf, int maxlen, FILE* outstream);
void get_cmds (int* cmds, int socket, char* buf, int maxlen, FILE* outstream);
void event_to_file(const ipcio_t* db, FILE* evfd);
int conn(Connection *c);
void disconn(Connection *c);
int openRawSocket(const char *device, int *deviceIndex);
void change_file_owner (FILE* fp, char* user);
VFASTConfig** parse_vfast_config (char* config_file, int* nconfig);
struct timespec get_ms_ts (int ms);
char* print_vfast_config (VFASTConfig* vc, FILE* fp);
time_t vdif_to_unixepoch (vdif_header*);
double vdif_to_dunixepoch (vdif_header* vdhdr, time_t* seconds);
int dump_check_name (char*,char*);
int voltage_check_name (char*,char*);
void* buffer_dump (void* mem);
int open_mc_socket (const char* group, int port, char* name, int* maxnsock, multilog_t* log);

#ifdef _cplusplus
}
#endif

#endif
