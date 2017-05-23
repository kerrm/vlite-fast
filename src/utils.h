#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdio.h>
#include "Connection.h"
//#include "def.h"
#include "ipcio.h"

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

int serve(int port, Connection* c);
int wait_for_cmd(Connection* c);
void event_to_file(const ipcio_t* db, FILE* evfd);
int conn(Connection *c);
void disconn(Connection *c);
int openRawSocket(const char *device, int *deviceIndex);
void change_file_owner (FILE* fp, char* user);
VFASTConfig** parse_vfast_config (char* config_file, int* nconfig);
char* print_vfast_config (VFASTConfig* vc, FILE* fp);

#ifdef _cplusplus
}
#endif

#endif
