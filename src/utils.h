#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdio.h>
#include "Connection.h"
//#include "def.h"
#include "ipcio.h"

#ifdef _cplusplus
extern "C" {
#endif
int serve(int port, Connection* c);
int wait_for_cmd(Connection* c);
void event_to_file(const ipcio_t* db, FILE* evfd);
int conn(Connection *c);
void disconn(Connection *c);
int openRawSocket(const char *device, int *deviceIndex);
void change_file_owner (FILE* fp, char* user);
#ifdef _cplusplus
}
#endif

#endif
