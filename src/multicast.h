#ifndef __MULTICAST_H__
#define __MULTICAST_H__

int MulticastSend(const char *group, int port, const char *message, int length);
int openMultiCastSocket(const char *group, int port);
int closeMultiCastSocket(int sock);
size_t MultiCastReceive(int sock, char *message, int maxlen, char *from);

#endif
