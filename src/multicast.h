#ifndef __MULTICAST_H__
#define __MULTICAST_H__

int MulticastSend(const char *group, int port, const char *message, int length);
int openMultiCastSocket(const char *group, int port);
int closeMultiCastSocket(int sock);
size_t MultiCastReceive(int sock, char *message, int maxlen, char *from);
//
//Multicast group IPs
static const char mc_testgrp[] = "239.199.3.2";
static const char mc_antpropgrp[] = "239.192.3.1";
static const char mc_obsinfogrp[] = "239.192.3.2";
static const char mc_alertgrp[] = "239.192.2.3";
static const char mc_vlitegrp[] = "224.3.29.71";

#define MC_READER_PORT 20000
#define MC_WRITER_PORT 20001
#define MC_INFO_PORT 20002
#define MC_TRIGGER_PORT 20003
#define MC_DUMPER_PORT 20004

//Multicast ports for executor packets
#define MULTI_OBSINFO_PORT 53001
#define MULTI_ANTPROP_PORT 53000
#define MULTI_ALERT_PORT 20011
#define MULTI_TEST_PORT 53901

#endif
