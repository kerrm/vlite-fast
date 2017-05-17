#ifndef __ALERT_H__
#define __ALERT_H__

#define VLA_LOCATION_SIZE	32
#define VLA_DEVICE_SIZE		24
#define VLA_MONITOR_NAME_SIZE	32
#define MAX_ALERT_LENGTH	1024

#include "options.h"

/* related to parsing of EVLAMessage style alerts */
typedef struct
{
	char location[VLA_LOCATION_SIZE];
	char *locationSuffix;	/* points into location at the point after first "-" if cannonical */
	int vlaAnt;		/* (1 through 28) */
	int subarrayAntId;	/* filled in by openFlagFiles */
	char deviceName[VLA_DEVICE_SIZE];
	char monitorName[VLA_MONITOR_NAME_SIZE];
	double timeStamp;
	int alertState;		/* 1 = alert active; 0 = alert inactive */
} AlertDocument;

int expandEntityReferences(char *dest, const char *src, int maxLength);

int alertReceiveOpen(const CommandLineOptions *opts);

int alertReceiveClose(int sock);

int parseAlertDocument(AlertDocument *A, const char *message);

void printAlertDocument(const AlertDocument *A);

void fprintAlertDocument(const AlertDocument *A, FILE *fd);

#endif
