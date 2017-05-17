#ifndef __OPTIONS_H__
#define __OPTIONS_H__

#include "defaults.h"

typedef struct
{
	int verbose;		/* terminal debugging verbosity level */
	int usage;		/* if true, print usage information and quit */
	const char *program;	/* set to argv[0] */
	int alertPort;
	char alertGroup[MULTICAST_GROUP_SIZE];
	int scanInfoPort;	/* for multicast observation, tcal documents */
	char scanInfoGroup[MULTICAST_GROUP_SIZE];	/* '' */
	int subarrayPort;	/* for multicast subarray, antenna properties documents */
	char subarrayGroup[MULTICAST_GROUP_SIZE];	/* '' */
	int allowPband;		/* disable P-band kill switch */
	int doSniff;		/* enable sniffing of data */
	int testMode;		/* 0 = no, 1 = test mode */
} CommandLineOptions;

void resetCommandLineOptions(CommandLineOptions *opts);

CommandLineOptions *newCommandLineOptions(int argc, char **argv);

void deleteCommandLineOptions(CommandLineOptions *opts);

void printCommandLineOptions(const CommandLineOptions *opts);

#endif
