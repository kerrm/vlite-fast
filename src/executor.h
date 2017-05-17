#ifndef __EXECUTOR_H__
#define __EXECUTOR_H__

#include "options.h"
#include "eop.h"
#include "defaults.h"
#include "vlaant.h"

#define EXECUTOR_CONFIGID_SIZE		64
#define EXECUTOR_ACTIVATIONID_SIZE	64
#define OBSERVATION_NAME_SIZE		32
#define SUBARRAY_NAME_SIZE		32
#define PRIMARY_BAND_SIZE		8


/* Note! Keep this enum in sync with ScanInfoTypeString[][24] in executor.c */
enum ScanInfoDocumentType
{
	SCANINFO_UNKNOWN	= 0,
	SCANINFO_OBSERVATION,
	SCANINFO_ANTPROP,
	SCANINFO_SUBARRAY,
	NUM_SCANINFO_TYPES		/* this needs to end the list of enumerations */
};

extern const char ScanInfoTypeString[][24];



/* Note! Keep this enum in sync with ArrayConfigurationString[][8] in executor.c */
enum ArrayConfiguration
{
	CONFIGURATION_UNKNOWN	= 0,
	CONFIGURATION_A,
	CONFIGURATION_AnB,
	CONFIGURATION_B,
	CONFIGURATION_BnC,
	CONFIGURATION_C,
	CONFIGURATION_CnD,
	CONFIGURATION_D,
	NUM_CONFIGURATIONS		/* this needs to end the list of enumerations */
};

extern const char ArrayConfigurationString[][8];

enum ArrayConfiguration parseArrayConfiguration(const char *arrayName);




/* Note! Keep this enum in sync with SubarrayActionString[][16] in executor.c */
enum SubarrayAction
{
	SUBARRAY_NOOP	= 0,
	SUBARRAY_CREATE,
	NUM_SUBARRAY_ACTIONS		/* this needs to end the list of enumerations */
};

extern const char SubarrayActionString[][16];

enum SubarrayAction parseSubarrayAction(const char *action);

typedef struct
{
	double startTime;	/* UT MJD */
	char datasetId[EXECUTOR_DATASETID_SIZE];
	char configId[EXECUTOR_CONFIGID_SIZE];
	char name[OBSERVATION_NAME_SIZE];		/* name of source */
	double ra;
	double dec;
	double dra;
	double ddec;
	double azoffs;
	double eloffs;
	double startLST;
	int scanNo;
	int subscanNo;
	char primaryBand[PRIMARY_BAND_SIZE];
	char scanIntent[SCAN_INTENT_SIZE];
	int usesPband;
} ObservationDocument;

typedef struct
{
	double creationTime;	/* UT MJD */
	char datasetId[EXECUTOR_DATASETID_SIZE];
	enum ArrayConfiguration arrayConfiguration;
	VLAAntenna antenna[VLA_ANTENNA_COUNT+1];	/* indexed by VLA antenna number */
	EOP eop[STATE_N_EOP];
} AntPropDocument;

typedef struct
{
	double timeStamp;
	enum SubarrayAction action;
	char activationId[EXECUTOR_ACTIVATIONID_SIZE];
	char configId[EXECUTOR_CONFIGID_SIZE];
	char name[SUBARRAY_NAME_SIZE];
	int antennaMask[VLA_ANTENNA_COUNT+1];	/* indexed by VLA antenna number; 0 -> not in subarray, 1 -> in subarray */
} SubarrayDocument;

typedef struct
{
	enum ScanInfoDocumentType type;
	union
	{
		ObservationDocument observation;
		AntPropDocument antProp;
		SubarrayDocument subarray;
	} data;
} ScanInfoDocument;

int subarrayReceiveOpen(const CommandLineOptions *opts);

int subarrayReceiveClose(int sock);

int mergeSubarray(VLAAntenna *vlaAntennas, const AntPropDocument *ap);

void updateAntennaPositions(VLAAntenna *vlaAntennas, const AntPropDocument *ap);

int scanInfoReceiveOpen(const CommandLineOptions *opts);

int scanInfoReceiveClose(int sock);

int parseScanInfoDocument(ScanInfoDocument *D, const char *message);

void printScanInfoDocument(const ScanInfoDocument *D);

void fprintScanInfoDocument(const ScanInfoDocument *D, FILE *fd);

#endif
