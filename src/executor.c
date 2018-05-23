#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <difxmessage.h>
#include "vlite_xml.h"
#include "executor.h"


/* these values correspond to enum ScanInfoDocumentType */
const char ScanInfoTypeString[][24] =
{
	"Unknown",
	"Observation",
	"AntennaProperties",
	"Subarray",
	"#"		/* list terminator; should never be printed */
};

/* these values correspond to enum ArrayConfiguration */
const char ArrayConfigurationString[][8] =
{
	"Unknown",
	"A",
	"AnB",
	"B",
	"BnC",
	"C",
	"CnD",
	"D",
	"#"		/* list terminator; should never be printed */
};

/* these values correspond to enum SubarrayAction */
const char SubarrayActionString[][16] =
{
	"noop",		/* an unnatural state */
	"create",
	"#"		/* list terminator; should never be printed */
};

enum ArrayConfiguration parseArrayConfiguration(const char *arrayName)
{
	enum ArrayConfiguration i;

	for(i = CONFIGURATION_UNKNOWN; i < NUM_CONFIGURATIONS; ++i)
	{
		if(strcmp(arrayName, ArrayConfigurationString[i]) == 0)
		{
			return i;
		}
	}

	return CONFIGURATION_UNKNOWN;
}

enum SubarrayAction parseSubarrayAction(const char *action)
{
	enum SubarrayAction i;

	for(i = SUBARRAY_NOOP; i < NUM_SUBARRAY_ACTIONS; ++i)
	{
		if(strcmp(action, SubarrayActionString[i]) == 0)
		{
			return i;
		}
	}

	return SUBARRAY_NOOP;
}

int subarrayReceiveOpen(const CommandLineOptions *opts)
{
	return openMultiCastSocket(opts->subarrayGroup, opts->subarrayPort);
}

int subarrayReceiveClose(int sock)
{
	return closeMultiCastSocket(sock);
}

/* takes subarray document, identifies and removes any subarrays with conflict, then assigns antennas to new subarray as appropriate */
/* returns number of antennas in vlaAntennas[] array that were removed from an existing subarray */
int mergeSubarray(VLAAntenna *vlaAntennas, const AntPropDocument *ap)
{
	int removed = 0;
	int i;

	/* first remove conflicting subarrays, if any */
	for(i = 1; i <= VLA_ANTENNA_COUNT; ++i)
	{
		if(vlaAntennas[i].datasetId[0] && ap->antenna[i].datasetId[0])
		{
			int j;

			for(j = 0; j < VLA_ANTENNA_COUNT; ++j)
			{
				if(i != j && strcmp(vlaAntennas[i].datasetId, vlaAntennas[j].datasetId) == 0)
				{
					vlaAntennas[j].datasetId[0] = 0;
					++removed;
				}
			}
			vlaAntennas[i].datasetId[0] = 0;
			++removed;
		}
	}

	/* then assign new subarray */
	for(i = 1; i <= VLA_ANTENNA_COUNT; ++i)
	{
		if(ap->antenna[i].datasetId[0])
		{
			strcpy(vlaAntennas[i].datasetId, ap->antenna[i].datasetId);
		}
	}

	return removed;
}

/* This routine adjusts for the VLA center */
void updateAntennaPositions(VLAAntenna *vlaAntennas, const AntPropDocument *ap)
{
	int i;

	for(i = 1; i <= VLA_ANTENNA_COUNT; ++i)
	{
		if(ap->antenna[i].datasetId[0])
		{
			vlaAntennas[i].X = ap->antenna[i].X + VLA_CENTER_X;
			vlaAntennas[i].Y = ap->antenna[i].Y + VLA_CENTER_Y;
			vlaAntennas[i].Z = ap->antenna[i].Z + VLA_CENTER_Z;
			vlaAntennas[i].axisOffset = ap->antenna[i].axisOffset;
		}
	}
}

int scanInfoReceiveOpen(const CommandLineOptions *opts)
{
	return openMultiCastSocket(opts->scanInfoGroup, opts->scanInfoPort);
}

int scanInfoReceiveClose(int sock)
{
	return closeMultiCastSocket(sock);
}

static int parseAntennaId(const char *str)
{
	if(str[0] == 'e' && str[1] == 'a')
	{
		return atoi(str+2);
	}
	else
	{
		return 0;
	}
}

static void XMLCALL startElement(void *userData, const char *name, const char **atts)
{
	ScanInfoDocument *D;
	VliteXML *X;

	X = (VliteXML *)userData;
	D = (ScanInfoDocument *)X->userData;

	X->string[0] = 0;

	if(D->type == SCANINFO_UNKNOWN && X->level == -1)
	{
		if(strcmp(name, "Observation") == 0)
		{
			int i;

			D->type = SCANINFO_OBSERVATION;

			for(i = 0; atts[i]; i += 2)
			{
				if(strcmp(atts[i], "datasetId") == 0)
				{
					snprintf(D->data.observation.datasetId, EXECUTOR_DATASETID_SIZE, "%s", atts[i+1]);
				}
				else if(strcmp(atts[i], "configId") == 0)
				{
					snprintf(D->data.observation.configId, EXECUTOR_CONFIGID_SIZE, "%s", atts[i+1]);
				}
				else if(strcmp(atts[i], "startTime") == 0)
				{
					D->data.observation.startTime = atof(atts[i+1]);
				}
			}
		}
		else if(strcmp(name, "AntennaPropertyTable") == 0)
		{
			int i;

			D->type = SCANINFO_ANTPROP;
			X->antId = 0;
			X->eopId = -1;

			for(i = 0; atts[i]; i += 2)
			{
				if(strcmp(atts[i], "creation") == 0)
				{
					D->data.antProp.creationTime = atof(atts[i+1]);
				}
				else if(strcmp(atts[i], "datasetID") == 0)
				{
					snprintf(D->data.antProp.datasetId, EXECUTOR_DATASETID_SIZE, "%s", atts[i+1]);
				}
				else if(strcmp(atts[i], "configuration") == 0)
				{
					D->data.antProp.arrayConfiguration = parseArrayConfiguration(atts[i+1]);
				}
			}
		}
#if 0
		else if(strcmp(name, "ns2:subArray") == 0)
		{
			int i;

			D->type = SCANINFO_SUBARRAY;

			for(i = 0; atts[i]; i += 2)
			{
				if(strcmp(atts[i], "action") == 0)
				{
					D->data.subarray.action = parseSubarrayAction(atts[i+1]);
				}
				else if(strcmp(atts[i], "activationId") == 0)
				{
					snprintf(D->data.subarray.activationId, EXECUTOR_ACTIVATIONID_SIZE, "%s", atts[i+1]);
				}
				else if(strcmp(atts[i], "configId") == 0)
				{
					snprintf(D->data.subarray.configId, EXECUTOR_CONFIGID_SIZE, "%s", atts[i+1]);
					
				}
				else if(strcmp(atts[i], "name") == 0)
				{
					snprintf(D->data.subarray.name, SUBARRAY_NAME_SIZE, "%s", atts[i+1]);
				}
				else if(strcmp(atts[i], "timeStamp") == 0)
				{
					time_t t;
					char flag;
					char *str;

					str = strdup(atts[i+1]);
					parseISO8601(str, &t, &flag);
					free(str);

					D->data.subarray.timeStamp = MJD_UNIX0 + t/86400.0;
				}
			}
		}
#endif
	}
	else if(D->type == SCANINFO_OBSERVATION)
	{
		if(strcmp(name, "sslo") == 0)
		{
			int i;

			for(i = 0; atts[i]; i += 2)
			{
				if(strcmp(atts[i], "Receiver") == 0)
				{
					if(strcmp(atts[i+1], "300MHz") == 0)
					{
						D->data.observation.usesPband = 1;
					}
					snprintf(D->data.observation.primaryBand, PRIMARY_BAND_SIZE, "%s", atts[i+1]);
				}
			}
		}
	}
	else if(D->type == SCANINFO_ANTPROP)
	{
		if(strcmp(name, "AntennaProperties") == 0)
		{
			int i;

			for(i = 0; atts[i]; i += 2)
			{
				if(strcmp(atts[i], "name") == 0)
				{
					X->antId = parseAntennaId(atts[i+1]);
					if(X->antId > 0 && X->antId <= VLA_ANTENNA_COUNT)
					{
						D->data.antProp.antenna[X->antId].number = X->antId;
						strcpy(D->data.antProp.antenna[X->antId].datasetId, D->data.antProp.datasetId);
					}
				}
			}
		}
		else if(strcmp(name, "eopday") == 0)
		{
			++X->eopId;
		}
	}
	else if(D->type == SCANINFO_SUBARRAY)
	{
		if(strcmp(name, "ns2:station") == 0)
		{
			int i;

			for(i = 0; atts[i]; i += 2)
			{
				if(strcmp(atts[i], "name") == 0)
				{
					int vlaAnt;

					vlaAnt = parseAntennaId(atts[i+1]);
					if(vlaAnt > 0 && vlaAnt <= VLA_ANTENNA_COUNT)
					{
						D->data.subarray.antennaMask[vlaAnt] = 1;
					}
				}
			}
		}
	}

	++X->level;
	strncpy(X->element[X->level], name, VLITE_XML_ELEMENT_SIZE-1);
	X->element[X->level][VLITE_XML_ELEMENT_SIZE-1] = 0;
}

static void XMLCALL endElement(void *userData, const char *name)
{
	ScanInfoDocument *D;
	VliteXML *X;
	const char *elem;

	X = (VliteXML *)userData;
	D = (ScanInfoDocument *)X->userData;

	elem = X->element[X->level];

	if(X->string[0] != 0)
	{
		const char *str = X->string;

		switch(D->type)
		{
		case SCANINFO_OBSERVATION:
			if(strcmp(elem, "name") == 0)
			{
				snprintf(D->data.observation.name, OBSERVATION_NAME_SIZE, "%s", str);
			}
			else if(strcmp(elem, "ra") == 0)
			{
				D->data.observation.ra = atof(str);
			}
			else if(strcmp(elem, "dec") == 0)
			{
				D->data.observation.dec = atof(str);
			}
			else if(strcmp(elem, "dra") == 0)
			{
				D->data.observation.dra = atof(str);
			}
			else if(strcmp(elem, "ddec") == 0)
			{
				D->data.observation.ddec = atof(str);
			}
			else if(strcmp(elem, "azoffs") == 0)
			{
				D->data.observation.azoffs = atof(str);
			}
			else if(strcmp(elem, "eloffs") == 0)
			{
				D->data.observation.eloffs = atof(str);
			}
			else if(strcmp(elem, "startLST") == 0)
			{
				D->data.observation.startLST = atof(str);
			}
			else if(strcmp(elem, "scanNo") == 0)
			{
				D->data.observation.scanNo = atoi(str);
			}
			else if(strcmp(elem, "subscanNo") == 0)
			{
				D->data.observation.subscanNo = atoi(str);
			}
			else if(strcmp(elem, "intent") == 0)
			{
				if(strncmp(str, "ScanIntent=", 11) == 0)
				{
					snprintf(D->data.observation.scanIntent, SCAN_INTENT_SIZE, "%s", str+12);
					D->data.observation.scanIntent[strlen(D->data.observation.scanIntent)-1] = 0;
				}
			}
			break;
		case SCANINFO_ANTPROP:
			if(strcmp(elem, "AntennaProperties") == 0)
			{
				X->antId = 0;
			}
			else if(X->antId > 0 && X->antId <= VLA_ANTENNA_COUNT)
			{
				if(strcmp(elem, "X") == 0)
				{
					D->data.antProp.antenna[X->antId].X = atof(str);
				}
				else if(strcmp(elem, "Y") == 0)
				{
					D->data.antProp.antenna[X->antId].Y = atof(str);
				}
				else if(strcmp(elem, "Z") == 0)
				{
					D->data.antProp.antenna[X->antId].Z = atof(str);
				}
				else if(strcmp(elem, "axisOffset") == 0)
				{
					D->data.antProp.antenna[X->antId].axisOffset = atof(str);
				}
			}
			else if(X->eopId >= 0 && X->eopId < STATE_N_EOP)
			{
				if(strcmp(elem, "epoch") == 0)
				{
					D->data.antProp.eop[X->eopId].mjd = (int)(atof(str)+0.5);
				}
				else if(strcmp(elem, "tai_utc") == 0)
				{
					D->data.antProp.eop[X->eopId].tai_utc = atof(str);
				}
				else if(strcmp(elem, "ut1_utc") == 0)
				{
					D->data.antProp.eop[X->eopId].ut1_utc = atof(str);
				}
				else if(strcmp(elem, "x_pole") == 0)
				{
					D->data.antProp.eop[X->eopId].xPole = atof(str);
				}
				else if(strcmp(elem, "y_pole") == 0)
				{
					D->data.antProp.eop[X->eopId].yPole = atof(str);
				}
			}
			break;
		default:
			break;
		}
	}

	--X->level;
}

int parseScanInfoDocument(ScanInfoDocument *D, const char *message)
{
	memset(D, 0, sizeof(ScanInfoDocument));

	return parseVliteXML(D, startElement, endElement, message);
}

void printScanInfoDocument(const ScanInfoDocument *D)
{
	const ObservationDocument *od;
	const AntPropDocument *ap;
	const SubarrayDocument *sa;
	int a;
	int n;

	printf("ScanInfoDocument\n");
	printf("  type = %d = %s\n", D->type, ScanInfoTypeString[D->type]);
	switch(D->type)
	{
	
	case SCANINFO_OBSERVATION:
		od = &(D->data.observation);
		printf("    datasetId = %s\n", od->datasetId);
		printf("    configId = %s\n", od->configId);
		printf("    startTime = %10.8f\n", od->startTime);
		printf("    name = %s\n", od->name);
		printf("    ra = %10.8f\n", od->ra);
		printf("    dec = %10.8f\n", od->dec);
		printf("    dra = %10.8f\n", od->dra);
		printf("    ddec = %10.8f\n", od->ddec);
		printf("    azoffs = %10.8f\n", od->azoffs);
		printf("    eloffs = %10.8f\n", od->eloffs);
		printf("    startLST = %10.8f\n", od->startLST);
		printf("    scanNo = %d\n", od->scanNo);
		printf("    subscanNo = %d\n", od->subscanNo);
		printf("    primaryBand = %s\n", od->primaryBand);
		printf("    usesPband = %d\n", od->usesPband);
		break;
	
	case SCANINFO_ANTPROP:
		ap = &(D->data.antProp);
		printf("    datasetId = %s\n", ap->datasetId);
		printf("    creationTime = %10.8f\n", ap->creationTime);
		printf("    array configuration = %d = %s\n", ap->arrayConfiguration, ArrayConfigurationString[ap->arrayConfiguration]);
		for(a = 1; a <= VLA_ANTENNA_COUNT; ++a)
		{
			if(ap->antenna[a].number > 0)
			{
				printVLAAntenna(&ap->antenna[a], a, 4);
			}
		}
		for(a = 0; a < STATE_N_EOP; ++a)
		{
      printEOP(&ap->eop[a], 4);
		}
		break;

	case SCANINFO_SUBARRAY:
		sa = &(D->data.subarray);
		printf("    name = %s\n", sa->name);
		printf("    configId = %s\n", sa->configId);
		printf("    activationId = %s\n", sa->activationId);
		printf("    action = %d = %s\n", sa->action, SubarrayActionString[sa->action]);
		printf("    timeStamp = %10.8f\n", sa->timeStamp);
		printf("    members =");
		n = 0;
		for(a = 1; a <= VLA_ANTENNA_COUNT; ++a)
		{
			if(sa->antennaMask[a])
			{
				printf(" %d", a);
				++n;
			}
		}
		printf("    members = %d\n", n);
	
	default:
		break;
	
	}
	printf("\n");
}


void fprintScanInfoDocument(const ScanInfoDocument *D, FILE *fd)
{
	const ObservationDocument *od;
	const AntPropDocument *ap;
	const SubarrayDocument *sa;
	int a;
	int n;

	fprintf(fd,"ScanInfoDocument\n");
	fprintf(fd,"  type = %d = %s\n", D->type, ScanInfoTypeString[D->type]);
	switch(D->type)
	{
	
	case SCANINFO_OBSERVATION:
		od = &(D->data.observation);
		fprintf(fd,"    datasetId = %s\n", od->datasetId);
		fprintf(fd,"    configId = %s\n", od->configId);
		fprintf(fd,"    startTime = %10.8f\n", od->startTime);
		fprintf(fd,"    name = %s\n", od->name);
		fprintf(fd,"    ra = %10.8f\n", od->ra);
		fprintf(fd,"    dec = %10.8f\n", od->dec);
		fprintf(fd,"    dra = %10.8f\n", od->dra);
		fprintf(fd,"    ddec = %10.8f\n", od->ddec);
		fprintf(fd,"    azoffs = %10.8f\n", od->azoffs);
		fprintf(fd,"    eloffs = %10.8f\n", od->eloffs);
		fprintf(fd,"    startLST = %10.8f\n", od->startLST);
		fprintf(fd,"    scanNo = %d\n", od->scanNo);
		fprintf(fd,"    subscanNo = %d\n", od->subscanNo);
		fprintf(fd,"    primaryBand = %s\n", od->primaryBand);
		fprintf(fd,"    usesPband = %d\n", od->usesPband);
		break;
	
	case SCANINFO_ANTPROP:
		ap = &(D->data.antProp);
		fprintf(fd,"    datasetId = %s\n", ap->datasetId);
		fprintf(fd,"    creationTime = %10.8f\n", ap->creationTime);
		fprintf(fd,"    array configuration = %d = %s\n", ap->arrayConfiguration, ArrayConfigurationString[ap->arrayConfiguration]);
		for(a = 1; a <= VLA_ANTENNA_COUNT; ++a)
		{
			if(ap->antenna[a].number > 0)
			{
			  fprintVLAAntenna(fd,&ap->antenna[a], a, 4);
			}
		}
		for(a = 0; a < STATE_N_EOP; ++a)
		{
      printEOP(&ap->eop[a], 4);
		}
		break;

	case SCANINFO_SUBARRAY:
		sa = &(D->data.subarray);
		fprintf(fd,"    name = %s\n", sa->name);
		fprintf(fd,"    configId = %s\n", sa->configId);
		fprintf(fd,"    activationId = %s\n", sa->activationId);
		fprintf(fd,"    action = %d = %s\n", sa->action, SubarrayActionString[sa->action]);
		fprintf(fd,"    timeStamp = %10.8f\n", sa->timeStamp);
		fprintf(fd,"    members =");
		n = 0;
		for(a = 1; a <= VLA_ANTENNA_COUNT; ++a)
		{
		  if(sa->antennaMask[a])
		    {
		      fprintf(fd," %d", a);
		      ++n;
		    }
		}
		fprintf(fd,"    members = %d\n", n);
	
	default:
		break;
	
	}
	fprintf(fd,"\n");
}
