#include <stdio.h>
#include <string.h>
#include "vlaant.h"


/* note plurality of antenna arguments.  In this file, "vlaAntennas" refers to an array of length VLA_ANTENNA_COUNT+1 but "vlaAntenna" refers to a pointer to a single antenna structure */


/* returns number of antennas in the provided subarray (datasetId) */
int vlaSubarraySize(const VLAAntenna *vlaAntennas, const char *datasetId)
{
	int a;
	int size = 0;

	for(a = 1; a <= VLA_ANTENNA_COUNT; ++a)
	{
		if(strcmp(datasetId, vlaAntennas[a].datasetId) == 0)
		{
			++size;
		}
	}

	return size;
}

void fprintVLAAntenna(FILE *out, const VLAAntenna *vlaAntenna, int antennaNumber, int indent)
{
	int i;

	for(i = 0; i < indent; ++i)
	{
		fprintf(out, " ");
	}
	fprintf(out, "VLA antenna %02d  X = %5.3f Y = %5.3f Z = %5.3f axisOffset = %5.3f\n", vlaAntenna->number, vlaAntenna->X, vlaAntenna->Y, vlaAntenna->Z, vlaAntenna->axisOffset);
	if(vlaAntenna->datasetId[0])
	{
		for(i = 0; i < indent+2; ++i)
		{
			fprintf(out, " ");
		}
		fprintf(out, "datasetId = %s\n", vlaAntenna->datasetId);
	}
}

void printVLAAntenna(const VLAAntenna *vlaAntenna, int antennaNumber, int indent)
{
	fprintVLAAntenna(stdout, vlaAntenna, antennaNumber, indent);
}
