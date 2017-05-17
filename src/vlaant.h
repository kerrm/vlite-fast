#ifndef __VLAANT_H__
#define __VLAANT_H__

#include <stdio.h>
#include "defaults.h"

typedef struct
{
	int number;					/* this is kept set equal to array index number */
	char datasetId[EXECUTOR_DATASETID_SIZE];	/* set by antenna properties document -- null for not set */
	double X, Y, Z;					/* [m] antenna positions */
	double axisOffset;				/* [m] */
} VLAAntenna;

int vlaSubarraySize(const VLAAntenna *vlaAntennas, const char *datasetId);

void fprintVLAAntenna(FILE *out, const VLAAntenna *vlaAntenna, int antennaNumber, int indent);

void printVLAAntenna(const VLAAntenna *vlaAntenna, int antennaNumber, int indent);

#endif
