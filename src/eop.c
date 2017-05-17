#include <stdio.h>
#include <string.h>
#include "eop.h"

void fprintEOP(FILE *out, const EOP *eop, int indent)
{
	int i;

	for(i = 0; i < indent; ++i)
	{
		fprintf(out, " ");
	}
	fprintf(out, "EOP %d : TAI-UTC = %8.6f UT1-UTC = %8.6f xPole = %8.6f yPole = %8.6f\n", eop->mjd, eop->tai_utc, eop->ut1_utc, eop->xPole, eop->yPole);
}

void printEOP(const EOP *eop, int indent)
{
	fprintEOP(stdout, eop, indent);
}

void copyEOP(EOP *dest, EOP *src, int n)
{
	int i;

	for(i = 0; i < n; ++i)
	{
		memcpy(dest+i, src+i, sizeof(EOP));
	}
}
