#ifndef __EOP_H__
#define __EOP_H__

#include <stdio.h>

/* Earth Orientation Parameters */
typedef struct
{
	int mjd;
	double tai_utc;			/* [s] TAI minus UTC; the leap-second count */
	double ut1_utc;			/* [s] UT1 minus UTC; Earth rotation phase */
	double xPole;			/* [asec] X component of spin axis offset */
	double yPole;			/* [asec] Y component of spin axis offset */
} EOP;

void fprintEOP(FILE *out, const EOP *eop, int indent);

void printEOP(const EOP *eop, int indent);

void copyEOP(EOP *dest, EOP *src, int n);

#endif
