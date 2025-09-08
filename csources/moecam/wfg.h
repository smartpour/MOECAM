#ifndef _WFG_H_
#define _WFG_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef double OBJECTIVE;

typedef struct
{
	OBJECTIVE *objectives;
} POINT;

typedef struct
{
	int nPoints;
	int n;
	POINT *points;
} FRONT;

typedef struct
{
	int nFronts;
	FRONT *fronts;
} FILECONTENTS;

#ifdef __cplusplus
extern "C" { double mainHV(int n, int m, double* objectives, double* nadir);
}
#endif

FILECONTENTS *readFile(int n, int m, double* objectives);

extern void printContents(FILECONTENTS *);

#endif
