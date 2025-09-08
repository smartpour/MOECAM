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

extern "C" FILECONTENTS *readFile(char[]);
extern "C" void printContents(FILECONTENTS *);




class WFG {

public:
	POINT ref; // the reference point 

	int maxm;
	int m_n;

	int fr ;   // current depth 
	FRONT *fs;    // memory management stuff 


	int safe;     // the number of points that don't need sorting
	int** nextp;
	int** prevp;
	int* firstp;
	int* lastp;
	int* psize;
	int maxdepth;


	double mainentryWFG(int m, int n, double* ref, double* pts);
	double mainentryWFG(int m, int n, FRONT* front);

	double hv(FRONT, int* n);

	void makeDominatedBit(FRONT ps, int p, int n);
	double hv2(FRONT ps, int k);
	double inclhv(POINT p, int n);
	double inclhv2(POINT p, POINT q, int n);
	double inclhv3(POINT p, POINT q, POINT r, int n);
	double inclhv4(POINT p, POINT q, POINT r, POINT s, int n);
	double exclhv(FRONT ps, int p, int n);

	void PrepareWFG(int maxn, int maxm);
	void FreeWFG();
};
#else
extern void printContents(FILECONTENTS *);
extern FILECONTENTS *readFile(char[]);

#endif

#endif
