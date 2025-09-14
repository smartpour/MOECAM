#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "wfg.h"

#define BEATS(x,y)   (x > y)
#define WORSE(x,y)   (BEATS(y,x) ? (x) : (y))

int n;     // the number of objectives
POINT ref; // the reference point

FRONT *fs;    // memory management stuff
int fr = 0;   // current depth
int maxm = 0; // identify the biggest fronts in the file
int maxn = 0;
int safe;     // the number of points that don't need sorting
int** nextp;
int** prevp;
int* firstp;
int* lastp;
int* psize;

double totaltime;

double hv(FRONT);

int greater(const void *v1, const void *v2)
// this sorts points worsening in the last objective
{
	POINT p = *(POINT*)v1;
	POINT q = *(POINT*)v2;
	for (int i = n - 1; i >= 0; i--) {
		if BEATS(p.objectives[i],q.objectives[i]) {
			return -1;
		}
		else if BEATS(q.objectives[i],p.objectives[i]) {
			return  1;
		}
	}
	return 0;
}

int greaterabbrev(const void *v1, const void *v2)
// this sorts points worsening in the penultimate objective
{
	POINT p = *(POINT*)v1;
	POINT q = *(POINT*)v2;
	for (int i = n - 2; i >= 0; i--) {
		if BEATS(p.objectives[i],q.objectives[i]) {
			return -1;
		}
		else if BEATS(q.objectives[i],p.objectives[i]) {
			return  1;
		}
	}
	return 0;
}

int dominates2way(POINT p, POINT q, int k)
// checks if p dominates q or vice versa, in the first k objectives
// returns -1 if p dominates q, 1 if q dominates p, 2 if incomparable
{
	int i;
	int p_worse_count = 0;
	int q_worse_count = 0;

	for (i = k - 1; i >= 0; i--) {
		if BEATS(p.objectives[i],q.objectives[i]) {
			q_worse_count++;
		}
		else if BEATS(q.objectives[i],p.objectives[i]) {
			p_worse_count++;
		}
		if (p_worse_count > 0 && q_worse_count > 0) {
			return 2;
		}
	}
	if (q_worse_count > 0) {
		return -1;
	}
	else if (p_worse_count > 0) {
		return 1;
	}
	else {
		return 2;
	}
}

int dominates(POINT p, POINT q, int k)
// checks if p dominates q or vice versa, in the first k objectives
// returns -1 if p dominates q, 1 if q dominates p, 0 if incomparable
{
	int i;
	int p_worse_count = 0;
	int q_worse_count = 0;

	for (i = k - 1; i >= 0; i--) {
		if BEATS(p.objectives[i],q.objectives[i]) {
			q_worse_count++;
		}
		else if BEATS(q.objectives[i],p.objectives[i]) {
			p_worse_count++;
		}
		if (p_worse_count > 0 && q_worse_count > 0) {
			return 0;
		}
	}
	if (q_worse_count > 0) {
		return -1;
	}
	else if (p_worse_count > 0) {
		return 1;
	}
	else {
		return 0;
	}
}

void makeDominatedBit(FRONT ps, int p)
// creates the front ps[0 .. p-1] in fs[fr], with each point bounded by ps[p] and dominated points removed
{
	int i, j;
	int l = 0;
	for (i = p - 1; i >= 0; i--) {
		if (dominates(ps.points[i], ps.points[p], n - 1) >= 0) {
			fs[fr].points[l] = ps.points[i];
			l++;
		}
	}
	fs[fr].nPoints = l;
}

double hv2(FRONT ps, int k)
// returns the hypervolume of ps[0 .. ps.nPoints - 1] in 2D
// assumes ps is sorted improving
{
	int i;
	double volume = fabs((ps.points[0].objectives[0] - ref.objectives[0]) * (ps.points[0].objectives[1] - ref.objectives[1]));
	for (i = 1; i < ps.nPoints; i++) {
		if (ps.points[i].objectives[1] > ps.points[i - 1].objectives[1]) {
			volume += fabs((ps.points[i].objectives[0] - ref.objectives[0]) * (ps.points[i].objectives[1] - ps.points[i - 1].objectives[1]));
		}
	}
	return volume;
}

double inclusive_hv(FRONT ps)
// returns the inclusive hypervolume of ps[0 .. ps.nPoints - 1]
{
	double volume = 0.0;
	int i;

	if (ps.nPoints == 0) {
		return 0.0;
	}

	if (n == 0) {
		return 0.0;
	}
	else if (n == 1) {
		return fabs(ps.points[0].objectives[0] - ref.objectives[0]);
	}
	else if (n == 2) {
		for (i = ps.nPoints - 1; i >= 0; i--) {
			if (i == 0) {
				volume += fabs((ps.points[i].objectives[0] - ref.objectives[0]) * (ps.points[i].objectives[1] - ref.objectives[1]));
			}
			else {
				volume += fabs((ps.points[i].objectives[0] - ref.objectives[0]) * (ps.points[i].objectives[1] - ps.points[i - 1].objectives[1]));
			}
		}
	}
	else {
		n--;
		for (i = ps.nPoints - 1; i >= 0; i--) {
			POINT tmpref = ref;
			int j;
			for (j = 0; j < n; j++) {
				ref.objectives[j] = WORSE(ref.objectives[j], ps.points[i].objectives[j]);
			}
			if (i == 0) {
				fs[fr].nPoints = 0;
				fr++;
			}
			else {
				makeDominatedBit(ps, i);
				fr++;
			}
			if (fs[fr - 1].nPoints < safe) {
				volume += fabs(ps.points[i].objectives[n] - tmpref.objectives[n]) * hv(fs[fr - 1]);
			}
			else {
				volume += fabs(ps.points[i].objectives[n] - tmpref.objectives[n]) * inclusive_hv(fs[fr - 1]);
			}
			ref = tmpref;
			fr--;
		}
		n++;
	}
	return volume;
}

double hv(FRONT ps)
// returns the hypervolume of ps[0 .. ps.nPoints - 1]
{
	// process small fronts with the IEA
	if (ps.nPoints < safe || n <= 2) {
		return inclusive_hv(ps);
	}
	// these points need IWFG
	else {
		qsort(&ps.points[safe], ps.nPoints - safe, sizeof(POINT), greater);
		return inclusive_hv(ps);
	}
}
