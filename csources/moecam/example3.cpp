#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <set>

using namespace std;

#include "gansoc.h"

#include "wfg.h"

#include "direct-master/src/direct.h"




int method = 0;
int func = 0;
int ToPrint = 0;
int NumObj = 1;
double Pi = 3.141592;
double sqr(double a) { return a*a; }


//typedef void (*USER_FUNCTION)(int *, double *, double *);
double CC[10] = { -6.089, -17.164, -34.054, -5.9514, -24.721, -14.986, -24.1, -10.708, -26.662, -22.179};

double recordedmin=10e10;
double recordedmin2 = 10e10;

vector<double> Obj1, Obj2, Obj3;

ofstream out1;
//FILE* out1;
//= fopen("outobj.txt", "w");

void myfunction(int* n, double* x, double* f)
{
	int i;
	*f=0;
//	for(i=0;i<*n;i++)
//		*f += x[i]*x[i];

      double t1=0.0;
	  for(i=0;i<*n;i++)
		  t1 +=x[i];
//t1+=exp(x[i]);


      double objf=0.0;
      for(i=0;i<*n;i++)
          objf +=x[i]*(CC[i]+log(fabs(x[i]/t1)) );
//		objf += exp(x[i])*(CC[i] +x[i] -t1);

	  *f=objf;

}


void myfunction1(int* n, double* x, double* f)
{
	double g;
	int i;

// Griewank
	*f=0;
	g=1;
	for( i=0;i<*n;i++)
	{
		*f += sqr((x[i]-03.2-0))/400.;
		g *= cos((x[i]-03.2-00)/sqrt(i+1.));
	}
	*f +=-g +1;
	*f -=3;
	if(*f<recordedmin) {recordedmin=*f; cout << "1 "<<recordedmin<<endl;}
}

void myfunction2(int* n, double* x, double* f)
{
	double g;
	int i;

	// Griewank
	*f = 0;
	g = 1;
	for (i = 0;i < *n;i++)
	{
		*f += sqr(((1-x[i]) - 03.2 - 0)) / 400.;
		g *= cos(((1-x[i]) - 03.2 - 00) / sqrt(i + 1.));
	}
	*f += -g + 1;
	*f -= 3;
	if (*f < recordedmin2) { recordedmin2 = *f; cout << "2 "<<recordedmin2 << endl; }
}

#define MMmax(a,b) ((a)>(b)? (a):(b))
#define sqr(a) ((a)*(a))

void MAD1(int* n, double* x, double* f)
{
	int i;
	*f=0;
	double f1,f2,f3;

	f1= x[0]*x[0] +x[1]*x[1] + x[0]*x[1] -1;
	f2 =sin(x[0]);
	f3 = - cos(x[1]);

	*f=MMmax(f1,MMmax(f2,f3));
}

void Prob28(int* n, double* x, double* f)
{
	int i;
	*f=0;
	double f1,f2,f3;

	f1 = sqr(x[0]*x[0] + x[1] -11);
	f2 = sqr(x[0] + x[1]*x[1] -7);

	*f=f1+f2;
}

double g(double x, double y){
 return (sqr(x-2*y+1)+sqr(3*x+y))*exp((-x*x-y*y)/2);
}
double g1(double x, double y){
 return (sqr(x-4*y+1)+sqr(1*x+2*y))*exp((-x*x-y*y)/2);
}
double g2(double x, double y){
 return -1/(sqr(1+x*x+y*y));
}


void Test4_1(int* n, double* x, double* f)
{
	*f= x[0] * x[0]-x[1] ;
}

void Test4_2(int* n, double* x, double* f)
{
	*f=  -0.5*x[0]  - x[1]-1;
}


void Kursawe_1(int* n, double* x, double* f)
{
	*f = 0;
	for (int i = 0;i < *n-1;i++) *f += -10 * exp(-0.2*sqrt(x[i] * x[i] + x[i + 1] * x[i + 1]));
}

void Kursawe_2(int* n, double* x, double* f)
{
	*f = 0;
	for (int i = 0;i < *n;i++) *f += pow(fabs(x[i]), 0.8) + 5 * sin(x[i] * x[i] * x[i]);
}


void Fon_1(int* n, double* x, double* f)
{
	*f = 0;
	for (int i = 0;i < *n;i++) *f +=sqr( (x[i]-1./sqrt(*n))) ;
	*f = 1 - exp(-*f);
}

void Fon_2(int* n, double* x, double* f)
{
	*f = 0;
	for (int i = 0;i < *n;i++) *f += sqr((x[i] + 1. / sqrt(*n)));
	*f = 1 - exp(-*f);
}

void ZDT1_1(int* n, double* x, double* f)
{	*f = x[0]; }

double h_ZDT(double a, double b) { return 1 - sqrt(a/b); }
void ZDT1_2(int* n, double* x, double* f)
{
	double g = 0;
	for (int i = 1; i < *n; i++) g += x[i];
	g *= (9. / (*n-1));
	g += 1;

	*f = g*h_ZDT(x[0],g);
}

void ZDT2_1(int* n, double* x, double* f)
{
	*f = x[0];
}

double h2_ZDT(double a, double b) { return 1 - sqr(a / b); }
void ZDT2_2(int* n, double* x, double* f)
{
	double g = 0;
	for (int i = 1; i < *n; i++) g += x[i];
	g *= (9. / (*n - 1));
	g += 1;

	*f = g * h2_ZDT(x[0], g);
}
void ZDT3_1(int* n, double* x, double* f)
{
	*f = x[0];
}

double h3_ZDT(double a, double b) { return 1 - sqrt(a / b) -(a/b)*sin(10*3.1415*a); }
void ZDT3_2(int* n, double* x, double* f)
{
	double g = 0;
	for (int i = 1; i < *n; i++) g += x[i];
	g *= (9. / (*n - 1));
	g += 1;

	*f = g * h3_ZDT(x[0], g);
}

void LTDZ1_1(int* n, double* x, double* f)
{
	*f = 3 - (1 + x[2])*cos(x[0] * Pi / 2)*cos(x[1] * Pi / 2);
	*f = -*f;
}
void LTDZ1_2(int* n, double* x, double* f)
{
	*f = 3 - (1 + x[2])*cos(x[0] * Pi / 2)*sin(x[1] * Pi / 2);
	*f = -*f;
}
void LTDZ1_3(int* n, double* x, double* f)
{
	*f = 3 - (1 + x[2])*sin(x[1] * Pi / 2);
	*f = -*f;
}

void SK2_1(int* n, double* x, double* f)
{
	*f = -(-sqr(x[0] - 2) - sqr(x[1] - 3) - sqr(x[2] - 5) - sqr(x[3] - 4) + 5);
}
void SK2_2(int* n, double* x, double* f)
{
	*f = sin(x[0])+ sin(x[1])+sin(x[2])+sin(x[3]);
	*f /= (1 + 0.01*(sqr(x[0]) + sqr(x[1]) + sqr(x[2]) + sqr(x[3]) )  );
	*f = - *f;
}
void TKLY1_1(int* n, double* x, double* f)
{
	*f = (x[0]+0.1);
}
void TKLY1_2(int* n, double* x, double* f)
{
	double g=1;
	for (int i = 1;i < 4;i++)  g *= (2.0 - exp(-sqr((x[i]-0.1)/0.004)) - 0.8*exp(-sqr((x[i]-0.9)/0.4)));
	*f = 1/(x[0]+0.1)   * g;
}

void VU1_1(int* n, double* x, double* f)
{
	*f = 1.0/(sqr(x[0])+sqr(x[1]) + 1);
}
void VU1_2(int* n, double* x, double* f)
{
	*f = sqr(x[0]) + 3*sqr(x[1]) + 1;
}
void VU2_1(int* n, double* x, double* f)
{
	*f = (x[0]) + (x[1]) + 1;
}
void VU2_2(int* n, double* x, double* f)
{
	*f = sqr(x[0]) + 2 * (x[1]) - 1;
}

vector<double> Pareto;
set<long int> Paretoindex;
list<long int> Erase;

double CurrentVolume;
int countParetoObj;

void ParetoFront1(int* n, int m, double* x, double* f)
{
	int d = 0;
	int nd;
	Erase.clear();
	*f = 0;

	double TotalVolume = CurrentVolume;

	for (auto it = Paretoindex.begin(); it != Paretoindex.end(); it++) {
		nd = 0;
		if (Obj1[m] <= Obj1[*it] )
			nd++;
		if (Obj2[m] <= Obj2[*it] )
			nd++;

		if (nd == 0) {
			d = 1; // dominated by it, we add it to the list of points dominating obj1, to reduce the associated volume
			Erase.push_back(*it); // you cann not have at the same time dominated or dominating
			
		}

		else if (nd == 2 ) // dominates it,  it needs tobe erased for the list.  numobjectives replaces 2 later on
		{
			Erase.push_back(*it); // how to reset the queue?
		}
	} //otherwise nee pareto point which does not dominate any other pareto

	if (!d) {  // d==0, means this point is on Pareto front. recalculate the volume, and remove the dominated items
		double* objectives = new double[2 * Erase.size()];
		double *nadir = new double[2];
		double *ideal = new double[2];
		int k = 0;

		ideal[0] = Obj1[m];
		ideal[1] = Obj2[m];
		nadir[0] = ideal[0];
		nadir[1] = ideal[1];


		for (auto it = Erase.begin(); it != Erase.end(); it++) { // see if we need to replace /excude any paretoindex dominated by i
			Paretoindex.erase(*it);  // remove dominated points if any

			objectives[k * 2]= Obj1[*it];
			objectives[k * 2 + 1] = Obj2[*it];
			k++;

			if (nadir[0] <= Obj1[*it])nadir[0] = Obj1[*it];
			if (nadir[1] <= Obj2[*it])nadir[1] = Obj2[*it];
			// calculate 
			//*f -= fabs(Obj3[*it]*10 );
			//*f -= 1;
		}



		double tempvol= mainHV(2, 1, nadir, ideal);


		TotalVolume -= mainHV(2, Erase.size(), objectives, nadir);  // the volume has increased
		TotalVolume += tempvol;


		// the case of no dominated points, what to do with the volume?
		CurrentVolume = TotalVolume;

		Paretoindex.insert(m);  // ad dthis point to Pareto
		delete[] objectives;
		delete[] nadir;
		delete[] ideal;
	}
	else {
		// dominated, still calculate the missing volume, and reduce the total vomume, but do not touch CurrentVolume, it has not changed
		double* objectives = new double[2 * Erase.size()];
		double *nadir = new double[2];
		int k = 0;
		for (auto it = Erase.begin(); it != Erase.end(); it++) { // see if we need to replace /excude any paretoindex dominated by i
			
			objectives[k * 2] = Obj1[*it];
			objectives[k * 2 + 1] = Obj2[*it];
			k++;
		}
		nadir[0] = Obj1[m];
		nadir[1] = Obj2[m];

		TotalVolume -= mainHV(2, Erase.size(), objectives, nadir);  // subtract the reduced volume
		delete[] objectives;
		delete[] nadir;
	}

	*f = -TotalVolume/10;

	//*f -= Erase.size();
	//*f /= 10.;
//	*f = *f / pow((Obj1.size() + 1), 0.8) * 100;
}


void ParetoFront(int* n, double* x, double* f)
{
	int i, j;
	*f = 0.0;
	for (i = 0;i < Obj1.size() - 1;i++)
		if ((Obj1[Obj1.size() - 1] >= Obj1[i]) || (Obj2[Obj2.size() - 1] >= Obj2[i])) *f += 1;

	*f = *f / pow((Obj1.size() + 1), 0.8)*1;
//	*f =- *f /( pow( (Obj1.size() + 1) ,0.8))*5;

	//*f = 1.0/(*f+1);//  /(Obj1.size()+1);
}

void myMOfunction(int* n, double* x, double* f, int* t)
{
	switch (*t) {
	case 0:  Kursawe_1(n, x, f); break;
	case 1:  Kursawe_2(n, x, f); break;
//	case 0:  Fon_1(n, x, f); break;
//	case 1:  Fon_2(n, x, f); break;

	case 2:  ParetoFront1(n, countParetoObj, x, f); countParetoObj++; break;
	}

	if (*t == 0) Obj1.push_back(*f);
	if (*t == 1) Obj2.push_back(*f);
	if (*t == 2) Obj3.push_back(*f);

	
//	if (*t == 0) fprintf(out1, "\n %f", *f);
//	else if(*t==1) fprintf(out1, " %f", *f);
//	else if (*t == 2) fprintf(out1, ", %f", *f);
}



void ObjFunTKLY1(int* n, double* x, double* f, int* t)
{
	switch (*t) {
	case 0:  TKLY1_1(n, x, f); break;
	case 1:  TKLY1_2(n, x, f); break;
	}
}
void ObjFunVU1(int* n, double* x, double* f, int* t)
{
	switch (*t) {
	case 0:  VU1_1(n, x, f); break;
	case 1:  VU1_2(n, x, f); break;
	}
}
void ObjFunVU2(int* n, double* x, double* f, int* t)
{
	switch (*t) {
	case 0:  VU2_1(n, x, f); break;
	case 1:  VU2_2(n, x, f); break;
	}
}
void ObjFunSK2(int* n, double* x, double* f, int* t)
{
	switch (*t) {
	case 0:  SK2_1(n, x, f); break;
	case 1:  SK2_2(n, x, f); break;
	}
}
void ObjFunLTDZ1(int* n, double* x, double* f, int* t)
{
	switch (*t) {
	case 0:  LTDZ1_1(n, x, f); break;
	case 1:  LTDZ1_2(n, x, f); break;
	case 2:  LTDZ1_3(n, x, f); break;
	}
}

void ObjFunKur(int* n, double* x, double* f, int* t)
{
	switch (*t) {
	case 0:  Kursawe_1(n, x, f); break;
	case 1:  Kursawe_2(n, x, f); break;
	}
}
void ObjFunZDT1(int* n, double* x, double* f, int* t)
{
	switch (*t) {
	case 0:  ZDT1_1(n, x, f); break;
	case 1:  ZDT1_2(n, x, f); break;
	}
}
void ObjFunZDT2(int* n, double* x, double* f, int* t)
{
	switch (*t) {
	case 0:  ZDT2_1(n, x, f); break;
	case 1:  ZDT2_2(n, x, f); break;
	}
}
void ObjFunZDT3(int* n, double* x, double* f, int* t)
{
	switch (*t) {
	case 0:  ZDT3_1(n, x, f); break;
	case 1:  ZDT3_2(n, x, f); break;
	}
}

void myMOfunctionGen(int* n, double* x, double* f, int* t)
{
	if (*t < NumObj) {
		switch (func) {
		case 1:	ObjFunKur(n, x, f, t); break;
		case 2: ObjFunZDT1(n, x, f, t); break;
		case 3: ObjFunZDT2(n, x, f, t); break;
		case 4: ObjFunZDT3(n, x, f, t); break;
		case 5: ObjFunVU1(n, x, f, t); break;
		case 6: ObjFunVU2(n, x, f, t); break;
		case 7: ObjFunSK2(n, x, f, t); break;
		case 8: ObjFunTKLY1(n, x, f, t); break;
		case 9: ObjFunLTDZ1(n, x, f, t); break;
		}

		if (*t == 0) Obj1.push_back(*f);
		if (*t == 1) Obj2.push_back(*f);
		if (*t == 2) Obj3.push_back(*f);
	//	if (*t == 3) Obj4.push_back(*f);

	//	if (*t == 0) out1 << x[0] << " " << x[1] << " ";
		

		if (ToPrint) out1 << *f << " ";
		if (ToPrint && *t== NumObj-1)
			out1 << endl;
		//	fprintf(out1, "%f ", *f);
	}
	else {

		if(method!=2)
		 ParetoFront1(n, countParetoObj, x, f); countParetoObj++;

			//fprintf(out1, "\n");
	}
	//out1 << *n;
}

void testgleb(int* n, double* z, double* f)
{
	int i;
	*f=0;
	double f1,x=z[0],y=z[1];

	f1=sin(cos(0.2*x*x + y))*sin(sin(1*x+0.2*y*y))+0.04*(x*x+y*x)+
		2.9*cos(x/10+y/10)-
		0.9*(g(x,y)+g1(x-4,y+5)) + 20*g2(x-2,y-8) + 20*g2(x+8,y-5);

	*f=f1;
}

double tst_obj(int n,  double *xy, int *undefined_flag, void *unused)
{
	double x, y, f;
	int t = 0;
	for(t=0;t<=NumObj;t++)
		myMOfunctionGen(&n, xy, &f, &t);


//	t = 1;
//	myMOfunction(&n, xy, &f, &t);
//	t = 2;
//	myMOfunction(&n, xy, &f, &t);
	return f;
}
void myMOfunctiontest(int* n, double* x, double* f, int* t)
{
	int m = 0;
	 *f = tst_obj(*n, x, &m, NULL);
}


#include <string>
int main(int argc, char **argv) {


	long int budget = 125;
	int dim = 2;

	string outputString;

	if (argc != 6)
	{
		cout << "Usage: " << argv[0] << "  outputfile method func dim budget\n";
		system("pause");
		exit(0);
	}


	outputString = string(argv[1]);
	method = atoi(argv[2]);
	func= atoi(argv[3]);
	dim = atoi(argv[4]);
	budget = atol(argv[5]);
	//out1	= fopen(outputString, "w");
	out1.open(outputString);
	

	//process functions
	if (func<= 8) NumObj = 2;
	if (func == 9) NumObj = 3;



	Ganso G;
	
//	double aa[3] = {-1,-1,0};
//	double bb[3] = {-0.5,0,0};
	double aa[3] = {3,1,0};
	double bb[3] = {-2.5,0,0};


	double A[30] = {1,2,2,0,0, 1,0,0,0,1,
					0,0,0,1,2, 1,1,0,0,0,
					0,0,1,0,0, 0,1,1,2,1
	};
	double B[10] = {2,1,1,0,0, 0,0,0,0,0};
	double C[10] = {0,0,1,0,0, 0,1,0,0,0};

	double D[10] = {0.2,0,0,0,0, 0,0,0,0,0};

	int index[10]= {1,4,5,6,7,8,9,0,0,0};

	double boundsL[30];
	double boundsU[30];
	double x0[30] ;
	double val;

	int i,j;
	for(i=0;i<30;i++) {x0[i]=-1.0;
	boundsL[i]=0.0001;
	}
	for(i=0;i<30;i++) {
	boundsL[i]=-5.05;
	boundsU[i]= 5;
	}

	if(func>1)
	for (i = 0;i < 30;i++) {
		boundsL[i] = 0;
		boundsU[i] = 1;
	}
	boundsU[0] = 1;



	if(func==5 || func == 6)
	for (i = 0;i < 30;i++) {
		boundsL[i] = -3; boundsU[i] = 3;
	}
	if (func == 7)
		for (i = 0;i < 30;i++) {
			boundsL[i] = -10; boundsU[i] = 10;
		}
	if (func == 7)
		for (i = 0;i < 30;i++) {
			boundsL[i] = -10; boundsU[i] = 10;
		}

	int retcode;

	G.m_objectives = NumObj+1;

	//Paretoindex.insert(0);

	CurrentVolume = 0;
	countParetoObj = 0;

//	retcode = G.MinimizeRandomStart(3, x0, &val, myMOfunction, 0, 0, NULL, NULL, NULL, NULL, boundsL, boundsU, NULL, 500);

	long int maxits = budget;
	int info;
	double minf;
	int force_stop = 0;

	if (method == 1) ToPrint = 1;
	if (method == 0) {
		ToPrint = 0; maxits += dim + 2;
	}
	if (method == 2) ToPrint = 1;

	if(method==0)
		retcode=G.MinimizeECAM(dim,x0,&val,myMOfunctionGen,0,0,70,NULL,NULL,NULL,NULL, boundsL, boundsU, NULL,-maxits);

//	retcode=G.MinimizeRandomStart(3,x0,&val, myMOfunction,0,0,NULL,NULL,NULL,NULL, boundsL, boundsU, NULL,55000);
	else if(method==1)

		info = direct_optimize(tst_obj, NULL, dim, boundsL, boundsU, x0, &minf,
			maxits, 500,
		0, 0, 0, 0,
		0.0, -1.0,
		&force_stop,
		DIRECT_UNKNOWN_FGLOBAL, 0,
		//stdout 
		NULL, DIRECT_ORIGINAL);

	else if(method==2)
		retcode = G.MinimizeRandomStart(dim, x0, &val, myMOfunctionGen, 0, 0, NULL, NULL, NULL, NULL, boundsL, boundsU, NULL, maxits);

//	double r;
//	G.ObjectiveF(D,&r);
//	cout<<val<<endl;
//	cout<<x0[0]<<" "<<x0[1]<<endl;

	out1.close();
	//fclose(out1);
	return 0;
}


/*
/* has two global minima at (0.09,-0.71) and (-0.09,0.71), plus
   4 additional local minima 
static int cnt = 0;
double tst_obj(int n, const double *xy, int *undefined_flag, void *unused)
{
	double x, y, f;
	x = xy[0];
	y = xy[1];
	f = ((x*x)*(4 - 2.1*(x*x) + ((x*x)*(x*x)) / 3) + x * y + (y*y)*(-4 + 4 * (y*y)));
	printf("feval:, %d, %g, %g, %g\n", ++cnt, x, y, f);
	return f;
}

int main(int argc, char **argv)
{
	int n = 2;
	double x[2], l[2], u[2];
	long int maxits = 0;
	int info;
	double minf;
	int force_stop = 0;

	maxits = argc < 2 ? 500 : atoi(argv[1]);

	l[0] = -3; l[1] = -3;
	u[0] = 3; u[1] = 3;

	info = direct_optimize(tst_obj, NULL, n, l, u, x, &minf,
		maxits, 500,
		0, 0, 0, 0,
		0.0, -1.0,
		&force_stop,
		DIRECT_UNKNOWN_FGLOBAL, 0,
		stdout, DIRECT_ORIGINAL);

	printf("min f = %g at (%g,%g) after %d evals, return value %d\n",
		minf, x[0], x[1], cnt, info);

	getchar();
	return EXIT_SUCCESS;
}


*/