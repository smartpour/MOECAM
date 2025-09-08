
#include "heap.h"

#include <cstdlib>
#include <vector>


typedef void ( *USER_FUNCTION_P)(double *, double *, int *);



#pragma optimize( "", off )
#pragma OPTIMIZE OFF

#ifdef _MSC_VER  
// if the compiler does not recognise this type, change it to another int type 8 bytes long
// like long long int
typedef  __int64 ULINT; //this type myst be 8 bytes long
#else
typedef unsigned long long int ULINT; //this type myst be 8 bytes long
#endif



#define PACKFLOAT(f) ( (f  & 0xFFFFFFFF))
#define f2ulint(a)  (*(unsigned int*) &(a))
#define ulint2f(a) (*(float*) &(a))


KEY_TYPE _getKeyP(node_t theNode)
{	KEY_TYPE ftP; // global
	unsigned int UL1P;
	UL1P=theNode.data>>32;
	ftP=ulint2f(UL1P) ;
	// use the first 4 bytes
	return -ftP;
}

//INDEX_TYPE _getIndex(node_t theNode) {	return 0; }
//void _setIndex(node_t theNode, INDEX_TYPE I) {}




INDEX_TYPE Merge1 (bheap_t *h, float f, unsigned short i, unsigned short j) 
{	
	ULINT ULP;
	unsigned int UI;
	f=-f;
	ULP=f2ulint(f);//<<32;
#ifdef _MSC_VER
       ULP = (ULP << 32) & 0xFFFFFFFF00000000UL ;
#else
       ULP = (ULP << 32) & 0xFFFFFFFF00000000ULL ;
#endif
    UI = i;
	UI = ((UI  << 16) & 0xFFFF0000) | j;
	ULP = ULP + UI;
	return bh_insert(h, ULP); 
}

void GetIndices(DATA_TYPE Node, unsigned short int *i, unsigned short int *j)
{
	//ULINT ULP;
	unsigned int UI;
#ifdef _MSC_VER
       UI = Node & 0x00000000FFFFFFFFUL ;
#else
       UI = Node & 0x00000000FFFFFFFFULL ;
#endif
	   *j = UI & 0x0000FFFF;
	   *i = (UI >> 16) & 0x0000FFFF;
}

#pragma optimize( "", on ) 


void ComputeMin(double x1, double x2, double f1, double f2, double M, double* t, float* f)
{
	*t= 0.5*(x1+x2) + 0.5/M*(f1-f2);
	*f= 0.5*(f1+f2) + 0.5*M*(x1-x2);
}

int Pijavski(double* x0, double *val, USER_FUNCTION_P F, double* Lip, double* Xl, double* Xu, double* precision, int* maxiter)
{
	bheap_t* HeapP = bh_alloc(); 
	DATA_TYPE Node;	
	std::vector<double> x, f;

	
	if (*maxiter>=0xFFFD) *maxiter=0xFFFD; // cannot use short indices!

	int Iter=0;
	double Best=10e10;
	double CurrPrecision=10e10;


	double t,t1, t2,t3,f1,f2, f3, M=*Lip;
	float ff;

	unsigned short int i,j,k;
	int fun = 0;

	t1=*Xl;
	F(&t1,&f1,&fun);     if(Best>f1) { Best=f1; *x0=t1;}
	x.push_back(t1);
	f.push_back(f1);

	t2=*Xu;
	F(&t2,&f2, &fun);		if(Best>f2) { Best=f2; *x0=t2;}
	x.push_back(t2);
	f.push_back(f2);

	i=0; j=1;
	ComputeMin(t1,t2,f1,f2,M,&t,&ff);
	CurrPrecision = Best - ff;

	Merge1(HeapP, -ff, i,j);

	Iter=1;
	while(Iter < *maxiter && CurrPrecision > *precision ) {
		Iter++;
		Node=bh_delete_min(HeapP);
		GetIndices(Node, &i, &j);

		t1=x[i]; t2=x[j]; f1=f[i]; f2=f[j];

		ComputeMin( t1,t2,f1,f2,  M, &t3, &ff);

		if(ff>= *val) // cannot achieve required minimum
		{
			goto exitL;
		}

		F(&t3,&f3, &fun);
		//// GB here do multioblective stuff
		// we ned to have several x, f , etc...


		x.push_back(t3);
		f.push_back(f3);
		k=Iter;

		if(Best>f3) { Best=f3; *x0=t3;}
		CurrPrecision = Best - ff;
		
// two new minima
		ComputeMin(t1,t3,f1,f3,M,&t,&ff);
		Merge1(HeapP, -ff, i,k);

	//	ComputeMin(t3,t2,f3,f2,M,&t,&ff);  // always the same value ff
		Merge1(HeapP, -ff, k,j);
	
	}

exitL:

	*val=Best;

	bh_free(HeapP);
	x.clear(); f.clear();

	*precision=CurrPrecision;
	*maxiter=Iter;
	if(CurrPrecision<0) return -1;  // too small Lip const


	return 0;
}
