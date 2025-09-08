
#include "ganso.h"


//#if _MSC_VER > 1000  // windows library conflicts with TNT
#undef min
#undef max
//#endif

// uncomment this  for demo version
//#define DEMO_VERSION

#include "src/qrandom.h"


#include "tnt/tnt.h"
#include "tnt/jama_lu.h"


typedef JAMA::LU<double> LUP;
typedef TNT::Array1D<double> A1D;
typedef TNT::Vector<int> IVec;
typedef TNT::Vector<double> DVec;
typedef TNT::Matrix<double> DMat;

extern int ToPrint;

#ifdef USEECAM
#include "src/ecam.h"
#include "src/pijavski.h"
#endif



#define INdex (*((IVec*)Index))
#define BforECAM (*((DVec*)bforECAM))
#define AforECAM (*((DMat*)aforECAM))
#define A2 (*((DMat*)a2))
#define A1 (*((DMat*)a1))
#define Bvec (*((DVec*)bvec))

// casted to
//	Matrix<double> A2;
//	Matrix<double> A1;
//	Vector<double> B;
//	Vector<int>	Index;
//	Matrix<double> AforECAM;
//	Vector<double> BforECAM;


Ganso::Ganso()
{
	bvec=new DVec;
	a1= new DMat;
	a2= new DMat;
	aforECAM= new DMat;
	bforECAM=new DVec;
	Index=new IVec;
	DSOPenalty=100;

	Basic=NULL; X=NULL;
	m_Xl=NULL; m_Xu=NULL;
	m_BoxConstraints=0;
	Fortran=0;

	ECAMLowerBound=0;
	FunCallsWhenFound=FunctionCounter=0;
	m_objectives = 1;
	MySTC = NULL;

}
	
Ganso::~Ganso()
{
	delete reinterpret_cast<DVec*>(bvec);
	delete reinterpret_cast<DMat*>(a1);
	delete reinterpret_cast<DMat*>(a2);
	delete reinterpret_cast<DMat*>(aforECAM);
	delete reinterpret_cast<DVec*>(bforECAM);
	delete reinterpret_cast<IVec*>(Index);
}
	

// MO: decide which of the objecties is most promising
int Ganso::DecideObjectives()
{
	int i, k=0;
	double crit=-1.0;
	double gap;
	for (i = 0;i < m_objectives;i++) {
		gap=MySTC[i].RequestGap();
		if (gap > crit) { crit = gap; k = i; }
	}

	k = 2;
//	printf("\n%d %f", k,crit);
//	if(MySTC[i].LastLabel>= MySTC[i].Iters/300 && MySTC[i].LastLabel <= 2*MySTC[i].Iters )  return 2;
	return k;
}

void Ganso::MakeStepbefore()
{
	int i;
	for (i = 0;i < m_objectives;i++) 
		MySTC[i].MakeStepbefore();
}

void Ganso::MakeStep(int objchosen)
{
	int i;
	dVec X = MySTC[objchosen].temp_x;

	for (i = 0;i < m_objectives;i++)
		MySTC[i].MakeStep(X);
}



/* The following are the main methods providing interface to the users program.
   
They all have similar structure: 
 1. call SetConstraints to decide on the basic/nonbasic variables and coordinate
	transformation.
 2. SetObjFunction to set the address of user's function
 3. Call one of the methods with m_basic as the dimension and other parameters
 4. Improve the solution by local DFBM method
 5. Convert the result back to original coordinates.

If there are no linear constraints, then no coordinate transformation is performed 
(i.e. the identity transformation is still used) 


Hint: setting maxiter to a NEGATIVE value skips step 4. Useful for checking what
a particular method returned without improvement by a local method
*/

int		Ganso::MinimizeDFBM(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter)
	{
		int ret;
		ret=SetConstraints(lineq,linineq,dim,AE,AI,RHSE,RHSI,Xl,Xu,basic);
		if(ret<0) return ret;

		SetObjFunction(f);
		// get startpoint from x0
		double* start=(double*) malloc(sizeof(double)*m_basic);

		int i;

		GetBasicVariables(x0,start);

		ret = CallDFBM(m_basic,start, val, 0, maxiter);
		// get back to the original coords
		BuildFullVector(start);
		for(i=0;i<dim;i++) x0[i]=X[i];

		free(start);
		return ret;
	}


	int		Ganso::MinimizeRandomStart(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter)
	{
		int ret;
		ret=SetConstraints(lineq,linineq,dim,AE,AI,RHSE,RHSI,Xl,Xu,basic);
		if(ret<0) return ret;

		SetObjFunction(f);
		// get startpoint randomly
		double* start=(double*) malloc(sizeof(double)*m_basic);
		double* best=(double*) malloc(sizeof(double)*m_basic);
		double* Left=(double*) malloc(sizeof(double)*m_basic);
		double* Right=(double*) malloc(sizeof(double)*m_basic);
		int i,j,saveret=0;

		double curval,minval=10e20;

		RandomSearch	RS;

		double* LC = new double[3];
		double* History = new double[dim*maxiter];
		double* distances = new double[ maxiter];
		double* HistoryZ=new double[3*maxiter];  //GB for now fix at 3
		double dist,mindist, temp;
		int distpos;

		RS.Parent=this;

		for(i=0;i<m_basic;i++) {
			Left[i] = GetLowBound(i);
			Right[i]= GetUpBound(i);
		}

		for (i = 0;i < 3;i++) {
			LC[i] = 0;
		}


		for(j=0;j<maxiter;j++) {
			
			// randomly generate starting point
			RS.GenerateOnCube(m_basic,start,Left,Right);

			BuildFullVector(start);
			if(m_box2<=1.0e-4) { // otherwise reject this point

				//ret= CallDFBM(m_basic,start, &curval, 0, 1000);

				for (int kk = 0; kk < dim;kk++) History[dim*j + kk] = start[kk];

				mindist = 10e10; distpos = 0;
				for (int jj = 0; jj < j;jj++) {
					dist = 0;
					for (int kk = 0; kk < dim;kk++)
						dist += fabs(History[dim*jj + kk] - start[kk]);

					distances[jj] = dist;

					if (dist < mindist) {
						mindist = dist;distpos = jj;
					}
				}



				i = 0;
				CurrentFunction(m_basic, &(start[0]), &curval, &i);				
				HistoryZ[3 * j] = curval;

				for (int jj = 0; jj < j;jj++) {
					temp = fabs(curval - HistoryZ[3*jj]) / (distances[jj] + 0.000001);
					if (LC[0] < temp)  LC[0] = temp;
				}
//etimate Lipschtz
				i = 1;
				CurrentFunction(m_basic, &(start[0]), &curval, &i);
				HistoryZ[3 * j+i] = curval;
				for (int jj = 0; jj < j;jj++) {
					temp = fabs(curval - HistoryZ[3*jj+i]) / (distances[jj] + 0.000001);
					if (LC[1] < temp)  LC[1] = temp;
				}


				i = 2;
				CurrentFunction(m_basic, &(start[0]), &curval, &i);
				HistoryZ[3 * j + i] = curval;
				for (int jj = 0; jj < j;jj++) {
					temp = fabs(curval - HistoryZ[3*jj+i]) / (distances[jj] + 0.000001);
					if (LC[2] < temp)  LC[2] = temp;
				}

				if(ret<0) saveret=ret;
				if(minval>curval) { // keep the smallerst value 
						for(i=0;i<m_basic;i++) best[i]=start[i];
						minval=curval;
				}
			} else j--;
		}

	//	ret = CallDFBM(m_basic,best, &curval, 0, 10000); // last iteration higher precision

		// get back to the original coords
		BuildFullVector(best);
		for(i=0;i<dim;i++) x0[i]=X[i];
		*val=curval;

// return in x0


		free(LC);
		free(History);
		free(HistoryZ);

		free(start); free(best);
		free(Left); free(Right);
		if(saveret<0) return -saveret; // warning
		return ret;
	}

int		Ganso::MinimizeDFBMECAM(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter, int iterECAM , int dimECAM)
	{
		int ret;
		ret=	SetConstraints(lineq,linineq,dim,AE,AI,RHSE,RHSI,Xl,Xu,basic);
		SetObjFunction(f);
		if(ret<0) return ret;

		// get startpoint from x0
		double* start=(double*) malloc(sizeof(double)*m_basic);

		GetBasicVariables(x0,start);

		ret = CallDFBM(m_basic,start,val,1,maxiter,iterECAM,dimECAM);
		// get back to the original coords
		BuildFullVector(start);
		for(int i=0;i<dim;i++) x0[i]=X[i];

		free(start);
		return ret;
	}

#ifdef USEECAM
int		Ganso::MinimizeECAM(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double LC, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter, int iterLocal)
	{
		int ret,i,k;
		ret=SetConstraints(lineq,linineq,dim,AE,AI,RHSE,RHSI,Xl,Xu,basic);
		if(ret<0) return ret;

		SetObjFunction(f);

		double* start=(double*) malloc(sizeof(double)*m_basic);
	 	double* Best=(double*) malloc(sizeof(double)*m_basic);

		ret = CallECAM(m_basic,start,val,LC,abs(maxiter),iterLocal,0);


		// improve the solution by DGM
		if(maxiter>=0 && ret == 0)
			ret = CallDFBM(m_basic, start, val, 0, DEFAULT_IterDGM);

		double Best_f=*val;

		if(maxiter>=0 && StartingPointsSize>0) { // try other starting points
			for(k=0;k<StartingPointsSize;k++) {
				 for( i=0;i<m_basic;i++) Best[i]=StartingPoints[k*m_basic+i];
				 CallDFBM(m_basic, Best, val, 0, DEFAULT_IterDGM);

				 if(*val<Best_f) { // save the best optimum
					 Best_f=*val;
					 for( i=0;i<m_basic;i++) start[i]=Best[i];
				 }
			}
			*val=Best_f;
			StartingPointsSize=0;
			free(StartingPoints);
		}

		// get back to the original coords
		BuildFullVector(start);
		for( i=0;i<dim;i++) x0[i]=X[i];

		free(start); free(Best);
		return ret;
	}


int		Ganso::MinimizeECAMDFBM(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double LC, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter)
	{
		int ret;
		ret= SetConstraints(lineq,linineq,dim,AE,AI,RHSE,RHSI,Xl,Xu,basic);
		if(ret<0) return ret;
		SetObjFunction(f);

		double* start=(double*) malloc(sizeof(double)*m_basic);

		ret = CallECAM(m_basic,start,val,LC,abs(maxiter),0,1);

		// improve the solution by DGM
		if(maxiter>=0 && ret==0)
		ret = CallDFBM(m_basic, start, val, 0, DEFAULT_IterDGM);

		// get back to the original coords
		BuildFullVector(start);
		for(int i=0;i<dim;i++) x0[i]=X[i];

		free(start);
		return ret;
	}
#endif




#ifdef USEDSO
	int		Ganso::MinimizeDSO(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq,  double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, double penalty,int speed, int precision)
	{
		int ret;
		ret=SetConstraints(lineq,linineq,dim,AE,AI,RHSE,RHSI,Xl,Xu,basic);
		if(ret<0) return ret;
		SetObjFunction(f);

		DSOPenalty=penalty;
		double* start=(double*) malloc(sizeof(double)*m_basic);
		double* Left=(double*) malloc(sizeof(double)*m_basic);
		double* Right=(double*) malloc(sizeof(double)*m_basic);
		int i;

		for(i=0;i<m_basic;i++) {
			Left[i] = GetLowBound(i);
			Right[i]= GetUpBound(i);
			start[i]=(Left[i]+Right[i])*0.5;  // middle of the box
		}

		ret = CallDSO(m_basic,start, Left, Right, val,abs(speed),precision);

		// improve the solution by DGM
		if(speed>=0 && ret==0)
			ret = CallDFBM(m_basic, start, val, 0, DEFAULT_IterDGM);

		// get back to the original coords
		BuildFullVector(start);
		for( i=0;i<dim;i++) x0[i]=X[i];

		free(start);
		free(Left); free(Right);
		return ret;
	}


	int		Ganso::MinimizeIterativeDSO(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq,  double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, double penalty, int speed, int precision, int rounds)
	{
		int i,ep,ret;
		ret=SetConstraints(lineq,linineq,dim,AE,AI,RHSE,RHSI,Xl,Xu,basic);
		if(ret<0) return ret;

		SetObjFunction(f);

		DSOPenalty=penalty;
		double* start=(double*) malloc(sizeof(double)*m_basic);
		double* Left=(double*) malloc(sizeof(double)*m_basic);
		double* Right=(double*) malloc(sizeof(double)*m_basic);


		for(i=0;i<m_basic;i++) {
			Left[i] = GetLowBound(i);
			Right[i]= GetUpBound(i);
			start[i]=(Left[i]+Right[i])*0.5;  // middle of the box
		}

		if(rounds>0)
			ret = CallDSO(m_basic,start, Left, Right, val,abs(speed),precision, rounds);
		else
			ret = CallDSO(m_basic,start, Left, Right, val,abs(speed),precision, 1);
		double ratio=2;

		if(rounds<0)
		 for(ep=1;ep<abs(rounds);ep++) {
			// new box
			for(i=0;i<m_basic;i++) {
				Left[i] = GetLowBound(i);
				Right[i]= GetUpBound(i);
				Right[i] = (Right[i]-Left[i])/ratio /2;   // half the previous sizesize
				Left[i] =start[i] - Right[i];  // box centered at start
				Right[i]=start[i] + Right[i];
				if(Left[i]< GetLowBound(i)) Left[i]= GetLowBound(i);
				if(Right[i]> GetUpBound(i)) Right[i]= GetUpBound(i);
			}

			ratio *= 2;
			ret = CallDSO(m_basic,start, Left, Right, val,abs(speed),precision);
		 }
		 // if

		// improve the solution by DGM
		if(speed>=0 && ret==0)
			ret = CallDFBM(m_basic, start, val, 0, DEFAULT_IterDGM);

		// get back to the original coords
		BuildFullVector(start);
		for( i=0;i<dim;i++) x0[i]=X[i];

		free(start);
		free(Left); free(Right);
		return ret;
	}

#ifdef USEECAM
	int		Ganso::MinimizeECAMDSO(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double LC, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter, int iterDSO)
	{
		int i,ret;
		ret=SetConstraints(lineq,linineq,dim,AE,AI,RHSE,RHSI,Xl,Xu,basic);
		if(ret<0) return ret;

		SetObjFunction(f);

		LeftT=(double*) malloc(sizeof(double)*m_basic);
	 	RightT=(double*) malloc(sizeof(double)*m_basic);
	 	BestX=(double*) malloc(sizeof(double)*m_basic);
		BestF=10e10;

		double* start=(double*) malloc(sizeof(double)*m_basic);

		ret = CallECAM(m_basic,start,val,LC,abs(maxiter),iterDSO,2); // obj function calls DSO

// retrieve best F
		*val=BestF;
		for( i=0;i<dim;i++) start[i]=BestX[i];

// run DSO as usual in the whole domain?

		// improve the solution by DGM
		if(maxiter>=0 && ret==0)
			ret = CallDFBM(m_basic, start, val, 0, DEFAULT_IterDGM);

		// get back to the original coords
		BuildFullVector(start);
		for( i=0;i<dim;i++) x0[i]=X[i];

		free(start); free(BestX);
		free(LeftT); free(RightT);
		return ret;

	}
#endif
#endif

/*  Simpfied versions for unconstrained problems */
	int		Ganso::MinimizeDFBM_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl, double* Xu, int maxiter)
	{ 	return MinimizeDFBM(dim, x0, val,  f, 0,0,0,0,0,0, Xl, Xu, 0,  maxiter);}

	int		Ganso::MinimizeDFBMECAM_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl, double* Xu,  int maxiter, int iterECAM , int dimECAM)
	{ 	return MinimizeDFBMECAM(dim, x0, val,  f, 0,0,0,0,0,0, Xl, Xu, 0, maxiter,iterECAM,dimECAM);}

	int		Ganso::MinimizeRandomStart_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl, double* Xu, int maxiter)
	{ 	return MinimizeRandomStart(dim, x0, val,  f, 0,0,0,0,0,0, Xl, Xu, 0,  maxiter);}


#ifdef USEECAM
	int		Ganso::MinimizeECAM_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double LC, double* Xl, double* Xu, int maxiter, int iterLocal)
	{ 	return MinimizeECAM(dim, x0, val,  f, 0,0,LC,0,0,0,0, Xl, Xu, 0,  maxiter, iterLocal);}

	int		Ganso::MinimizeECAMDFBM_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		 double LC, double* Xl, double* Xu,  int maxiter)
	{ 	return MinimizeECAMDFBM(dim, x0, val,  f, 0,0,LC,0,0,0,0, Xl, Xu, 0,  maxiter);}
#endif


#ifdef USEDSO
	int		Ganso::MinimizeDSO_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl, double* Xu, int maxiter, int steps)
	{ 	return MinimizeDSO(dim, x0, val,  f, 0,0,0,0,0,0, Xl, Xu, 0,  maxiter,steps);}

	int		Ganso::MinimizeIterativeDSO_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl, double* Xu, int maxiter, int steps, int epoch)
	{ 	return MinimizeIterativeDSO(dim, x0, val,  f, 0,0,0,0,0,0, Xl, Xu, 0,  maxiter,steps,epoch);}

#ifdef USEECAM
	int		Ganso::MinimizeECAMDSO_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double LC, 	double* Xl, double* Xu, int maxiter, int iterDSO)
	{ 	return MinimizeECAMDSO(dim, x0, val,  f, 0,0,LC,0,0,0,0, Xl, Xu, 0,  maxiter,iterDSO);}
#endif
#endif

/*================= these are internal methods ===================================================*/

/* in the CallSomeMethod we set the type of the objective function transformation
   define an instance of the corresponding class and call its minimization routine
*/

int		Ganso::CallDFBM(int dim, double* startpoint, double* val, int UseECAM, 
				int maxiter, int iterECAM, int dimECAM )
{
	int i=0,savedfunctiontype=m_CurrentFunctionType;
	m_CurrentFunctionType=0;
//	CDGM MyDFBM;

	if(dimECAM>dim) dimECAM=dim;

//	MyDFBM.Parent=this;

//	i=MyDFBM.minimise_(dim,startpoint, val, maxiter,  (UseECAM==0)?1:2,  iterECAM, dimECAM);

	m_CurrentFunctionType=savedfunctiontype;

	return i;
}

Ganso * GlobalParent;
#ifdef USEECAM

void AuxUserFunction(double* x, double* v, int* t) {
	GlobalParent->CurrentFunction(1, x,  v,t);
}
int	Ganso::CallECAM (int dim,   double* sol, double* val,  double LC, int iterECAM, int iterLocal,
					 int useLocal)
{
	int i,ret=0;
	MySTC=new SawToothCover[m_objectives];

	for (i = 0;i < m_objectives;i++) {
		MySTC[i].MyParent = this; MySTC[i].m_myobj = i;
	}

	GlobalParent=this;

	if(useLocal==1) m_CurrentFunctionType=3; else  // dfbm
	if(useLocal==2) m_CurrentFunctionType=6; else  // dso
		m_CurrentFunctionType=2;



	double precision=1e-8;
	double xL, xU;
	for(i=0;i<m_basic;i++) sol[i]=m_Xl[ INdex[i]  ]; // shift the origin
	PrepareForECAM(m_basic, sol);  // set up transformation matrices

	if(m_basic==1) {

		// set up bounds	
		xL=0; xU=m_Xu[ INdex[0]  ] - m_Xl[ INdex[0]  ]; 
		*val=1e10;
		ret=Pijavski(sol, val, AuxUserFunction, &LC,  &xL, &xU, &precision, &iterECAM);
		ComputeXBasicShift(m_basic,sol,sol);
		StartingPointsSize=0; //iterLocal - proably this is not needed in 1-d case?
		ECAMLowerBound=*val-precision;
		FunCallsWhenFound=iterECAM;

		delete[] MySTC;
		return ret;
	}

	for(i=0;i<m_basic;i++) sol[i]=m_Xu[INdex[i]]-m_Xl[INdex[i]]; // size of the cube

	for (i = 0;i < m_objectives;i++) {
		MySTC[i].SetBoundaries(m_basic, sol);
		MySTC[i].SetConstants(LC, m_basic);  // Lipschitz constant
		MySTC[i].SetConst(10000);			  // this is an aux parameter - just some big value
		MySTC[i].SetCut(10000);     // cut big values of f(x) (to avoid infinities)
		MySTC[i].SetListSize(iterLocal);
	}

// set these global parameters for DSO
	if(useLocal==2) {StartingPointsSize=iterLocal;Iterations=iterECAM;curriter=0;
		MySTC[0].SetListSize(0);
	}

//	sol[0]=1; sol[1]=0.5; *val=100;
//			printf("entered %d %d %f %g %g\n", iterLocal, m_BoxConstraints,(float) (*val), sol[0], sol[1]);
//	ObjFunction(&m_basic,sol,val);
//	ComputeXBasicShift(m_basic,sol,sol);
//			printf("shift %d %g %g %g\n", m_basic, *val, sol[0], sol[1]);
//		BuildFullVector(sol); 
//		ObjFunction(&m_basic,sol,val);  // users function
//	ObjectiveF(sol,val);
//	printf("exit %d %g %g\n", m_nvar, X[0], X[1]);
//	CurrentFunction(m_basic,sol,val);

	for (i = 0;i < m_objectives;i++)
		MySTC[i].InitialSteps(dim, iterECAM);

	ToPrint = 1;

	MySTC[0].Minimise(m_basic,iterECAM);  // main method
	ret=MySTC[0].LastError();
	StartingPointsSize=MySTC[0].ListSize;
	StartingPoints = (double *)malloc(sizeof(double)*m_basic*StartingPointsSize);
	for(i=0;i<StartingPointsSize*m_basic;i++)
		StartingPoints[i]=MySTC[0].FinalList[i];

	for(i=0;i<m_basic;i++) sol[i] = MySTC[0].xglob_x[i];
	ComputeXBasicShift(m_basic,sol,sol);
	*val=MySTC[0].xglob_val;

	ECAMLowerBound=MySTC[0].GVA;
	FunCallsWhenFound=MySTC[0].FunCallsWhenFound;
	delete[] MySTC;
	return ret;
}
#endif


#ifdef USEDSO

// need a global function wrapper, supplied to DSO
void   FunctionForDSO(long *n, double* x, double* val)
{
	// an instance of Ganso
	GlobalParent->CurrentFunction(*n,x,val);
}

int	Ganso::CallDSO (int dim,   double* sol, double* L,double* R, double* val, int speedDSO, int precisionDSO, int rounds)
{
	m_CurrentFunctionType = 4; // means original with penalty on nonbasic
	GlobalParent=this;

// here comes dll call, use some defaults =================================

	//variables defined
	int  printOut=0;
    
	int  retval = calculate(FunctionForDSO, sol, R, L, rounds, dim, precisionDSO, speedDSO, printOut );
	CurrentFunction(dim,sol,val);

/*	
	//DLL handle
	HINSTANCE hDLL = NULL;
	//DLL function "calculate" handle
	CALCULATE_FUNC calculate;

	//vectors defined
	double results[NUMBER_OF_VARIABLES];
	double xmax[NUMBER_OF_VARIABLES]; 
	double xmin[NUMBER_OF_VARIABLES];

	//return value
	int retval;

	//load the DLL
	hDLL = LoadLibrary("AGOP_v1");
	int i;
	for(i=0;i<dim;i++) {xmax[i]=R[i]; xmin[i]=L[i];}

	//check if load was OK
	if (hDLL == NULL)
	{
		printf("Could not load the DLL.\n");
		return -1;
	}
	else
	{
		//get a function pointer to the "calculate" function from within the DLL
		calculate = (CALCULATE_FUNC)GetProcAddress(hDLL, "calculate");
	}

	if(calculate != NULL) {
//	    retval = calculate(FunctionForDSO, sol, L, R, rounds, numVar, precision, speed, printOut );
	    retval = calculate(FunctionForDSO, results, xmin, xmax, rounds, numVar, precision, speed, printOut );

	} else return -1;

*/
// suppy this function FunctionForDSO

	return retval;
}
#endif


/* Every minimization method calls this method as teh objective function 
   Here we decide which way we transform coordinates, and how to call the users function
*/

int		Ganso::CurrentFunction(int n, double* x, double* val, int* t)
{
	switch(m_CurrentFunctionType) {
	case 0:	 return ObjectiveF(x,val,t);
	case 1:	 return FunctionECAMSubspace(n,x,val);
	case 2:  return FunctionECAM(n,x,val,t);
	case 3:	 return FunctionECAMLocal(n,x,val); 
	//case 4:	 return FunctionDSO(n,x,val); 
	//case 5:	 return FunctionDFBM(n,x,val); 
	//case 6:	 return FunctionECAMDSO(n,x,val); 

	default: return ObjectiveF(x,val,t);
	}
}

// build the full vector and call users function
int		Ganso::ObjectiveF(double * x, double* val, int* t)
{
		BuildFullVector(x); 
		ObjFunction(&m_nvar,X,val,t);  // users function

		FunctionCounter++;

		return 0;
}

// instead of f(x) calls DFBM
int		Ganso::FunctionDFBM(int n, double * x, double* val)
{
	int savedfunctiontype=m_CurrentFunctionType;

	m_CurrentFunctionType=0; // the actual users function
	CallDFBM(n, x,val);

	m_CurrentFunctionType=savedfunctiontype;
	return 0; //??
}
// aux call from ECAM
int		Ganso::FunctionECAMLocal(int n, double * x, double* val)
{	// obtain original variables
	ComputeXBasicShift(n,x,Basic);
	return FunctionDFBM(n,Basic,val);
}

// calls obj function mapping a lower dinemsional subspace to the basic variables
int		Ganso::FunctionECAMSubspace(int n, double * x, double* val)
{
	// obtain original variables
	int i = 0;
	ComputeXBasic(n,x,Basic);
	return ObjectiveF(Basic,val,&i);
}

// shift the origin
int		Ganso::FunctionECAM(int n, double * x, double* val, int* t)
{
	// obtain original variables
	ComputeXBasicShift(n,x,Basic);
	ObjectiveF(Basic,val,t);		
	*val += m_box2 * m_box2 * DSOPenalty;

	return 0;
}

int		Ganso::FunctionECAMDSO(int n, double * x, double* val)
{	// obtain original variables
	ComputeXBasicShift(n,x,Basic);
	int i=0; double curr;

	curriter++;
	if(curriter+StartingPointsSize<Iterations) 
		// too early to start DSO, just obj function
	{
		 ObjectiveF(Basic,val, &i);
		 *val += m_box2 * DSOPenalty;

		if(*val<BestF)
		{
			BestF = *val;
			for(i=0;i<m_basic;i++) BestX[i]=Basic[i];
		}

		 return 0;
	}
// the last StartingPointsSize iterations we call DSO

	// here val contains the size of the box

//	*val ;
	if(*val>0/* && *val< (GetUpBound(0)-GetLowBound(0)) / 3.0 */) // call only for small regions
	for(i=0;i<m_basic;i++) {
			LeftT[i] = Basic[i]-   *val;
			RightT[i]= Basic[i]+   *val;
			if(LeftT[i]< GetLowBound(i)) LeftT[i]= GetLowBound(i);
			if(RightT[i]> GetUpBound(i)) RightT[i]= GetUpBound(i);
	} else *val=0;

	i = 0;
	ObjectiveF(Basic,&curr, &i); // evaluate atnthis point
	*val += m_box2 * m_box2* DSOPenalty;

//	FILE* fff;
	if(*val>0) {
//		fff=fopen("outganso.txt","a");
//		fprintf(fff,"%d %f\n",m_basic,val);
//		fprintf(fff,"%f %f %f\n",LeftT[0],LeftT[1],LeftT[2]);
//		fprintf(fff,"%f %f %f\n",RightT[0],RightT[1],RightT[2]);
//		fprintf(fff,"%f %f %f\n",Basic[0],Basic[1],Basic[2]);

//		fflush(fff);
//		fclose(fff);

	//	CallDSO(m_basic, Basic, LeftT, RightT, val,3,3); // try to reach another minimum
		if(*val > curr) *val=curr;  
	} 
		else *val=curr;
// record best value			
		if(*val<BestF)
		{
			BestF = *val;
			for(i=0;i<m_basic;i++) BestX[i]=Basic[i];
		}

	m_CurrentFunctionType=6; // was reset in CallDSO
	return 0;
}


// Called from DSO. Adds penalty on constraint violation for nonbasic variables
int		Ganso::FunctionDSO(int n, double * x, double* val, int* t)
{
	BuildFullVector(x);
	ObjFunction(&m_nvar,X,val,t);
	*val += m_box2 * m_box2* DSOPenalty;
	return 0;
}





// Linear transformation of variables
void	Ganso::ComputeXBasic(int n, double* x, double* Bas)
{
	// Bas = AforECAM * x + BforECAM
	int i,j;
	double t;
	for(i=0;i<m_basic;i++) {
		t=BforECAM[i];
		for(j=0;j<n;j++) t += AforECAM[i][j] * x[j];
		Bas[i]=t;
	}
}
	
//using n+1 points define a subspace spanned by them and use their CH as the domain 
// Here we define the matrix of the transformation
int		Ganso::PrepareSubspaceForECAM(int n, double* xb)
{
	// the original dim is m_basic
	AforECAM.newsize(m_basic,n);
	BforECAM.newsize(m_basic);
	int i,j;
	for(i=0;i<m_basic;i++) BforECAM[i]=xb[i];
	for(j=0;j<n;j++)
		for(i=0;i<m_basic;i++)
			AforECAM[i][j] = xb[j*m_basic + i] - xb[i];

// now ECAM should work in the unit cube/simplex - its basic domain

	return 0;
}

// shift the origin only
void	Ganso::ComputeXBasicShift(int n, double* x, double* Bas)
{
	// Bas =  x + BforECAM
	for(int i=0;i<m_basic;i++) 	Bas[i] = BforECAM[i] + x[i];
}
// calculate the shift
int		Ganso::PrepareForECAM(int n, double* xa)
{
	// the original dim is m_basic
	BforECAM.newsize(m_basic);
	for(int i=0;i<m_basic;i++) BforECAM[i]=xa[i];
	return 0;
}

/* This is an important step: define basic/nonbasic variables, index them, define
transformation matrices
*/

int		Ganso::SetConstraints(int lineq, int linineq, int nvar, double* AE, double* AI, 
							  double* RHSE, double * RHSI, double* Xl, double* Xu, int* basic)
{
#ifdef DEMO_VERSION
	if(nvar>100 || lineq>100 || linineq>100) return -200;
#endif
	// check that params >0
	if(lineq<0 || linineq<0) return -100;
	if(lineq>0 && (AE==NULL || RHSE==NULL) && Fortran==0) return -101;
	if(linineq>0 && (AI==NULL || RHSI==NULL) && Fortran==0) return -102;

	// construct the matrices of constraints
	m_lineq = lineq; m_linineq = linineq;  m_linconstr =  m_lineq + m_linineq; m_nvar=nvar;

	m_total = m_nvar + m_linineq; // slack variables
	m_basic = m_total - m_linconstr;  // basic vars
	m_nonbasic = m_total - m_basic;

	int flagbasic=0;
	if(m_basic<1) return -103;

	// we assume AE and AI are non-degenerate and the dimensions are correct

	TNT::Array2D<double> CONSTR(m_total,m_linconstr);
	TNT::Array2D<double> A2t(m_linconstr,m_linconstr);

	A1.newsize(m_linconstr,m_basic);
	A2.newsize(m_linconstr,m_linconstr);
	Bvec.newsize(m_linconstr);

	int i,j,k;

	for(i=0;i<m_linineq;i++) Bvec[i] = RHSI[i];
	for(i=0;i<m_lineq;i++) Bvec[i+m_linineq] = RHSE[i];  // the first are slack vars

	for(i=0;i<m_linineq;i++) {
		for(j=0;j<m_linineq;j++)
			if(i!=j) CONSTR[j][i] = 0; else CONSTR[j][i]=1;

		for(j=0;j<nvar;j++) CONSTR[j+m_linineq][i] = AI[i*m_nvar + j];
	}
	for(i=0;i<m_lineq;i++) {
		for(j=0;j<m_linineq;j++) CONSTR[j][i+m_linineq] = 0;
		for(j=0;j<nvar;j++) CONSTR[j + m_linineq][i + m_linineq] = AE[i*m_nvar + j];
	}


//	Array2D<double> CONSTRS(m_linconstr,m_total);
/*	for(i=0;i<m_total;i++)
		for(j=0;j<m_linconstr;j++) 
			 CONSTRS[j][i]=CONSTR[i][j];  // save a copy, transposed
*/

	TNT::Array1D<int> Pivot(m_total);
	if(m_linconstr == 0 ) {  // no linear constraints
		INdex.newsize(m_total);
		for(i=0;i<m_total;i++) INdex[i]=i;
		goto L2; // break
	} 

	if(Fortran==1) {
		if(basic==NULL || *basic==0) flagbasic=1;
		else for(i=0;i<m_basic;i++) basic[i]--;  // fortran indices changed to 0-based

	} else flagbasic = (basic==NULL);

	if(flagbasic &&  m_lineq>0) {
		LUP TempLUP(CONSTR);
		Pivot = TempLUP.getPivot();
		// check for errors?
	} else if(flagbasic &&  m_lineq==0) { // no equality constraints
		for(i=0;i<m_total;i++) Pivot[i]=i; 
	} else // index given
	{
		for(i=0;i<m_basic;i++) 
			Pivot[i+m_linconstr]=basic[i]+m_linineq;
		// populate the remaining
		k=0;
		for(i=0;i<m_total;i++) 
		{
			for(j=0;j<m_basic;j++) if(i == Pivot[j+m_linconstr]) goto L1;
			Pivot[k++]=i;
L1:;
		}
	}

	for(i=0;i<m_linconstr;i++) { // columns
		for(j=0;j<m_linconstr;j++)
			A2t[j][i] = CONSTR[Pivot[i]][j];
	}
	for(i=0;i<m_basic;i++) { // columns
		for(j=0;j<m_linconstr;j++)  // rows
			A1[j][i] = CONSTR[Pivot[i+m_linconstr]][j];
	}

	// factorize A2t

	MyLU = (void*) new LUP(A2t);
	if(((LUP*)MyLU)->isNonsingular() !=1) return -105; // singular transformation matrix

// now save the index
	INdex.newsize(m_total);
// currently I have slack, slack slack, x1 x2 ..x_nvar
// and the basic are in Pivot[m_linconstr + ...]
	for(i=0;i<m_basic;i++) {
		INdex[i] = Pivot[m_linconstr + i ] - m_linineq;
	    if(INdex[i]<0) // it was slack
			INdex[i]=nvar + Pivot[m_linconstr + i ]; // careful using these - the array dim
	}
	for(i=0;i<m_nonbasic;i++) {
		INdex[i + m_basic] = Pivot[ i ] - m_linineq;
	    if(INdex[i + m_basic]<0) // it was slack
			INdex[i + m_basic]=nvar + Pivot[ i ]; // careful using these - the array dim
	}

L2:
	Basic = new double[m_basic];
	X = new double[m_total];


// now box constraints --------------------------------------------
	m_BoxConstraints=0;
	if(m_linineq>0) m_BoxConstraints=1;
	for(i=0;i<m_nvar;i++) {
		if(Xl !=NULL) if(Xl[i]>-Infty) m_BoxConstraints=1;
		if(Xu !=NULL) if(Xu[i]<Infty) m_BoxConstraints=1;
		if(m_BoxConstraints) break;
	}
// if there are box constraints:
	if(m_BoxConstraints) {
		m_Xl=new double[m_total];
		m_Xu=new double[m_total];

		if(Xl !=NULL)
			for(i=0;i<m_nvar;i++) { if(Xl[i]>-Infty) m_Xl[i]=Xl[i]; else m_Xl[i]=-Infty;}
			else for(i=0;i<m_nvar;i++) m_Xl[i]=-Infty;
		if(Xu !=NULL)
			for(i=0;i<m_nvar;i++) { if(Xu[i]<Infty) m_Xu[i]=Xu[i]; else m_Xu[i]=Infty;}
			else for(i=0;i<m_nvar;i++) m_Xu[i]=Infty;
		

		for(i=m_nvar;i<m_total;i++) {m_Xu[i]=Infty; m_Xl[i]=0; }
	}


	return 0;
}

// useful for penalty functions 
int		Ganso::ConstraintViolations(double* x_basic, double* boxbasic, double* boxnonbasic, 
									double* slack)
{
	int i=BuildFullVector(x_basic);

	*boxbasic=m_box1;
	*boxnonbasic=m_box2;
	*slack=m_box3;
	return i;
}
	
void	Ganso::GetBasicVariables(double* x0, double* start)
{	int i;
	for(i=0;i<m_basic;i++) start[i]=x0[INdex[i]];
}

double max__(double a ,double b)
{ if(a>b) return a; else return b;}

int		Ganso::BuildFullVector(double* x_basic) // builds vector X
{
//#define max__(a,b) ((a)>(b)? (a):(b))

	A1D	Temp1(m_linconstr);

	int i,j;
	double t;
	// now if I have x_basic I populate X as
	for(i=0;i<m_basic;i++) 	X[INdex[i]] = x_basic[i];
	// generate non-basic

	if(m_linconstr<=0) goto L1;

	for(i=0;i<m_linconstr;i++) { // multiply by A1
		t=0;	
		for(j=0;j<m_basic;j++) 
			t += A1[i][j] * x_basic[j];
		Temp1[i] = Bvec[i] - t;
	}
	Temp1=((LUP*)MyLU)->solve(Temp1);
	for(i=0;i<m_linconstr;i++) 
		X[INdex[i+m_basic]] = Temp1[i];

L1:
// calculate box constraints violations
	 m_box1=m_box2=m_box3=0; // basic, non-basic, slack

	if(m_BoxConstraints) {
		for(i=0;i<m_basic;i++) {
			j=INdex[i];
			m_box1 += max__(0,-(X[j]-m_Xl[j])) + max__(0,X[j]-m_Xu[j]);
		}
		for(i=0;i<m_linconstr;i++) {  
			j=INdex[i+m_basic];
			m_box2 +=  max__(0, -(X[j] - m_Xl[j])) + max__(0,X[j] - m_Xu[j]);
		}

// total box constr violations is m_box1+m_box2

		for(i=m_nvar;i<m_total;i++)
			m_box3 += max__(0, -(X[i] - m_Xl[i])) + max__(0,X[i] - m_Xu[i]);
	}


	return 0;
}


// this is an auxiliary class implementing random search
RandomSearch::RandomSearch(){
	p_QRan=new qrandom;
}

RandomSearch::~RandomSearch(){
	delete (qrandom*)p_QRan;
}

int	RandomSearch::SearchSimplex(int dim, double * sol, double * val, int iter)
{
	int i,j;
	double  v,best=10e20;
	double* temp=(double*)malloc(sizeof(double)*(dim+1));
	qrandom	MyQRan;

	for(i=0;i<iter;i++) {
			MyQRan.RandomVecSimplex(dim+1,temp);
			Parent->CurrentFunction(dim, temp, &v, 0);
			if(v<best) {
					best=v;
					for(j=0;j<dim;j++) sol[j]=temp[j];
			}
	}
	*val=best;
	return 0;
}
	
int RandomSearch::GenerateOnCube(int dim, double* sol, double* left, double* right)
{
	int i;
	if(dim==1) {
		sol[0]=rand()/(double)RAND_MAX;
	} else
		((qrandom*)(p_QRan))->RandomVec(dim,sol);

	for(i=0;i<dim;i++) sol[i] = sol[i]*(right[i]-left[i]) + left[i];
	return 0;
}


/*=======service routines========================*/
	double	Ganso::GetWidthMin()
	{	int i;
		double t,r=2e6;
		if(m_BoxConstraints)
		for(i=0;i<m_basic;i++) {
			t=m_Xu[INdex[i]]-m_Xl[INdex[i]];
			if(r>t) r=t;
		}
		return r;
	};

	double Ganso::GetLowBound(int i) {
		if(m_BoxConstraints)
		return m_Xl[INdex[i]];
		else return -10e10;
	};
	double Ganso::GetUpBound(int i) {
		if(m_BoxConstraints)
		return m_Xu[INdex[i]];
		else return 10e10;
	};


/*===   procedural interface   ====================================*/

GANSODLL_API int  MinimizeDFBM(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter)
{	Ganso G;
	return G.MinimizeDFBM( dim,  x0, val,  f, 
		 lineq,  linineq,  AE,  AI,  RHSE,  RHSI, Xl,  Xu, basic, maxiter);}

GANSODLL_API int  MinimizeDFBMECAM(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter, int iterECAM , int dimECAM)
{	Ganso G;
	return G.MinimizeDFBMECAM( dim,  x0, val,  f, 
		 lineq,  linineq,  AE,  AI,  RHSE,  RHSI, Xl,  Xu, basic, maxiter,  iterECAM ,  dimECAM);}

GANSODLL_API int  MinimizeRandomStart(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter)
{	Ganso G;
	return G.MinimizeRandomStart( dim,  x0, val,  f, 
		 lineq,  linineq,  AE,  AI,  RHSE,  RHSI, Xl,  Xu, basic, maxiter);}

GANSODLL_API int  MinimizeECAM(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double LC, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter, int iterLocal)
{	Ganso G;
	return G.MinimizeECAM( dim,  x0, val,  f, 
		 lineq,  linineq,  LC, AE,  AI,  RHSE,  RHSI, Xl,  Xu, basic, maxiter,  iterLocal);}


GANSODLL_API int  MinimizeECAMDFBM(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double LC, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter)
{	Ganso G;
	return G.MinimizeECAMDFBM( dim,  x0, val,  f, 
		 lineq,  linineq,  LC, AE,  AI,  RHSE,  RHSI, Xl,  Xu, basic, maxiter);}


#ifdef USEDSO

GANSODLL_API int  MinimizeDSO(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic,  double penalty, int maxiter, int steps)
{	Ganso G;
	return G.MinimizeDSO( dim,  x0, val,  f, 
		 lineq,  linineq,  AE,  AI,  RHSE,  RHSI, Xl,  Xu, basic, penalty, maxiter,  steps);}

GANSODLL_API int  MinimizeIterativeDSO(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic,  double penalty, int maxiter, int steps, int epoch)
{	Ganso G;
	return G.MinimizeIterativeDSO( dim,  x0, val,  f, 
		 lineq,  linineq,  AE,  AI,  RHSE,  RHSI, Xl,  Xu, basic, penalty, maxiter,  steps,  epoch);}


GANSODLL_API int  MinimizeECAMDSO(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double LC, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter, int iterDSO)
{	Ganso G;
	return G.MinimizeECAMDSO( dim,  x0, val,  f, 
		 lineq,  linineq,  LC, AE,  AI,  RHSE,  RHSI, Xl,  Xu, basic, maxiter, iterDSO);}

#endif



GANSODLL_API int  MinimizeDFBM_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl, double* Xu, int maxiter)
{ 	Ganso G;
	return G.MinimizeDFBM_0( dim,  x0, val,  f,  Xl,  Xu,  maxiter);}


GANSODLL_API 	int		MinimizeDFBMECAM_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl, double* Xu,  int maxiter, int iterECAM , int dimECAM)
{	Ganso G;
	return G.MinimizeDFBMECAM_0( dim,  x0, val,  f,  Xl,  Xu,  maxiter,iterECAM,dimECAM);}

GANSODLL_API 	int		MinimizeRandomStart_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl, double* Xu, int maxiter)
{	Ganso G;
	return G.MinimizeRandomStart_0( dim,  x0, val,  f,  Xl,  Xu,  maxiter);}


GANSODLL_API 	int		MinimizeECAM_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double LC, double* Xl, double* Xu, int maxiter, int iterLocal)
{	Ganso G;
	return G.MinimizeECAM_0( dim,  x0, val,  f, LC,  Xl,  Xu,  maxiter,  iterLocal);}

GANSODLL_API 	int		MinimizeECAMDFBM_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		 double LC, double* Xl, double* Xu,  int maxiter)
{	Ganso G;
	return G.MinimizeECAMDFBM_0( dim,  x0, val, f,LC,   Xl,  Xu,  maxiter);}


#ifdef USEDSO
GANSODLL_API 	int		MinimizeDSO_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl, double* Xu, int maxiter, int steps)
{	Ganso G;
	return G.MinimizeDSO_0( dim,  x0, val,  f,  Xl,  Xu,  maxiter,steps);}

GANSODLL_API 	int		MinimizeIterativeDSO_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl, double* Xu, int maxiter, int steps, int epoch)
{	Ganso G;
	return G.MinimizeIterativeDSO_0( dim,  x0, val,  f,  Xl,  Xu,  maxiter,steps, epoch);}

GANSODLL_API 	int		MinimizeECAMDSO_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double LC, 	double* Xl, double* Xu, int maxiter, int iterDSO)
{	Ganso G;
	return G.MinimizeECAMDSO_0( dim,  x0, val,  f, LC, Xl,  Xu,  maxiter, iterDSO);}
#endif




USER_FUNCTION_C UserF = 0;
#ifdef _MSC_VER
void mystdcall STDCALLUSER(int* n, double* x, double* val, int* t){ UserF(n,x,val, t);}
#else
inline void  STDCALLUSER(int* n, double* x, double* val, int* t){ UserF(n,x,val,t);}
#endif
//=========================================================================
// interface for calls from FORTRAN
GANSODLL_API int  minimizedfbm(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		int* lineq, int* linineq, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int* maxiter, int* boxconstraints)
{	Ganso G;
	G.Fortran=1;
	UserF=f;

	if(*maxiter==0 )*maxiter=DEFAULT_IterDGM;

	if(*boxconstraints==0)
	return G.MinimizeDFBM( *dim,  x0, val,  STDCALLUSER, 
		 *lineq,  *linineq,  AE,  AI,  RHSE,  RHSI, NULL,  NULL, basic, *maxiter);
	else if (*boxconstraints==1)
	return G.MinimizeDFBM( *dim,  x0, val,  STDCALLUSER, 
		 *lineq,  *linineq,  AE,  AI,  RHSE,  RHSI, Xl,  NULL, basic, *maxiter);
	else  if( *boxconstraints==2)
	return G.MinimizeDFBM( *dim,  x0, val,  STDCALLUSER, 
		 *lineq,  *linineq,  AE,  AI,  RHSE,  RHSI, NULL,  Xu, basic, *maxiter);
	else 
	return G.MinimizeDFBM( *dim,  x0, val,  STDCALLUSER, 
		 *lineq,  *linineq,  AE,  AI,  RHSE,  RHSI, Xl,  Xu, basic, *maxiter);
}

GANSODLL_API 	int		minimizedfbmecam(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		int* lineq, int* linineq, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int* maxiter, int* iterECAM , int* dimECAM, int* boxconstraints)
{	Ganso G;
	G.Fortran=1;
	UserF=f;

	if(*maxiter==0 )*maxiter=DEFAULT_IterDGM;

	if(*boxconstraints==0)
	return G.MinimizeDFBMECAM( *dim,  x0, val,  STDCALLUSER, 
		 *lineq,  *linineq,  AE,  AI,  RHSE,  RHSI, NULL,  NULL, basic, *maxiter, *iterECAM ,  *dimECAM);
	else if (*boxconstraints==1)
	return G.MinimizeDFBMECAM( *dim,  x0, val,  STDCALLUSER, 
		 *lineq,  *linineq,  AE,  AI,  RHSE,  RHSI, Xl,  NULL, basic, *maxiter, *iterECAM ,  *dimECAM);
	else if( *boxconstraints==2)
	return G.MinimizeDFBMECAM( *dim,  x0, val,  STDCALLUSER, 
		 *lineq,  *linineq,  AE,  AI,  RHSE,  RHSI, NULL,  Xu, basic, *maxiter, *iterECAM ,  *dimECAM);
	else
	return G.MinimizeDFBMECAM( *dim,  x0, val,  STDCALLUSER, 
		 *lineq,  *linineq,  AE,  AI,  RHSE,  RHSI, Xl,  Xu, basic, *maxiter, *iterECAM ,  *dimECAM);
}

GANSODLL_API 	int	 	minimizerandomstart(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		int* lineq, int* linineq, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int* maxiter)
{	Ganso G;
	UserF=f;
	G.Fortran=1;

	if(*maxiter==0 )*maxiter=1;

	return G.MinimizeRandomStart( *dim,  x0, val,  STDCALLUSER, 
		 *lineq,  *linineq,  AE,  AI,  RHSE,  RHSI, Xl,  Xu, basic, *maxiter);}

GANSODLL_API 	int	 	minimizeecam(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		int* lineq, int* linineq, double* LC, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int* maxiter, int* iterLocal)
{	Ganso G;
	UserF=f;
	G.Fortran=1;

	if(*maxiter==0 )*maxiter=DEFAULT_IterECAM;

	return G.MinimizeECAM( *dim,  x0, val,  STDCALLUSER, 
		 *lineq,  *linineq,  *LC, AE,  AI,  RHSE,  RHSI, Xl,  Xu, basic, *maxiter, *iterLocal);}

GANSODLL_API 	int	 	minimizeecamdfbm(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		int* lineq, int* linineq, double* LC, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int* maxiter)
{	Ganso G;
	UserF=f;
	G.Fortran=1;

	if(*maxiter==0 )*maxiter=DEFAULT_IterECAM;

	return G.MinimizeECAMDFBM( *dim,  x0, val,  STDCALLUSER, 
		 *lineq,  *linineq,  *LC, AE,  AI,  RHSE,  RHSI, Xl,  Xu, basic, *maxiter);}

#ifdef USEDSO
GANSODLL_API 	int	 	minimizedso(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		int* lineq, int* linineq,  double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic,  double* penalty, int* maxiter, int* steps)
{	Ganso G;
	UserF=f;
	G.Fortran=1;
	return G.MinimizeDSO( *dim,  x0, val,  STDCALLUSER, 
		 *lineq,  *linineq,  AE,  AI,  RHSE,  RHSI, Xl,  Xu, basic, *penalty, *maxiter,  *steps);}

GANSODLL_API 	int	 	minimizeiterativedso(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		int* lineq, int* linineq,  double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic,  double* penalty, int* maxiter, int* steps, int* epoch)
{	Ganso G;
	UserF=f;
	G.Fortran=1;
	return G.MinimizeIterativeDSO( *dim,  x0, val,  STDCALLUSER, 
		 *lineq,  *linineq,  AE,  AI,  RHSE,  RHSI, Xl,  Xu, basic, *penalty, *maxiter,  *steps,  *epoch);}

GANSODLL_API 	int	 	minimizeecamdso(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		int* lineq, int* linineq, double* LC, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int* maxiter, int *iterDSO)
{	Ganso G;
	UserF=f;
	G.Fortran=1;

	if(*maxiter==0 )*maxiter=DEFAULT_IterECAM;

	return G.MinimizeECAMDSO( *dim,  x0, val,  STDCALLUSER, 
		 *lineq,  *linineq,  *LC, AE,  AI,  RHSE,  RHSI, Xl,  Xu, basic, *maxiter,  *iterDSO);}

#endif


GANSODLL_API int  minimizedfbm_0(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		double* Xl, double* Xu, int* maxiter, int* boxconstraints)
{ 	Ganso G;
	UserF=f;
	G.Fortran=1;
	if(*maxiter==0 )*maxiter=DEFAULT_IterDGM;

	if(*boxconstraints==0)
	return G.MinimizeDFBM_0( *dim,  x0, val,  STDCALLUSER,   NULL,  NULL, *maxiter);
	else if (*boxconstraints==1)
	return G.MinimizeDFBM_0( *dim,  x0, val,  STDCALLUSER,   Xl,  NULL, *maxiter);
	else if (*boxconstraints==2)
	return G.MinimizeDFBM_0( *dim,  x0, val,  STDCALLUSER,   NULL,  Xu, *maxiter);
	else 
	return G.MinimizeDFBM_0( *dim,  x0, val,  STDCALLUSER,   Xl,  Xu, *maxiter);
//	return G.MinimizeDFBM_0( *dim,  x0, val,  UserF,  Xl,  Xu,  *maxiter);
}



GANSODLL_API 	int		minimizedfbmecam_0(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		double* Xl, double* Xu,  int* maxiter, int* iterECAM , int* dimECAM, int* boxconstraints)
{	Ganso G;
	UserF=f;
	G.Fortran=1;
	if(*maxiter==0 )*maxiter=DEFAULT_IterDGM;

	if(*boxconstraints==0)
	return G.MinimizeDFBMECAM_0( *dim,  x0, val,  STDCALLUSER,   NULL,  NULL,  *maxiter, *iterECAM ,  *dimECAM);
	else if (*boxconstraints==1)
	return G.MinimizeDFBMECAM_0( *dim,  x0, val,  STDCALLUSER,   Xl,  NULL,  *maxiter, *iterECAM ,  *dimECAM);
	else if (*boxconstraints==2)
	return G.MinimizeDFBMECAM_0( *dim,  x0, val,  STDCALLUSER,   NULL,  Xu, *maxiter, *iterECAM ,  *dimECAM);
	else 
	return G.MinimizeDFBMECAM_0( *dim,  x0, val,  STDCALLUSER,   Xl,  Xu, *maxiter, *iterECAM ,  *dimECAM);
//	return G.MinimizeDFBMECAM_0( *dim,  x0, val,  f,  Xl,  Xu,  *maxiter,*iterECAM,*dimECAM);
}


GANSODLL_API 	int	 	minimizerandomstart_0(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		double* Xl, double* Xu, int* maxiter)
{	Ganso G;
	UserF=f;
	G.Fortran=1;
	if(*maxiter==0 )*maxiter=1;

	return G.MinimizeRandomStart_0( *dim,  x0, val,  STDCALLUSER,  Xl,  Xu,  *maxiter);}


GANSODLL_API 	int	 	minimizeecam_0(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		double* LC, double* Xl, double* Xu, int* maxiter, int* iterLocal)
{	
	Ganso G;
	UserF=f;
	G.Fortran=1;
	if(*maxiter==0 )*maxiter=DEFAULT_IterECAM;

	return G.MinimizeECAM_0( *dim,  x0, val,  STDCALLUSER, *LC,  Xl,  Xu,  *maxiter, *iterLocal);}



GANSODLL_API 	int	 	minimizeecamdfbm_0(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		 double* LC, double* Xl, double* Xu,  int* maxiter)
{	Ganso G;
	UserF=f;
	G.Fortran=1;
	if(*maxiter==0 )*maxiter=DEFAULT_IterECAM;

	return G.MinimizeECAMDFBM_0( *dim,  x0, val, STDCALLUSER,*LC,   Xl,  Xu,  *maxiter);}


#ifdef USEDSO
GANSODLL_API 	int	 	minimizedso_0(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		double* Xl, double* Xu, int* maxiter, int* steps)
{	Ganso G;
	UserF=f;
	G.Fortran=1;
	return G.MinimizeDSO_0( *dim,  x0, val,  STDCALLUSER,  Xl,  Xu,  *maxiter,*steps);}


GANSODLL_API 	int	 	minimizeiterativedso_0(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		double* Xl, double* Xu, int* maxiter, int* steps, int* epoch)
{	Ganso G;
	UserF=f;
	G.Fortran=1;
	return G.MinimizeIterativeDSO_0( *dim,  x0, val,  STDCALLUSER,  Xl,  Xu,  *maxiter,*steps, *epoch);}


GANSODLL_API 	int	 	minimizeecamdso_0(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		double* LC, 	double* Xl, double* Xu, int* maxiter, int *iterDSO)
{	Ganso G;
	UserF=f;
	G.Fortran=1;
	if(*maxiter==0 )*maxiter=DEFAULT_IterECAM;

	return G.MinimizeECAMDSO_0( *dim,  x0, val,  STDCALLUSER, *LC, Xl,  Xu,  *maxiter,  *iterDSO);}

#endif


