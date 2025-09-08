

#include <cstdlib>

// this is needed for compatibility with lf95
#ifdef _MSC_VER
typedef void (__cdecl *USER_FUNCTION_C)(int *, double *, double *, int *);
#else
typedef void ( *USER_FUNCTION_C)(int *, double *, double *, int *);
#endif

typedef void ( *USER_FUNCTION)(int *, double *, double *, int* );

#define  Infty 10e20
#define DEFAULT_IterDGM 10000
#define DEFAULT_IterECAM 50

//#define USEDSO
#define USEECAM
#define GANSODLL_API


class SawToothCover;

class GANSODLL_API  Ganso {

public:
//  the number of linear equality, intequality and total linear constraints
//  the number of original variables, basic and nonbasic variables, total (ie+ slack)
	int		m_lineq, m_linineq, m_linconstr, m_nvar, m_basic, m_nonbasic, m_total;

	double *X;				// vector to be passed to the user's obj function
	double *Basic;			// aux vector of basic variables
	double *m_Xl, *m_Xu;		// upper and lower bounds
	int	    m_BoxConstraints;   // 1 of there are box constraints, 0 otherwise
	double  m_box1, m_box2, m_box3; // basic, non-basic, slack
	double	DSOPenalty;			    // defaults to 100, may be changed by the user
	int		Fortran;		// flag, if called from Fortran, for parameter basic 

	int		FunctionCounter, FunCallsWhenFound;
	double  ECAMLowerBound;

	int		m_objectives;
	SawToothCover* MySTC;

private:
// these are placeholders for transformation matrices, used internally
	void *MyLU;
	void *a2, *a1, *bvec;
	void *Index;
	void *aforECAM, *bforECAM; 

	double *LeftT, *RightT;
	double *BestX, BestF;
	double *StartingPoints;
	int		StartingPointsSize, Iterations, curriter;

	USER_FUNCTION		ObjFunction;	// pointer to the user's objective function
//	USER_FUNCTION_C		ObjFunction_C;	// pointer to the user's objective function __cdecl

public:
	Ganso();
	~Ganso();
// The interface to GANSO optimization methods. See the manual for their description
	int		MinimizeDFBM(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double* AE=NULL, double* AI=NULL, double* RHSE=NULL, double * RHSI=NULL,
		double* Xl=NULL, double* Xu=NULL, int* basic=NULL, int maxiter=0);

	int		MinimizeDFBMECAM(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double* AE=NULL, double* AI=NULL, double* RHSE=NULL, double * RHSI=NULL,
		double* Xl=NULL, double* Xu=NULL, int* basic=NULL, int maxiter=0, int iterECAM =0, int dimECAM=0);

	int		MinimizeRandomStart(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double* AE=NULL, double* AI=NULL, double* RHSE=NULL, double * RHSI=NULL,
		double* Xl=NULL, double* Xu=NULL, int* basic=NULL, int maxiter=0);


#ifdef USEECAM
	int		MinimizeECAM(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double LC=100, double* AE=NULL, double* AI=NULL, double* RHSE=NULL, double * RHSI=NULL,
		double* Xl=NULL, double* Xu=NULL, int* basic=NULL, int maxiter=100, int iterLocal=1);

	int		MinimizeECAMDFBM(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double LC=100, double* AE=NULL, double* AI=NULL, double* RHSE=NULL, double * RHSI=NULL,
		double* Xl=NULL, double* Xu=NULL, int* basic=NULL, int maxiter=100);

#endif


#ifdef USEDSO
	int		MinimizeDSO(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq,  double* AE=NULL, double* AI=NULL, double* RHSE=NULL, double * RHSI=NULL,
		double* Xl=NULL, double* Xu=NULL, int* basic=NULL, double penalty=10,int speed=3, int precision = 4);

	int		MinimizeIterativeDSO(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq,  double* AE=NULL, double* AI=NULL, double* RHSE=NULL, double * RHSI=NULL,
		double* Xl=NULL, double* Xu=NULL, int* basic=NULL, double penalty=10,int speed=3, int precision = 4, int rounds=1);

#ifdef USEECAM
	int		MinimizeECAMDSO(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double LC=100, double* AE=NULL, double* AI=NULL, double* RHSE=NULL, double * RHSI=NULL,
		double* Xl=NULL, double* Xu=NULL, int* basic=NULL, int maxiter=100, int iterDSO=1);
#endif
#endif

// simplified versions of the above methods, with no linear constraints


	int		MinimizeDFBM_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl=NULL, double* Xu=NULL, int maxiter=0);

	int		MinimizeDFBMECAM_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl=NULL, double* Xu=NULL,  int maxiter=0, int iterECAM =0, int dimECAM=0);

	int		MinimizeRandomStart_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl=NULL, double* Xu=NULL, int maxiter=0);

#ifdef USEECAM
	int		MinimizeECAM_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double LC=100, double* Xl=NULL, double* Xu=NULL, int maxiter=100, int iterLocal=1);

	int		MinimizeECAMDFBM_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		 double LC=100, double* Xl=NULL, double* Xu=NULL,  int maxiter=100);
#endif

#ifdef USEDSO
	int		MinimizeDSO_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl=NULL, double* Xu=NULL,  int speed=3, int precision = 4);

	int		MinimizeIterativeDSO_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl=NULL, double* Xu=NULL,  int speed=3, int precision = 4, int rounds=1);

#ifdef USEECAM
	int		MinimizeECAMDSO_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double LC=100, 	double* Xl=NULL, double* Xu=NULL, int maxiter=100, int iterDSO=1);
#endif
#endif


// called by the minimization algorithms. Calls the appropriate obj fucntion, depending on the method used
	int		CurrentFunction(int n, double* x, double* val, int* t);
	int		ConstraintViolations(double* x_basic, double* boxbasic, double* boxnonbasic, double* slack);

/*=======service routines========================*/
	double GetWidthMin();
	double GetLowBound(int i);
	double GetUpBound(int i);

	int	DecideObjectives();
	void MakeStepbefore();
	void MakeStep(int objchosen);

// these methods are called internally

#ifdef USEECAM
	int		CallECAM ( int dim,   double* sol, double* val,  double LC, int iterECAM, int iterLocal,
					 int useLocal);
#endif

	int		CallECAMfromDFBM (
		int dim, double* basis, double LC, int iterECAM, double* sol, double* fval		);
// call ECAM with a special function in the unit domain

	int		CallDFBM(int dim, double* startpoint, double* val, 
		int UseECAM=0, int maxiter=0, int iterECAM=0, int dimECAM =0);

	int	CallDSO (int dim,   double* sol, double* L,double* R, double* val, int iterDSO=3, int stepsDSO=4, int rounds=1);

// some statistics 
	int		GetFunctionCalls() {return FunctionCounter;}
#ifdef USEECAM
	double		GetECAMLowerBound(){ return ECAMLowerBound;}
	int		GetFunctionCallsSolFound() {return FunCallsWhenFound;}
#endif

private:
// specifies what coordinate transformation will be done
	int		m_CurrentFunctionType;
// below are the options

// Called from within DFBM, x is in the subspace of dimension n
	int		FunctionECAMSubspace(int n, double * x, double* val);

// Called by DGM, transforms x
	int		FunctionDFBM(int n, double * x, double* val);

// Called from ECAM, shifts and then transforms x 
	int		FunctionECAM(int n, double * x, double* val, int* t);

// Called from ECAM, shifts x, then calls DFBM from this starting point
	int		FunctionECAMLocal(int n, double * x, double* val);

// Called from ECAM, shifts x, then calls DSO from this starting point
	int		FunctionECAMDSO(int n, double * x, double* val);


// Called from DSO, uses penalty
	int		FunctionDSO(int n, double * x, double* val, int* t);

// Calls the user supplied obj. function directly
	int		ObjectiveF(double * x, double* val, int* t);


// determines basic ans slack variables, creates matrices for transformation, and the index vector
	int		SetConstraints(int lineq, int linineq, int nvar, double* AE, double* AI, double* RHSE, double * RHSI, 
			double* Xl=NULL, double* Xu=NULL, int* basic=NULL);

	int		BuildFullVector(double* x_basic); // builds vector X from the basic variables
	void	GetBasicVariables(double* x0, double* start);

	int		SetObjFunction(USER_FUNCTION f) {FunctionCounter=0; 
	m_CurrentFunctionType=0; ObjFunction=f; 
	if(f!=NULL) return 0; else return -1;};

// prepares matrix AforECAM and vector BforECAM, n is the dimnsion of the subspace
	int		PrepareSubspaceForECAM(int n, double* xb);

// transforms x to Basic variables, by multiplying and shifting it Bas = AforECAM x + BforECAM
	void	ComputeXBasic(int n, double* x, double* Bas);

// prepares the shift vector BforECAM
	int		PrepareForECAM(int n, double* xa);

// transforms x to Basic variables, by shifting it
	void	ComputeXBasicShift(int n, double* x, double* Bas);

};

// auxiliary class that perform random search

class	RandomSearch 
{
public:
	RandomSearch();
	~RandomSearch();
	Ganso* Parent;
	void* p_QRan;
	int SearchSimplex(int dim, double * sol, double * val, int iter); 
	int GenerateOnCube(int dim, double* sol, double* left, double* right);

};
