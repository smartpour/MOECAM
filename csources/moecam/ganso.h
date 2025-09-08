// use this only if we are compiling a dll
#ifdef _MSC_VER 
//_USRDLL


#ifdef GANSODLL_EXPORTS
#define GANSODLL_API __declspec(dllexport)
#include "../gansodll/stdafx.h"
#else
#define GANSODLL_API __declspec(dllimport)
#endif

#else  // just an empty string if not using dll (ie for linux)
#define GANSODLL_API
#endif

#define GANSODLL_API

#define compilecdecl
#ifdef compilecdecl
#define mystdcall __cdecl
#else 
#define mystdcall __stdcall
#endif

#ifdef __cplusplus
#include "gansoc.h"

extern "C" {
#else
typedef void ( *USER_FUNCTION)(int *, double *, double *);

#ifdef _MSC_VER
typedef void (__cdecl *USER_FUNCTION_C)(int *, double *, double *);
#else
typedef void ( *USER_FUNCTION_C)(int *, double *, double *);
#endif

#endif

// procedural interface to the members of Ganso class


GANSODLL_API int  MinimizeDFBM(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter);

GANSODLL_API 	int		MinimizeDFBMECAM(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter, int iterECAM , int dimECAM);

GANSODLL_API 	int	 	MinimizeRandomStart(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter);

GANSODLL_API 	int	 	MinimizeECAM(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double LC, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter, int iterLocal);

GANSODLL_API 	int	 	MinimizeECAMDFBM(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double LC, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter);

#ifdef USEDSO
GANSODLL_API 	int	 	MinimizeDSO(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq,  double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, double penalty, int speed, int precision);

GANSODLL_API 	int	 	MinimizeIterativeDSO(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq,  double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, double penalty, int speed, int precision, int rounds);

GANSODLL_API 	int	 	MinimizeECAMDSO(int dim, double* x0, double *val, USER_FUNCTION f, 
		int lineq, int linineq, double LC, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int maxiter, int iterDSO);
#endif


GANSODLL_API int  MinimizeDFBM_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl, double* Xu, int maxiter);

GANSODLL_API 	int		MinimizeDFBMECAM_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl, double* Xu,  int maxiter, int iterECAM , int dimECAM);

GANSODLL_API 	int	 	MinimizeRandomStart_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl, double* Xu, int maxiter);


GANSODLL_API 	int	 	MinimizeECAM_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double LC, double* Xl, double* Xu, int maxiter, int iterLocal);


GANSODLL_API 	int	 	MinimizeECAMDFBM_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		 double LC, double* Xl, double* Xu,  int maxiter);

#ifdef USEDSO
GANSODLL_API 	int	 	MinimizeDSO_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl, double* Xu, int speed, int precision);

GANSODLL_API 	int	 	MinimizeIterativeDSO_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double* Xl, double* Xu, int speed, int precision, int rounds);

GANSODLL_API 	int	 	MinimizeECAMDSO_0(int dim, double* x0, double *val, USER_FUNCTION f, 
		double LC, 	double* Xl, double* Xu, int maxiter, int iterDSO);
#endif


//=========================================================================
// interface for calls from FORTRAN
GANSODLL_API int  minimizedfbm(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		int* lineq, int* linineq, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int* maxiter, int* boxconstraints);

GANSODLL_API 	int		minimizedfbmecam(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		int* lineq, int* linineq, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int* maxiter, int* iterECAM , int* dimECAM, int* boxconstraints);

GANSODLL_API 	int	 	minimizerandomstart(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		int* lineq, int* linineq, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int* maxiter);

GANSODLL_API 	int	 	minimizeecam(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		int* lineq, int* linineq, double* LC, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int* maxiter, int *iterLocal);

GANSODLL_API 	int	 	minimizeecamdfbm(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		int* lineq, int* linineq, double* LC, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int* maxiter);

#ifdef USEDSO
GANSODLL_API 	int	 	minimizedso(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		int* lineq, int* linineq,  double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic,  double* penalty, int* speed, int* precision);

GANSODLL_API 	int	 	minimizeiterativedso(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		int* lineq, int* linineq,  double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, double* penalty, int* speed, int* precision, int* rounds);

GANSODLL_API 	int	 	minimizeecamdso(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		int* lineq, int* linineq, double* LC, double* AE, double* AI, double* RHSE, double * RHSI,
		double* Xl, double* Xu, int* basic, int* maxiter, int *iterDSO);
#endif


GANSODLL_API int  minimizedfbm_0(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		double* Xl, double* Xu, int* maxiter, int* boxconstraints);

GANSODLL_API 	int		minimizedfbmecam_0(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		double* Xl, double* Xu,  int* maxiter, int* iterECAM , int* dimECAM, int* boxconstraints);

GANSODLL_API 	int	 	minimizerandomstart_0(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		double* Xl, double* Xu, int* maxiter);


GANSODLL_API 	int	 	minimizeecam_0(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		double* LC, double* Xl, double* Xu, int* maxiter, int *iterLocal);


GANSODLL_API 	int	 	minimizeecamdfbm_0(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		 double* LC, double* Xl, double* Xu,  int* maxiter);

#ifdef USEDSO
GANSODLL_API 	int	 	minimizedso_0(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		double* Xl, double* Xu, int* speed, int* precision);

GANSODLL_API 	int	 	minimizeiterativedso_0(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		double* Xl, double* Xu, int* speed, int* precision, int* rounds);

GANSODLL_API 	int	 	minimizeecamdso_0(int* dim, double* x0, double *val, USER_FUNCTION_C f, 
		double* LC, 	double* Xl, double* Xu, int* maxiter, int *iterDSO);
#endif





#ifdef __cplusplus
}
#endif


