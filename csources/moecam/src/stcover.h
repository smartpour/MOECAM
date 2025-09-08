// cadg with heap

#ifndef _STCOVER
#define _STCOVER

#include "ecam.h"



class SawToothCover;



// Class Saw tooth cover--------------------------------------------------
// this is our main object
/*---------------------------------
	Main algorithm:
	consists in growing the tree (forest) of SVSetNode by splitting leafs failing 
	cond (2). The top leaf (glob. minimum of aux function, largest SVSet.Dval), 
	although not necessarily, could be any other point, generates a new SV.
	The tree is processed from the root. Only branches that fail cond (2) are processed,
	because for a leaf to fail (2), parents must also fail (2)

	There can be virtual branches (when splitting is delayed, and possible killer is kept
	in the list killers). Before generating SV, ensure that we actually have a leaf,
	not virtual branch. ProcessFront will return only once this is ensured.

	1. Call ProcessFront to get the leaf with best Dval
	2. Call AssessCandidate to guarantee separation from the boundary (+ constraints)
	3. Call FormNewSV - evaluate f() and add new SV
	4a. Call ProcessTopwithSV to split the SVSet generated this SV (for efficiency)
	4b. Call ProcessAll with this SV to test cond (2) on the rest of the tree
	 This will generate virtual branches, which will be split once they are at front

-----------------------------*/

class	SawToothCover {
public:
// where to put the results
	FILE *fp;
	FILE *logfile;

	Problem *Prob;

	int retrievedump;
	char dumpfilename[200];
    SVList  PreparedPointsList;

	SawToothCover();
	~SawToothCover();

	void InitProblem(Problem * _prob) {Prob=_prob;}

	SVDeque				SVectors;	// here support vectors live
	SVDeque::iterator	iter;
	SVDeque::iterator	t_iter; // just iterator

	SVXDeque ValuesFromFile;


// SVSets live here
	Forest				HeapPossibleMin;

// aux staff
	siVec				m_lastindex, m_TempSiVec;
	real				tempFval;
	int					Match;

// solution sits here
	dVec	xglob_x, temp_x;
	real	xglob_val, temp_val;
	double*	FinalList;

	int			iteration, Iters, ListSize;

	int			LastLabel;
	int			Dim; // problem size

	dVec	m_Constants;    // here Lipschitz constants live
	real	C_minus;
	dVec	m_bounds;
	double	m_SimplexBound;
	int		m_myobj;

	double	GVA, frecord, flast;

	double Cut ;
	double Const ;

#ifdef MINLP
	iVec		m_intvariables;
	int			m_pureinteger;
	int			m_numberofintegers;
#endif
	double  radius;
	int		FunCallsWhenFound;

	clock_t clockS,clockF;
	double duration;

	double RequestGap();
	// value of the objectiove function
	double	Value(dVec &x);
	// value of the objectiove function, cut to Cut and with added Const
	double	ValueNormalise(double val);

	void	InitPopulate(); // vertices of simplex
	void	AddVertices();

	int		AllowVirtual();
// check for constraints / closeness to the boundary, and prepare for FormSVector
	int		AssessCandidate(SVSetNode* S, SVSetNodePtr thisnode);

// using the good candidate, compute f() (or retrieve it if in parallel)
	void FormNewSV();


	void InitialSteps(int d, int iters);
	int MakeStep(dVec&  X);
	int MakeStepbefore();



	int		StopCriteria();
	void	Minimise(int code=0);  //main routine

// 
	void	LoadAdditionalPoints();
	void	LoadAdditionalPoints(SVDeque& pts);


// the rest are auxiliary routines	
	void PrintSolution();
	void PrintSolution(char *filename);
	void CallLocalMinimiser(dVec& x, double& val, dVec& sol);
	void CallLocalMinimiserExt(dVec& x, double& val, dVec& sol);
	int PopNextCandidate(dVec& x);
	void PrepareFinalList(double* List, int* size);

	int PopNextCandidateExternal(dVec& x);
	void	ReturnSol(dVec &x);

	void	PrintPossibleMinima();

// dump into file routines
	void Dump(char *filename);
	void RetrieveDump(char *filename);
	void DumpBin(char *filename);
	void RetrieveDumpBin(char *filename);



// the aux. methods

//	void SetBoxConstraints(dVec&a, dVec& b);
	void SetConst(double newconst);
	void SetCut(double newcut);
	void SetIter(int iter);
	double GetConst();
	double GetCut();
	void  SetListSize(int m);	

	void	SetConstants(dVec& newconst);
	void	SetConstants(double newconst, int dim); // create a vector, the last one *= sqrt(n-1)
	void	ComputeCminus();
	void	ComputeConstant();
	void	SetBoundaries(dVec& bounds){m_bounds=bounds;}
//	void	SetSimplexBound(double bound){m_SimplexBound=bound;}
	void	UpdateBasisVectors();
	int		IsPenalty();

	void	ProcessAllDyn(support_vector* SV);
	real	ProcessAllNew(dVec& X);   // new query for the value at X
//	void	ProcessAllDynAggressive(support_vector* SV);
	int IsAggressive() {return Match;}



	int		CheckExistingSV(dVec&  X);
	int		EquivalentSV(dVec& X, dVec& Y);

	int PopNextCandidate(dVec& x, double& d, dVec& RHS);

#ifdef TRIANGULATION
	void	ComputeTriangulation();  // wrappers
	int		GetNumberOfSimplices();
	GBPolytope * GetPolytope(int i) {return HeapPossibleMin.GetPolytope(i); }  
#endif

#ifdef RANDOM_GENERATOR
	dMat	Savedbasis;
#endif
};




#endif
