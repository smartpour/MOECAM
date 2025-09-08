#if !defined ECAM_H
#define ECAM_H


// requires: prob.h heap.h forest.h memblock.h

// there is a lot of auxiliary staff here (various classes, etc.) 

//#define MPI
#define TNT_NO_BOUNDS_CHECK
#define real double // 8 bytes: choose one of these


#define SHORT_HEAP
#define USEHEAP

typedef unsigned int SVINDEX ;

//Do not use box constraints, the function is periodic
//#define PERIODIC

// if the fucntion is not periodic, it can be extended as an even periodic function
// this also avoids using constraints
#ifndef PERIODIC
#define ANTIPERIODIC
#endif

#define BOXCONSTRAINTS

#ifndef BOXCONSTRAINTS
#define SIMPLEXCONSTRAINTS
#endif

//#define MINLP  not used
/* ------------------MINLP-------------------------*/
#ifdef MINLP

#ifdef BOXCONSTRAINTS
#undef BOXCONSTRAINTS
#endif

#ifndef GENERALCONSTRAINTS
#define GENERALCONSTRAINTS
#endif

#endif
/* ------------------MINLP-------------------------*/

#define My_precision 0.02
#define limitworst 1e10

#include <cstdio>
#include <cstdlib>
#include <list>
#include <ctime>
#include <deque>
#include <set>
#include <vector>
//#include <algorithm>

using namespace std;

#include "../tnt/tnt.h"
#include "../tnt/jama_lu.h"
//using namespace TNT;
//using namespace JAMA;

#include "memblock.h"
#include "mymath.h"

//#include "memblockmpi.h"
//#include "probs.h"

#include "heap.h"
#include <math.h>


typedef TNT::Vector<float>		fVec;
typedef TNT::Vector<real>		dVec;
typedef TNT::Vector<int>		iVec;
typedef TNT::Matrix<real>		dMat;
typedef TNT::Matrix<int>		iMat;
typedef TNT::Vector<SVINDEX>	siVec; // short int (2 byte)

const real Infinity = 1.0e16;
const real Special1 = 1.0e20; // inverse special

#define sqr_(a) ((a)*(a))
#define min_(a,b) ((a)<(b)?(a):(b))
#define max_(a,b) ((a)>(b)?(a):(b))



#define EXP_LIMIT 20
// the epsilon separating from the boundary of the simplex
const real Eps_2 = 2.1e-3;
const real EXP_EPS=5; //ln(0.002)=ln(Eps_2)


// defaults
// max number of iterations of CA and local methods
const	int		MaxIter=200;
// how many points to use for local search at the end
const int LOCAL_SEARCHES=2;
// important constants: f = min(f,Cut) + Const
const double _Cut=5;
const double _Const = 10.0;

double ElapsedTime();
void   ResetTime() ;

#define sqr_(a) ((a)*(a))
#define min_(a,b) ((a)<(b)?(a):(b))

void	FormTransformMatrices(int Dim, dVec& a, dVec& b);
void	FormTransformMatricesScale(int Dim, dVec& a, dVec& b);
dVec	TransformCube2SimplexScale(dVec& x);
dVec	TransformSimplex2CubeScale(dVec& x);
void	TransformSystemOfConstraintsScale(dMat& Constr, dVec& RHS, dMat& ModConstr, dVec& ModRHS);

void	TransformSystemOfConstraintsIdentity(dMat& Constr, dVec& RHS, dMat& ModConstr, dVec& ModRHS);
dVec	TransformSimplex2SimplexInverse(dVec& x);
dVec	TransformSimplex2SimplexDirect(dVec& x);
void	FormTransformMatricesIdentity(int Dim, double beta);


#pragma optimize( "", off )
//-------------------------
class SawToothCover;
class support_vector;

class support_vector {
public:
	unsigned int label;
	dVec		 vec;
	real		 funvalue;

	void	SVForm(dVec& x, real val);
// inc funvalue by meps to break the ties.
	void	Increment();
	short int	ChangeF(real newval); // returnts old value < newval
	// returns the coordinates of the point
	void ReturnX(dVec& x);
	support_vector* This();
	void PrintMe(FILE *f);
};

// extended support vector
class support_vector_X:public support_vector {
public:
	real	Dval;
	bool operator<(const support_vector_X & x) const { return funvalue < x.funvalue;};
};


typedef deque <support_vector>	SVDeque;
typedef list <support_vector_X> SVList;
typedef list <support_vector_X>::iterator SVListIter;
typedef deque <support_vector_X>	SVXDeque;

typedef struct {
	unsigned int PossibleKiller;
} PossibleKillerStruc;
typedef deque <PossibleKillerStruc>  PossibleKillerList;

//#define EUINT unsigned int
#define SVSetNodePtr EUINT

#ifdef _MSC_VER  
// if the compiler does not recognise this type, change it to another int type 8 bytes long
// like long long int
typedef  __int64 ULINT; //this type myst be 8 bytes long
#else
typedef unsigned long long int ULINT; //this type myst be 8 bytes long
#endif


// just to pack these into 8 bytes
#define vecnumber SVSetNodeData[0]
#define children__ SVSetNodeData[1]
#define numchildren__ SVSetNodeData[2]

#ifndef SMALL_DIM
#define parent__ SVSetNodeData[3]
#else
#define parent__ SVSetNodeData[2]
#endif


// todo
// make distributed priority queue
// otherwise limited by 10^8 nodes about 800 MB ram
// shoudl work well

// when adding a killer no need to update remote nodes unless it is the first one
#define mynode SVSetNodeData[2]  // shares with children - 
//interpretation depends on numchildren, if 0 then it's leaf, so in heap


#define VEC_VEC(a) ((a) & 0x00FFFFFF)
#define POSVEC_VEC(a,b) ((a)<<24 | (b))
 
#define IS_ROOT(a) (((a) &0xFFFFFF) == 0xFFFFFF) 
#define SET_ROOT(a) ((a) |= 0xFFFFFF) 

// pos 3F would mean root (1,2,3,...,n), no more than 255 variables
#ifdef SMALL_DIM
#define POS_NUMCHLF(a) (((a)>>28))   // + 1 as 0 children does not make sense, so 0 means 1, and to ensure F is NOT used (special for root)
#define POS_VEC(a) (((a)>>24) & 0x0F)   // only 15 possible positions
#else
#define POS_VEC(a) (((a)>>24) & 0xFF)   // only the last 6 bits, the first 2 bits reserved
#endif


/* ---------- Aux classes------------------------*/
#define ChildrenBlockSize 4

class SVSetNode;
class Children;
EUINT CreateChildrenArray(EUINT size);
void DeleteChildrenArray(EUINT& first);
Children* GetMemBlockAt(EUINT b);

class Children {
public:
	EUINT m_children[ChildrenBlockSize]; // an array of pointers (of size Dim)
	EUINT NextBlock;
	EUINT& operator[](int i) {
		if(i<ChildrenBlockSize) return m_children[i];

		Children* block=this;
		div_t b=div(i, ChildrenBlockSize);
		int j;
		for(j=0;j<b.quot;j++) 
			if(block->NextBlock!= MB_BADINDEX) block=GetMemBlockAt(block->NextBlock);
			//else return MB_BADINDEX; // error!!!
		return (block->m_children[b.rem]);
	};
	Children() { NextBlock=MB_BADINDEX; };
} ;

class ChildrenList {
public:
	EUINT First; // element
	EUINT Last;  // element
	EUINT FirstChild;
	EUINT LastChild;

	ChildrenList();
	~ChildrenList();
	void push_back(EUINT v);
	EUINT pop_back();
	EUINT pop_front();
	EUINT AddBlock();
	inline EUINT& operator[](int i); 
	void SetAt(int i, EUINT val);
	EUINT erase(int i);
	int size() {return Last-First;};
	int empty() {return size()==0;};
	void copy(ChildrenList* source, int start);
};

/* ---------- Aux classes------------------------*/


/*-----------------------------------------
  Local minumum of aux. function. Contains the value in the minimum, the value of the function
  at this point (if evaluated preemptively) and the inddex of the SV, which differs it from parent.
  This is in the SVSetPosNumber, which is an int (4 bytes), the 1 byte is the position,
  the last 3 bytes is the SV number. Splitting and joining these variables is
  performed with the above macros.

  The killers is the pointer to the list of possible killers (for virtual branches).
  If this list is empty, it is deletes and killers=NULL.

  For roots, the first n elements are the actual index, and what follows are the
  possible killers

--------------------------------------*/
/*-----------------------------------------
	This class incorporates 2 data structures, a tree (rather forest) and
	Fibonacci heap, containing leaves. Leaves are the local minima of saw-tooth
	cover. There are two essential members: Heads (list of pointers to roots) and
	HeapP, a pointer to heap (created at initialisation (Init), destroyed in the 
	"destructor" EraseAll). perhaps couled be moved into proper contructor/destructor,
	although dimension of the problem Parent->Dim should be known.

	There are 2 types of methods, the routine insert/delete and problem-specific

	Forest is needed for minimise. Only one tree for serial and forest for parallel.
	There are also virtual branches. These are identified with a list of possible killers.

	Works like this:
	1 take top node
	2 is it virtual? if yes, split into children, and insert into the heap. goto 2
	3 it is the leaf. Then take this node and form SV
	4 process this node with the generated SV
	5 process the rest of the nodes

  The root (if there is more than 1 root), keeps its index vector in full, as the first
  N members of the list Killers. When processing virtual nodes, the first n elements are 
  skipped for the roots. For leaves children=NULL, for branches it becomes an array of pointers.

----------------------------------------*/
class SVSetNode {
public:

// it is possible to keep the maxvalue (averageF) instead of calculating it every time, at the expense of extra 4 bytes
// investigate this possibility.

#ifdef SMALL_DIM
	EUINT	SVSetNodeData[3]; // packs children and parent
#else
	EUINT	SVSetNodeData[4]; // packs children and parent
#endif
	ChildrenList*			killers; // for virtual branches   // perhaps cast them to children_ ?
	real Dval;				 // value of the minimum

// end data members--------------------


	SVSetNode(); // constructor, assigns NULL to pointers
	~SVSetNode(); 
	void Delete();
	void Init();
	SVSetNode* This() {return this;}
	int IsValid();
	SVSetNodePtr GetParent() ;

#ifdef SMALL_DIM
	int	GetNumChildren() { 
		EUINT a=POS_NUMCHLF(vecnumber);
		if(int((a) - 15) <= 0) return -1; 
		else return a; };
	void SetNumChildren(int ncld) { 
		if(ncld==-1) vecnumber|=0xF0000000; else {
			vecnumber &= 0x0FFFFFFF; vecnumber |= (ncld) << 28; }};
#else
	inline int	GetNumChildren() {if((numchildren__ & 0x00FFFFFF ) != 0x00FFFFFF) return (numchildren__ & 0x00FFFFFF); else return -1; };
	inline void SetNumChildren(int ncld) {if(ncld<0) numchildren__=0x00FFFFFF; else numchildren__ = (numchildren__& 0x80000000) + ncld; };
#endif
	inline int  GetHeapNode() {return (mynode>>31);}
	inline void SetHeapNode(int n) {if(n>0) mynode |= 0x80000000; else mynode &=0x7FFFFFFF; }

	// attaches a child "node" to this, at position pos
	void AddChild(SVSetNodePtr thisnode, SVSetNodePtr node, int pos);

	// deletes all children. Used to clear memory when destroying the tree
	void Clear();

	// removes just the reference
	void RemoveChild(SVSetNodePtr child, int pos);//	{	children[pos]=NULL; }

	// these two methods test cond (2) for SVector v
	// the first version is to test nodes other than root (index is not important)
	// the second version is to test root, in which case index should be the
	// list of indices of SV comprising this node
	// returns 0 if passed, 1 if failed (dominance), 2 if nonstrict dominance, and 3 if below best function value,
	// in which case branch is eliminated
	int TestVector(dVec& v, siVec* index);
	int TestVectorIndex(dVec& v, siVec* index);
	int TestVectorIndexQ(dVec& v, siVec* index);
	int TestVectorQ(dVec& v, siVec* index);
	int	TestVectorPenalty(dVec& v, siVec* index);

	// computes indices of SVectors by referring to the parents successively
	// ind is the reference to working memory, can contain any data.
	// returns vec, the vector of int
	void ComputeVectors(siVec& vec, siVec* ind, SVSetNodePtr thisnode);

	// for the root node generates the indices of SVectors. for ROOT returns 1,2,3,,,.n
	// otherwise returns the acural indices, stored in VectorPos
	void GenerateInitVector(siVec* initvec, SVSetNodePtr thisnode);

	// tests cond (1) with SV at position pos. Assumes that the parent
	// satisfies this condition, and hence tests only column pos
	// index contains the actual SV indices. Also returns the olddiag, the value
	// of the element on diagonal to be replaced. It will be used in updating DVal
	int TryNewVectorIndex(support_vector* SV, int pos, siVec* index, real &olddiag);

	// calculates the aux. minimum as the trace of L. Updates the old value by 
	// replacing only one diag. element
	void CalculateDval(SVSetNodePtr parentnode, SVSetNode* ParentNode, support_vector* SV, real olddiag);
	void CalculateDvalComplete(SVSetNodePtr parentnode,  SVSetNode* ParentNode, support_vector* SV, real olddiag, siVec& idx);
	void	CalculateDvalandMin(siVec& index,real& d, dVec& minpos);

	bool operator<(const SVSetNode& x) const { return Dval < x.Dval;}; // uses actual Dval=H(x)

	void CopyTo(SVSetNode* copy);

	void PrintMe(FILE* f);

	void	CalculateDvalModified(siVec& index,int i,real fmax, support_vector& SVT);
	real	CalculateDvalExplicit(siVec& index, support_vector& SVT);
	real	ComputeAverageF(siVec* index);
	real	ComputeFunValue(dVec& X, siVec* index);

};

struct HeadStruc {
	SVSetNodePtr Head;
	siVec*	p_index;
};

#include "forest.h"


// Class Saw tooth cover-------------------------------here starts the visible interface-----

// this is our main object
/*-------------------------------------------------------------------------------------------
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

------------------------------------------------------------------------------------------*/

class Ganso;
// most members should not be accessed from outside, they are public for now (development stage)
class	SawToothCover {
public:
// where to put the results
//	FILE *fp;
	FILE *logfile;


    SVList  PreparedPointsList;

	SawToothCover();
	~SawToothCover();

//	void InitProblem(Problem * _prob) {Prob=_prob;}
	Ganso*				MyParent;

	SVDeque				SVectors;	// here support vectors live
	SVDeque::iterator	iter;
	SVDeque::iterator	t_iter; // just iterator


// SVSets live here
	Forest				HeapPossibleMin;

// aux staff
	siVec				m_lastindex, m_TempSiVec;
	real				tempFval;
	int					Match;
	int					m_lasterror;

// solution sits here
	dVec	xglob_x, temp_x, savedtemp_x;
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

	double Cut;
	double Const;


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
	void	FormNewSV();


	void InitialSteps(int d, int iters);
	int MakeStep(dVec&  X);
	int MakeStepbefore();

	int		StopCriteria();
//	void	Minimise(int code=0);  //main routine
	void	Minimise(int d, int iters);
	int		LastError();

	void PrepareFinalList(double* List, int* size);


// the rest are auxiliary routines	
	int PopNextCandidate(dVec& x);

	int PopNextCandidateExternal(dVec& x);
	void	ReturnSol(dVec &x);


// the aux. methods

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
	void	SetBoundaries(int n, double* bounds){m_bounds.newsize(n); for(int i=0;i<n;i++) m_bounds[i]=bounds[i];}
//	void	SetSimplexBound(double bound){m_SimplexBound=bound;}
	void	UpdateBasisVectors();
	int		IsPenalty();

	int IsAggressive() {return Match;}

	int		CheckExistingSV(dVec&  X);
	int		EquivalentSV(dVec& X, dVec& Y);

	int PopNextCandidate(dVec& x, double& d, dVec& RHS);

	void ExtendPeriodic(dVec& x);
	void ExtendAntiPeriodic(dVec& x);


};
#pragma optimize( "", on )

#endif

