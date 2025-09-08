/*
	06.07.03. Gleb Beliakov

	This file contains two class declarations, SVSetNode and Forest

  SVSetNode is the class representing SV combination.
  SVSetNode has additional data members which allow keeping tree structure, 
  namely pointers to the parent and array of pointers to children, and
  also reference to heap node.

  The methods are routine for trees, plus methods to process SVector (cond(2))
  and to reconstruct the index vector from the family tree.


  children is a pointer to array of pointers. If this is a leaf, children is NULL.
  if not, an array of pointers is created and referred to by children. This saves memory

  perhaps merge mynode and children ? (how will I know if this is the leaf?)

*/


#define MAXDIMCONSTR 20 // value + dim
#define HEAPLIMITEXTRA 3

class ConstrMin {
public:
 float ConstMinType[MAXDIMCONSTR];
 real GetVal() {return ConstMinType[0];}
 void  SetVal(real r) { ConstMinType[0]=(float)r;}
 void  SetVec(dVec & p) {	for(int i=0;i<p.size(); i++) ConstMinType[i+1]=(float)p[i];}
 void  GetVec(dVec & p) {	for(int i=0;i<p.size(); i++) p[i]=ConstMinType[i+1];}
};

#define ConstMinPtr EUINT

#define MAXDIMINDEX 20 //  dim
class SavedIndex {
public:
 SVINDEX Index[MAXDIMINDEX];
 void  SetVec(siVec & p) {	for(int i=0;i<p.size(); i++) Index[i]=p[i];}
 void  GetVec(siVec & p) {	for(int i=0;i<p.size(); i++) p[i]=Index[i];}
};
#define IndexNodePtr EUINT


#define ChildrenBlockSize 4

typedef vector<char> shortindexvector;
typedef vector<char>::iterator shortindexvectoriter;

typedef set<EUINT> indexset;
typedef set<EUINT>::iterator indexsetiter;


class Forest {
public:
	int size, sizevirtual, sizemem; // aux. members for testing
	int sizepacked;
	int size_constrained,size_infeasible;

	siVec m_initvec, m_index, temp_index; // just not to create it in all functions
	// provide temp. storage passed to through pointer

	support_vector	SVT;

	deque<HeadStruc> Heads; // here we keep the roots of the trees

	bheap_t* HeapP;	

	SVSetNodePtr m_TempChildren;


	int		m_heaplimit, m_heapoverlimit;
	float	m_heapworst;
	float   GlobalWorstHeap;

public:

	void	SetHeapLimit(int n){m_heaplimit=n; m_heapoverlimit=n*HEAPLIMITEXTRA;}

	void	ComputeHeapWorst(int n, int extra);

	float	GetGlobalWorstHeap();

// constructor
	void Init(); // to create Heap and aux. storage
// destructor
	void	EraseAll();

// Routine methods
	int GetSizeMem()	{return sizemem; }; // returns the size of the forest
	int GetSize()	{return size; }; // returns the size of the forest
	int SizeRoot()	{return Heads.size();}; 
	int Size()		{return (SizeRoot() <<24) + HeapP->nodeCount;}; // the size of the heap (leaves only)

	int ComputeSize(SVSetNodePtr node);
	int ComputeSize(SVSetNodePtr node, int not_this_child);
	int	ComputeSizeRecursive(int req, SVSetNodePtr node, SVSetNodePtr& foundnode);


	// returns 1 if heap is empty
	int		Empty() {return ((HeapP==NULL) || (HeapP->nodeCount==0));};

	double	GetConst();

	SVSetNodePtr Front();	// returns the top of the heap, does not remove 
	void	RemoveTop(); // removes the top
	void	RemoveNodeFromHeap(SVSetNodePtr node); // only removes node from the heap, 
	void	RemoveNodeFromHeap(SVSetNodePtr node, SVSetNode* Node); // only removes node from the heap, 
										//	not from the tree
	void	AddRootNode(  SVSetNodePtr node); // starter: called in InitPopulate
	void	AddTree(SVSetNodePtr root); // add a branch
	siVec*	GetVecAddress() {return &m_initvec;}; // provides working memory

	void	AddLeaf(SVSetNodePtr node); // called recursively to find the leafs and insert into the heap
private:	
	void	ClearBranch(SVSetNodePtr branch); //like EraseBranch, but not removed from heap

public:
	void	EraseBranch(SVSetNodePtr branch, int processparent=-1); 
	void	EraseRootEntry(SVSetNodePtr branch); // like erase branch, but processes roots

// these are problem-specifis methods
private:
// called internally from ProcessAll. This is the working horse
// given new SV, and the root index vector *initvec (calculated before the first call)
// returns 2 if node is too big (needs to be erased), 1 if test (2) fails (needs to split),
// children processed if any, otherwise split.
// and 0 if not affected by SV. In this case processing stops (no children can fail (1))
	int		ProcessNode(SVSetNodePtr node, support_vector* SV, siVec* initvec);

// this one is used if virtual branches are allowed. Takes the top of the heap
// and processes it. Called from ProcessFront
// Simply ensures that node does not contain virtual branches. Executes the split and
// updates the tree and heap if there are. Retirns 0 if finished. initvec is working memory only
	int		ProcessTopNode(SVSetNodePtr node,SVSetNode* Node, siVec* initvec);

public:
// called from Minimise. Starts at roots and processes all trees in the forest. Splits
// and updates the tree automatically
	void	ProcessAll(support_vector* SV);


// enshures that the node is leaf and does not contain virtual branches. Calls ProcessTopNode
// to execute the split
	int		ProcessFront();
	int		ProcessFront(SVSetNodePtr &head,	SVSetNode** Head);

	void	ProcessNodeConstraintsSimple(SVSetNodePtr node, SVSetNode* Node, siVec& index, dVec& soln);


public:
// this piece is for queries only ------------------
	list<SVSetNodePtr> m_ListAffectedNodes;
	list<IndexNodePtr> m_ListAffectedNodesIdx;
	indexset m_indexset;

	real	delta;
// as above, but no tree modifications, just augment the list
	int		ProcessNodeQ(SVSetNodePtr node, support_vector* SV, siVec* initvec);
	void	ProcessAllQ(support_vector* SV);

};


