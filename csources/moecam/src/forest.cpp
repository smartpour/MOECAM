
#include "ecam.h"


extern SawToothCover*	Parent;   // points to the parent STCover class


#ifdef MPI
extern MemoryBlockMPI<SVSetNode>	MBSV;
#else
extern MemoryBlock<SVSetNode>	MBSV;
#endif

extern MemoryBlockE<Children>    MBCL; // blocks for SVSetNode and children array
extern MemoryBlockE<ConstrMin>   MBCM; // blocks for ConstMin array



int GlobPos;
//SVSetNode GlobNode;
//float	GlobalWorstHeap;

#define IS_VALID_PTR(a) ((a)!=0 && (a) !=0xFFFFFFFF)

// these are required for the heap ------------------------------------
#define LARGEHEAP
#ifdef LARGEHEAP


#pragma optimize( "", off )
#pragma OPTIMIZE OFF

//KEY_TYPE ft; // global
#define PACKFLOAT(f) ( (f  & 0xFFFFFFFF))
#define f2ulint(a)  (*(EUINT*) &(a))
#define ulint2f(a) (*(float*) &(a))



KEY_TYPE _getKey(node_t theNode)
{
	EUINT UL1; KEY_TYPE ft;
	UL1=theNode.data>>32;
	ft=ulint2f(UL1) ;
	// use the first 4 bytes
	return -ft;
	return -fabs(ft);
}

KEY_TYPE _getKey(DATA_TYPE Node)
{
	EUINT UL1; KEY_TYPE ft;
	UL1=Node>>32;
	ft=ulint2f(UL1) ;
	return -fabs(ft);
}


INDEX_TYPE _getIndex(node_t theNode) {	return 0; }

void _setIndex(node_t theNode, INDEX_TYPE I) {}

#else

inline KEY_TYPE _getKey(node_t theNode)
{
	SVSetNode* n=MBSV.GetAt(theNode.data, &GlobNode);
	return -fabs(n->Dval); // just to process negative (which is a flag for special case), "-" because we need min, heap uses descending sort
}

inline INDEX_TYPE _getIndex(node_t theNode)
{
	EUINT u;
	SVSetNode* n=MBSV.GetAt(theNode.data,&GlobNode);
	u= n->mynode;
	return u;
}
inline  void _setIndex(node_t theNode, INDEX_TYPE I)
{
	SVSetNode* n=MBSV.GetAt(theNode.data,&GlobNode);
//	if(n->mynode!=BADINDEX || n->vecnumber==MB_BADINDEX)
	n->mynode = I;
	MBSV.SetAt(theNode.data,n);
}
#endif

INDEX_TYPE Merge(bheap_t *h, SVSetNodePtr el) 
{	
	return bh_insert(h, el); 
}


/*
#ifdef _MSC_VER 
#define LARGECONST 0xFFFFFFFF00000000UL 
#else
#define LARGECONST 0xFFFFFFFF00000000ULL 
#endif
*/

INDEX_TYPE Merge1 (bheap_t *h, SVSetNodePtr el, float f) 
{	
	ULINT UL;
	UL=f2ulint(f);//<<32;
#ifdef _MSC_VER
       UL = (UL << 32) & 0xFFFFFFFF00000000UL ;
#else
       UL = (UL << 32) & 0xFFFFFFFF00000000ULL ;
#endif
//	UL = (UL << 32) & LARGECONST ;
	UL = UL + el;
	return bh_insert(h, UL); 
}

// ensure no children
#define InsertNodeIntoHeapOld(n,N) {/*(N)->mynode=*/Merge(HeapP,(n));  }

#define InsertNodeIntoHeap(n,N) { \
(N)->SetHeapNode(0); Merge1(HeapP,(n),(float)N->Dval);  \
GlobalWorstHeap=MMmax(GlobalWorstHeap,(float)fabs(N->Dval)); }


#pragma optimize( "", on ) 


// these are required for the heap ------------------------------------


SVSetNode::SVSetNode()   { 	Init(); }

void SVSetNode::Init() { 
		children__=MB_BADINDEX;
		vecnumber=MB_BADINDEX; 

		parent__=MB_BADINDEX;
		Dval=Infinity; 
		killers=NULL;
}

SVSetNode::~SVSetNode() { }
void SVSetNode::Delete() {
		if(killers !=NULL) {
//			if(Dval<0) {MBCM.FreeBlockC((EUINT)killers);} //delete (support_vector*)killers; // special case , casting
//				else delete killers; // clears itself
			killers=NULL;
		}
 }

void SVSetNode::CopyTo(SVSetNode* copy) { 
	copy->Dval=Dval;
	copy->vecnumber= vecnumber;
}

void SVSetNode::Clear()
{
	int i;
	if(!IsValid()) return;
	if(children__ != MB_BADINDEX)
	for(i=0;i< GetNumChildren();i++) 
	{
		MBSV.FreeBlock(children__ + i);  
	}
	children__=MB_BADINDEX; // do we care?

	if(killers !=NULL) {
//			if(Dval<0) {MBCM.FreeBlockC((EUINT)killers);}  // special case 
//			else 	delete killers;
			killers=NULL;
	}
}


void SVSetNode::RemoveChild(SVSetNodePtr child, int pos)	
{		
	SVSetNode GlobNode;
		SVSetNode* n=MBSV.GetAt(children__ + pos, &GlobNode);
		n->Clear();
		MBSV.SetAt(children__ + pos,n);
}

int SVSetNode::IsValid()
{
	if(POS_VEC(vecnumber) != 0xFF || parent__==MB_BADINDEX) return 1; else return 0;
	return 1;
}

void SVSetNode::PrintMe(FILE* f)
{
	if(f==NULL) return;

	short int pos = POS_VEC(vecnumber);
	int vecnumber_ = VEC_VEC(vecnumber);

	fprintf(f,"Node %f , vecpos %d %d\n", Dval, pos, vecnumber_);
}

inline SVSetNodePtr SVSetNode::GetParent() { 	return parent__; }

// compute vectors given the chain in the tree
void SVSetNode::ComputeVectors(siVec& vec, siVec* ind, SVSetNodePtr thisnode)
{		
		siVec* filled =&(Parent->m_TempSiVec);
		int allpos=Parent->Dim;

		(*filled)=0;
		SVSetNode* SNode;
		SVSetNodePtr SNodePtr;
		int i;
		short int pos = POS_VEC(vecnumber);
		int vecnumber_ = VEC_VEC(vecnumber);
		if(IS_ROOT(vecnumber)) {
			GenerateInitVector(ind, thisnode);
			if(ind==NULL)
				{for(i=0;i<vec.size();i++) vec[i]= i; return;}
			else
				{for(i=0;i<vec.size();i++) vec[i]=(*ind)[i]; return;}
		}

// else go up one level recursively
		vec[pos]=vecnumber_; 
		(*filled)[pos]=1;
		allpos--; // once filled all, stop

		SNodePtr=GetParent();
// this is recursion (implicit)
		while(SNodePtr!=MB_BADINDEX) {
			SVSetNode GlobNode;
			SNode = MBSV.GetAt(SNodePtr, &GlobNode);
			pos = POS_VEC(SNode->vecnumber);
			vecnumber_ = VEC_VEC(SNode->vecnumber);

		  if(IS_ROOT(SNode->vecnumber)) {// root
			SNode->GenerateInitVector(ind, SNodePtr);
			if(ind==NULL)
			{	for(i=0;i<vec.size();i++) if(!(*filled)[i]) vec[i]=i; 
				return;
			} else
			{	for(i=0;i<vec.size();i++) if(!(*filled)[i]) vec[i]=(*ind)[i]; 
				return;
			}
		  } // pos ==FF

			if(!(*filled)[pos]) 
			{	vec[pos]=vecnumber_;
				(*filled)[pos]=1;
				allpos--;
				if(allpos <=0) return; // all positions filled
			}

			SNodePtr=SNode->GetParent(); 
		}
}

int SVSetNode::TestVectorIndex(dVec& v, siVec* index)
{
	int i,flag=0;

	real u; // make them global?
	if(index==NULL) 
	for(i=0; i<v.size(); i++) {
		u=Parent->SVectors[ i ].vec[i] - v[i];
		if( u > 0) {return 0;}
		if(u==0) 
		{ flag +=1;
		} 
	}
	else
	for(i=0; i<v.size(); i++) {
		u=Parent->SVectors[ (*index)[i] ].vec[i] - v[i];
		if( u > 0) {return 0;}
		if(u==0) { 
			flag+=1; // indicates nonstrict dominance
		}	
	}
	if(flag==1) return 2;  // nonstrict dominance, test not failed but an extra minimum (copy) is needed
	else 
		return 1;  // strict dominance, test failed
}

int SVSetNode::TestVector(dVec& v, siVec* index)
{
		// remove these SVSets: they are above the best function value
		if(Dval != Special1 /*&& Dval>0*/) // but only those w/o constraints
			if(fabs(Dval) - Parent->GetConst()  > Parent->frecord ) return 3; 

		real u;
		short int pos = POS_VEC(vecnumber);
		int vecnumber_ = VEC_VEC(vecnumber);

		if(IS_ROOT(vecnumber)) { // root, use index explicitly
			return TestVectorIndex(v,index);
		} else { // test only one diagonal element, assumes parent has passed the test
			u = Parent->SVectors[ vecnumber_ ].vec[pos] -  v[pos];
			if( u >  0) { return 0;}  ////?????????????? >=
			if( u == 0) 
			{	return 2; }      // nonstrict dominance
			return 1; // strict dominance
		}
		return 1; //should not get here
}

// as above but strict inequality
int SVSetNode::TestVectorIndexQ(dVec& v, siVec* index)
{	int i;
	for(i=0; i<v.size(); i++) {
		if( Parent->SVectors[ (*index)[i] ].vec[i] > v[i]) {return 0;}
	}
	return 1;
}

int SVSetNode::TestVectorQ(dVec& v, siVec* index)
{
		short int pos = POS_VEC(vecnumber);
		int vecnumber_ = VEC_VEC(vecnumber);

		if(IS_ROOT(vecnumber)) { // root, use index explicitly
			return TestVectorIndex(v,index);
		} else { // test only one diagonal element, assumes parent has passed the test
			if( Parent->SVectors[ vecnumber_ ].vec[pos] > v[pos]) {return 0;} //<=
			return 1;
		}
		return 1; //should not get here
}

//returns 0 if insuccessful, 1 if cond 1 is satisfied
int  SVSetNode::TryNewVectorIndex(support_vector* SV, int pos, siVec* index, real &olddiag )
{
	// test of cond (1). Cond (2) supposedly passed. Only one column is tested. 
	// olddiag is returned for another method (to calculate Dval)
	olddiag= Parent->SVectors[ (*index)[pos]].vec[pos]; // old diagonal
	real val=SV->vec[pos];
	for(int i=0;i<Parent->Dim;i++) { // if diag not smaller than this column, exit 0
		if((pos!=i) && (val >=  Parent->SVectors[ (*index)[i]].vec[pos]) ) return 0;   //??????????????????///// <=
	}
	return 1;
}

void SVSetNode::GenerateInitVector(siVec* initvec, SVSetNodePtr thisnode)
{
	// this can be called only for one of the roots, to generate its
	// index set, for subsequent processing of children

	if(initvec==NULL) return;
// will retrieve init vector from the killers list
	if(IS_ROOT(vecnumber) ) // means root
	{		
// attempt to find this root in the list of heads
	deque<HeadStruc>::iterator iter;

	for(iter=Parent->HeapPossibleMin.Heads.begin(); iter!=Parent->HeapPossibleMin.Heads.end(); iter++)
		if((*iter).Head == thisnode) {
			(*initvec) = *((*iter).p_index);
			return;
		}
		 for(int i=0;i<Parent->Dim; i++) (*initvec)[i]=i;
	}
}

// various calculations of Dval-------------------------------------------------
void SVSetNode::CalculateDvalComplete(SVSetNodePtr parentnode,  SVSetNode* ParentNode, support_vector* SV, real olddiag, siVec& idx)
{// assume index is computed! This is for DVal of parent corrupted by SpecialValue
	real f=0;
	for(int i=0;i<Parent->Dim;i++) f+= Parent->SVectors[ idx[i]].vec[i] ; // trace of parent
	f=(f+1.0)/Parent->C_minus; // this is parent's value

	Dval = f + (- (olddiag) + SV->vec[ POS_VEC(vecnumber) ])/Parent->C_minus;
	Dval += Parent->GetConst();
}

// this is normal computation in unconstrained case, just one operation
void SVSetNode::CalculateDval(SVSetNodePtr parentnode, SVSetNode* ParentNode, support_vector* SV, real olddiag)
{
//	ConstrMin* cm;

	if(ParentNode->Dval>0)
		Dval=ParentNode->Dval + (- (olddiag) + SV->vec[ POS_VEC(vecnumber) ])/ Parent->C_minus;
	else { // was saved there
//		cm=MBCM.GetAt((EUINT) (ParentNode->killers));
//		Dval=cm->GetVal()  + (- (olddiag) + SV->vec[ POS_VEC(vecnumber) ])/Parent->C_minus;
	}
}


void	SVSetNode::CalculateDvalModified(siVec& index, int i,real fmax, support_vector& SVT)
{
	// computes local min with one of the SV (i-th) modified, so that its diagonal is fmax
	// result (coords of the minimum) returned in SVT and SVT.funvalue (true value)
	int j;
	real d=0;
	d=fmax;
	SVT.vec[i]=fmax;

	for(j=0;j<Parent->Dim;j++) {
		if(i!=j) {d += Parent->SVectors[ index[j]].vec[j]; SVT.vec[j] = Parent->SVectors[ index[j]].vec[j];}
	}
	// trace computed

	d=(d+1.0)/Parent->C_minus; 

	// compute the coordinates
	for(j=0;j<Parent->Dim;j++) {
		SVT.vec[j] = 	d/Parent->m_Constants[j] - SVT.vec[j] ; 
	}
	SVT.funvalue = d;
}

real	SVSetNode::CalculateDvalExplicit(siVec& index, support_vector& SVT)
{
	// compute as max min at the point SVT.vec
	real r=0;
	for(int i=0;i<Parent->Dim;i++)
		if(r < (SVT.vec[i] + Parent->SVectors[ index[i]].vec[i])*Parent->m_Constants[i]  ) 
			r = (SVT.vec[i] + Parent->SVectors[ index[i]].vec[i])*Parent->m_Constants[i] ;

	return r;
}

void	SVSetNode::CalculateDvalandMin(siVec& index,real& d, dVec& minpos)
{
	d=0;
	int j;
	for(j=0;j<Parent->Dim;j++) {
		d += Parent->SVectors[ index[j]].vec[j]; 
		minpos[j] = Parent->SVectors[ index[j]].vec[j];
	}
	// trace computed

	d=(d+1.0)/Parent->C_minus; 
	// compute the coordinates
	for(j=0;j<Parent->Dim;j++) {
		minpos[j] = 	d/Parent->m_Constants[j] - minpos[j] ; 
	}

}


// various calculations of Dval-------------------------------------------------


///////////////////////////////Forest////////////////////////

void Forest::Init() 
{
	HeapP = bh_alloc(); 
	// take one constrMin (otherwise it will start at 0 and it can't delete it: memory leak)
	MBCM.GetNextFree(); // forget this one for now, later track it and delete

	size=sizemem=0; sizevirtual=0;
	size_constrained=size_infeasible=0;

	m_initvec.newsize(Parent->Dim); 
	temp_index.newsize(Parent->Dim); 
	m_index.newsize(Parent->Dim);
	SVT.vec.newsize(Parent->Dim);
	m_ListAffectedNodesIdx.clear();
	m_ListAffectedNodes.clear();

	m_TempChildren=MBSV.GetNextFree(Parent->Dim);

	GlobalWorstHeap=-10e-10;
	m_heapworst=limitworst;
};

float	Forest::GetGlobalWorstHeap() {return GlobalWorstHeap;}

int Forest::ComputeSizeRecursive(int req, SVSetNodePtr node, SVSetNodePtr& foundnode)
{
	int sz=0;
	foundnode=MB_BADINDEX;

	int i;
	SVSetNode Node;
	SVSetNode* n=MBSV.GetAt(node, &Node);  
	if(n->IsValid()) {
	if( n->children__!=MB_BADINDEX && n->GetNumChildren()>0) {
		for(i=0; i<n->GetNumChildren(); i++) {
			sz += ComputeSizeRecursive(req, n->children__ + i, foundnode);
			if(foundnode!=MB_BADINDEX) return sz;
		}
	} else sz=1; 
	} else sz=0;

	if(sz>0) sz++;

	if(sz>=req) foundnode=node;
	return sz;
}

int Forest::ComputeSize(SVSetNodePtr node)
{
	int sz=1;
	int i;
	SVSetNode Node;
	SVSetNode* n=MBSV.GetAt(node, &Node);
	if(n->IsValid())
	if( n->children__!=MB_BADINDEX) {
		for(i=0; i<n->GetNumChildren(); i++) {
			sz += ComputeSize(n->children__ + i);
		}
	}
	return sz;
}

int Forest::ComputeSize(SVSetNodePtr node, int not_this_child)
{
	int sz=1;
	int i;
	SVSetNode Node;
	SVSetNode* n=MBSV.GetAt(node,&Node);
	if(n->IsValid())
	if( n->children__!=MB_BADINDEX) 
		for(i=0;i<n->GetNumChildren();i++) {
			if(i!=not_this_child) sz += ComputeSize(n->children__ + i);
	}
	return sz;
}


SVSetNodePtr Forest::Front()
{	return (SVSetNodePtr) bh_return_min(HeapP);
}

void	Forest::RemoveTop()
{ // use heap
	SVSetNodePtr head=(SVSetNodePtr)bh_delete_min(HeapP);
	SVSetNode Node;
	SVSetNode* Head=MBSV.GetAt(head, &Node);
	Head->SetHeapNode(1); //mynode=BADINDEX;
	MBSV.FreeBlock(head);
}

void	Forest::RemoveNodeFromHeap(SVSetNodePtr node)
{	// decrease key
	SVSetNode Node;
	SVSetNode* n=MBSV.GetAt(node,&Node);

	if(n->GetHeapNode()==1) return; //mynode==BADINDEX not in the heap

	n->SetHeapNode(1);//mynode=BADINDEX;
	MBSV.SetAt(node,n);
}

void	Forest::RemoveNodeFromHeap(SVSetNodePtr node, SVSetNode* Node)
{	Node->SetHeapNode(1);}


void	Forest::AddRootNode( SVSetNodePtr node)
{
// this method called only once in the initpopulate
	SVSetNode Node;
	SVSetNode* n=MBSV.GetAt(node,&Node);

	n->parent__=MB_BADINDEX;
	SET_ROOT(n->vecnumber);
	n->SetNumChildren(-1);

//	n->mynode=
	InsertNodeIntoHeap(node,n)
	n->SetHeapNode(0);
	MBSV.SetAt(node,n);

	HeadStruc Head;
	Head.Head=node;
	Head.p_index = new(siVec);
	*(Head.p_index) = temp_index;

	Heads.push_back(Head);
	size++;
}

void	Forest::AddLeaf(SVSetNodePtr node)
{
	// inserts the leaf into the heap, otherwise calls recursively
	int i;
	SVSetNode Node;
	SVSetNode* n= MBSV.GetAt(node, &Node);

	if(n->GetNumChildren()<=0 && n->IsValid()  /*==MB_BADINDEX*/) // means this is a leaf
	{
		//n->mynode = 
		Merge(HeapP,node);
		n->SetHeapNode(0);
		size++; // size of tree
		MBSV.SetAt(node,n);

		return;
	}
	else { // process children if any
		for(i=0;i<n->GetNumChildren();i++) {
			// the child can be empty
			 AddLeaf(n->children__ +i);   // where is memory allocation?
		}
		size++; 
		MBSV.SetAt(node,n);
	}
}

void	Forest::AddTree(SVSetNodePtr root)
{
	if(root==MB_BADINDEX) return;

	HeadStruc Head;
	Head.Head=root;
	Head.p_index = new(siVec);
	*(Head.p_index) = temp_index;
	Heads.push_back(Head);

	SVSetNode Node;
	SVSetNode* R= MBSV.GetAt(root, &Node);
	R->parent__=MB_BADINDEX;
	AddLeaf(root);
	MBSV.SetAt(root,R);
}



void	Forest::EraseBranch(SVSetNodePtr branch, int processparent)
{
	if(!IS_VALID_PTR((EUINT)branch)) return;

	SVSetNode Node;
	SVSetNode* Branch= MBSV.GetAt(branch, &Node);

	if(!Branch->IsValid()) return;

	int i;
	if(Branch->GetNumChildren()>0/*!=MB_BADINDEX*/) {
		for(i=0;i<Branch->GetNumChildren();i++)
			EraseBranch(Branch->children__ +i);
	}


	if(Branch->killers !=NULL  ) if(Branch->Dval>0 /*&& Branch->Dval!=Special1*/)	{delete Branch->killers; sizevirtual--;}
//	else if(IS_VALID_PTR((EUINT)(Branch->killers)))
//		 MBCM.FreeBlockC((EUINT)(Branch->killers));
	
	Branch->killers =NULL;
	size--;
	if(Branch->GetNumChildren()<0)
		RemoveNodeFromHeap(branch);
	if(processparent > -1) // if I want to remove reference from parent
		if(Branch->parent__!=MB_BADINDEX) {// if not root
			;//MBSV.FreeBlock(branch);
		}
		else // means root
			EraseRootEntry(branch);
	
	MBSV.FreeBlock(branch);	

	Branch->vecnumber=0xFFFFFFFF; // not root: indicates currupted
	Branch->parent__=0;
//	delete branch;
}

void	Forest::EraseRootEntry(SVSetNodePtr branch)
{
	// find the root in the list of roots and erase it.
	deque<HeadStruc>::iterator iter;

	for(iter=Heads.begin(); iter!=Heads.end(); iter++)
		if((*iter).Head == branch) {
			delete (*iter).p_index;
			iter=Heads.erase(iter); 
			return;
		}
}

void	Forest::ClearBranch(SVSetNodePtr branch)
{
// this method is called from destructor. It differes from EraseBranch in that
// the nodes are not deleted from the heap (to save time)
	if(!IS_VALID_PTR((EUINT)branch)) return;

	SVSetNode Node;
	SVSetNode* Branch= MBSV.GetAt(branch, &Node);
	if(!Branch->IsValid()) return;

	int i;
	if(Branch->GetNumChildren()>0/*!=MB_BADINDEX*/) {
		for(i=0;i<Branch->GetNumChildren();i++)
			ClearBranch(Branch->children__ +i);
	}
	if(Branch->killers !=NULL  ) if(Branch->Dval>0)	delete Branch->killers; 
//	else
//		if((EUINT)(Branch->killers) !=MB_BADINDEX)
//		 MBCM.FreeBlockC((EUINT)(Branch->killers));
	Branch->killers =NULL;

// destructor, keeps in the heap invalid reference
	size--;
}




void Forest::EraseAll()
{
	deque<HeadStruc>::iterator iter;

	if(MBSV.IsValid() )
		for(iter=Heads.begin(); iter!=Heads.end(); iter++) { 
			delete (*iter).p_index; ClearBranch((*iter).Head);  } //???????

//	we have to clear virtual nodes as well!!
//	MBCL.ClearAll();
	MBSV.ClearAll();
/**/
	Heads.clear();
	size=0;
	bh_free(HeapP); // release memory
}


// here interesting stuff--------------------------------
void	Forest::ProcessNodeConstraintsSimple(SVSetNodePtr node, SVSetNode* Node, siVec& index, dVec& soln)
{
	int i,j;
//	support_vector* SVPtr;
	support_vector *svt;
	ConstrMin* cm;
	ConstMinPtr cmP;
	RemoveNodeFromHeap(node, Node);

	double slack, f1,f2,f3,ft;
	slack=0;
	for(i=0;i<Parent->Dim-1;i++) // project to positive quadrant
	{
		soln[i]= max_(0,soln[i]);
		slack += soln[i];
	}
	soln[Parent->Dim-1]=1.0-slack;

	siVec* idxt =&(Parent->m_TempSiVec);

	(*idxt)=0;
	int last=0;

	f2=0;
	for(i=0;i<Parent->Dim-1;i++) // project to positive quadrant
	{
		if(soln[i] > Parent->m_bounds[i]) { 
			(*idxt)[i]=1;
			last++;
			f1=soln[i]-Parent->m_bounds[i];
			f2+=f1;
			soln[i]=Parent->m_bounds[i];
			for(j=0;j<Parent->Dim;j++)
				//if(j != i) {
				if((*idxt)[j]==0) {
					soln[j] += f2/(Parent->Dim - last);
			}
		}
	}


	// evaluate function value
	f1=-Infinity; f3=Infinity;
	for(i=0;i<Parent->Dim;i++) {
		svt = &(Parent->SVectors[index[i]]);
		f2=(svt->vec[i] + soln[i]) * Parent->m_Constants[i];
		if(f1<f2) f1=f2;
		if(f3>svt->funvalue) f3= svt->funvalue;
	}
		
	if(f1<f3 ) // OK, feasible  ??? was f2
	{
		ft=Node->Dval; // save it as well
		if(f1 < -GetConst()) f1= -GetConst()+1;
		Node->Dval= - (f1 + GetConst()); // negative means special case 

		cmP=MBCM.GetNextFree();
		cm=MBCM.GetAt(cmP);
		cm->SetVal(ft);
		cm->SetVec(soln);

		Node->killers = (ChildrenList*) cmP;

			 // insert into the heap only
		InsertNodeIntoHeap(node,Node);
		MBSV.SetAt(node,Node);
//		Node->mynode=Merge(HeapP,node);
		size_constrained++;
		return;
	}
	else // empty intersection
	{
		Node->Dval = Special1;
		Node->killers =0;
//		InsertNodeIntoHeap(node,Node);
		
		MBSV.SetAt(node,Node);
//		Node->mynode=Merge(HeapP,node);
		size_infeasible++;

		return;
	}
}




// this one ensures that the front does not have virtual children
int Forest::ProcessFront()
{
//	SVSetNodePtr head;
	SVSetNode Node;
//	SVSetNode* Head;
#ifdef VIRTUAL_NODES

	while(1) {
		head=(SVSetNodePtr) bh_return_min(HeapP);
		if (head==MB_BADINDEX) return 0;

		Head=MBSV.GetAt(head, &Node);
		if(Head->GetHeapNode() != 0) bh_delete_min(HeapP);
		else break;
	}
	if (Head->killers==NULL || Head->Dval<0 || Head->Dval==Special1) 
	{ 		return 0; 	}

		switch(ProcessTopNode(head,Head,&m_initvec)) {
			case 0:		return 0;
			case 1:
			case 2:
			default:	return 1; // try again
	}


#else  // no virtual nodes

	return 0;
#endif
}


int Forest::ProcessFront(SVSetNodePtr &head,	SVSetNode** Head)
{
#ifdef VIRTUAL_NODES

	head=(SVSetNodePtr) bh_return_min(HeapP);
	if (head==MB_BADINDEX) return 0;

	*Head=MBSV.GetAt(head, *Head);
	if((*Head)->GetHeapNode() != 0) { 
//		printf("delete heap %x %f\n",head, (*Head)->Dval);
		bh_delete_min(HeapP); return 1;
	}

	if ((*Head)->killers==NULL || (*Head)->Dval<0 || (*Head)->Dval==Special1) 
	{  		return 0; 	}

		switch(ProcessTopNode(head,*Head,&m_initvec)) {
			case 0:		return 0;
			case 1:
			case 2:
			default:	
				bh_delete_min(HeapP);  // OK, remove it from heap
				RemoveNodeFromHeap(head,*Head);
				MBSV.SetAt(head,*Head);
				return 1; // try again
	}


#else  // no virtual nodes

	head=(SVSetNodePtr) bh_return_min(HeapP);
	if (head==MB_BADINDEX) return 0;

	*Head=MBSV.GetAt(head, *Head);
	if((*Head)->GetHeapNode() != 0) { 
		bh_delete_min(HeapP); return 1;
	}

	return 0;
#endif
}

void CopyNode(SVSetNode* source , SVSetNode* dest)
{	memcpy(dest, source, sizeof(SVSetNode)); }

// splits the node
int Forest::ProcessTopNode(SVSetNodePtr node, SVSetNode* Node, siVec* initvec)
{
	// return: 0 means ready to eval f in it, 1 means it has been split, 2 means ==, special split 3 means it has been cut

	SVSetNode NodeObj;
	SVSetNode NodeObjTemp, NodeObjTemp1;
//	SVSetNode* Node=MBSV.GetAt(node,&NodeObj);

	int P,i,root, numchld;
	support_vector SV;
	real olddiag;
	SVSetNode* Newnode, *tempnode;
	SVSetNodePtr newnode;

	ChildrenList* inheritedkillers;

	deque <PossibleKillerStruc>::iterator iter;

	if(Node->killers->empty()) {
		delete Node->killers; 
		Node->killers=NULL;
		return 0;
	}

	// root can never be in the heap, unless it's ROOT. In this case killers is empty.
// this is questionable
	Node->ComputeVectors(m_index, initvec, node);

	int it=0;
// now special case: root, could have virtual branches: no worries, skip first N
	root=IS_ROOT(Node->vecnumber);
	
// now process possible killers
	while(it < Node->killers->size() )
	{
		SV=Parent->SVectors[   (*(Node->killers))[it]  ];

		P=Node->TestVectorIndex(SV.vec, &m_index); // explicitly test, do not remove!! list may be inherited
		switch(P) {
		case 3: // erase this branch
			// this case should never happen (at least for roots). Ensure this!
			sizevirtual -= 1;//Node->killers->size();

			delete Node->killers; Node->killers=NULL;
			EraseBranch(node,1); // this won't work for root  //?????
			MBSV.SetAt(node,Node);
			return 2;
			break;

		case 2:
			SV.Increment();  // very questionable
			// special split
			//break;

		case 1: // split into children
			it  = Node->killers->erase(it); 
			sizevirtual--;
			// start splitting

			numchld=0;

			for(i=0;i<Parent->Dim;i++)
				if(Node->TryNewVectorIndex(&SV, i, &m_index, olddiag)) { // cond (1)
					// add new node
					size++;
					newnode = m_TempChildren + numchld;
					numchld++;

					Newnode = MBSV.GetAt(newnode,  &NodeObjTemp);
					Newnode->Init();
					Newnode->vecnumber=POSVEC_VEC(i, SV.label);
					Newnode->SetNumChildren(-1);
					Newnode->parent__=node;

			Newnode->CalculateDval(node, Node, &SV, olddiag);

					if( ! Node->killers->empty() ) { // inhetit possible killers
						inheritedkillers = new ChildrenList();
						// special case: root ????? Supposedly root cannot have more than 1 killer ?
						// copy
						inheritedkillers->copy(Node->killers,it);
						Newnode->killers=inheritedkillers;
						sizevirtual++;
					} else Newnode->killers=NULL;

					MBSV.SetAt(newnode,Newnode);

				} 


				// no need to remember constr. minimum

			/*	if(Node->Dval<0 && Node->killers!=NULL) {
					MBCM.FreeBlockC((EUINT)(Node->killers));
					Node->killers=0;
				}*/

				Node->SetNumChildren(numchld);
				Node->children__= MBSV.GetNextFree(numchld);

				sizemem+=numchld;
				for(i=0;i<numchld;i++) { // copy to the actual nodes
					newnode=Node->children__ + i;
					tempnode= MBSV.GetAt(m_TempChildren + i, &NodeObjTemp1);
					Newnode = MBSV.GetAddress(newnode, &NodeObjTemp);
					CopyNode(tempnode, Newnode);

					InsertNodeIntoHeap(newnode,Newnode);

					MBSV.SetAt(newnode,Newnode);
				}


			if(!root) // not root
			{ 
				delete (Node->killers);  Node->killers=NULL;	
//				MBSV.SetAt(node,Node); // will be done on exit, with erasing from the heap
			}
			return 1;
			break;

		default:
			it=Node->killers->erase(it); // try next killer, but discard this one
			//sizevirtual--;
			break;
		}
	}


	//if none succeeded, clear the list and split this node
	if(!root) {
		//sizevirtual -= node->killers->size();
		// decreased in default!
		delete Node->killers;
		Node->killers=NULL;
		sizevirtual--;
	}

	MBSV.SetAt(node,Node);
	return 0; // means we process this leaf
}


// global
	SVSetNodePtr newnode;
	real olddiag;
	SVSetNode *Newnode,*tempnode;
	EUINT Killer;
	int AtLeastOneFound;
// global (not in recursion)
int Globvecnumber_;

int Forest::ProcessNode(SVSetNodePtr node, support_vector* SV, siVec* initvec)
{
	// recursive calls
	SVSetNode* Node;
	int P,i,numchld;

	int NC;
	SVSetNodePtr CHLD;

	SVSetNode NodeObj;
	SVSetNode NodeObjTemp, NodeObjTemp1;

	if(node!=MB_BADINDEX)	Node=MBSV.GetAt(node,&NodeObj);
	else return 0;
	if(!Node->IsValid()) return 0;

	// here I can keep indexvector (starting from the top), so no computevectors is necessary.
	P=Node->TestVector(SV->vec, initvec);

	switch(P) {
	case 2: // special split, very questionable
			// remember position to restore after return
//			Pos=GlobPos;
		
	case 1: // dominance, split this node
		if(Node->GetNumChildren()<=0 /*Node->children__ == MB_BADINDEX*/) { // means leaf
			AtLeastOneFound++;
#ifdef VIRTUAL_NODES
			if(Node->Dval>0 && Node->Dval!=Special1 && Parent->AllowVirtual() ) {
				// delay splitting, even if this is root
				if(Node->killers==NULL) {Node->killers=new ChildrenList;  sizevirtual++;}
		
				Killer = SV->label;
				Node->killers->push_back(Killer);
				
				goto L1;
			}
#endif


#ifdef SHORT_HEAP
			if(fabs(Node->Dval) >= m_heapworst) // do nothing
				break;
#endif


			if(P==2 /*|| ISIDXSET()*/) 
				SV->Increment(); // questionable ??????????

			// must be no children__ at this stage
			// create children__ if any


			Node->ComputeVectors(m_index,initvec,node);
			numchld=0;

			for(i=0;i<Parent->Dim;i++)
				if(Node->TryNewVectorIndex(SV, i, &m_index, olddiag)) { // cond (1)

					// add new node
					size++;
					newnode = m_TempChildren + numchld;
					numchld++;

					Newnode = MBSV.GetAt(newnode,&NodeObjTemp);
					Newnode->Init();
					Newnode->vecnumber=POSVEC_VEC(i, SV->label);
					Newnode->SetNumChildren(-1);
					Newnode->parent__=node;

					if(Node->Dval!=Special1)  Newnode->CalculateDval(node, Node, SV, olddiag);
					else Newnode->CalculateDvalComplete(node, Node, SV, olddiag, m_index);

					MBSV.SetAt(newnode,Newnode);
				} 

				RemoveNodeFromHeap(node,Node); // node became branch, remove from heap.

				// free constrained minimum of any (should be no killers by now)
//				if(Node->Dval<0 && IS_VALID_PTR((EUINT)(Node->killers)) && numchld>0) {
//					MBCM.FreeBlockC((EUINT)(Node->killers));
//					Node->killers=0;
//				}


				Node->SetNumChildren(numchld);
				Node->children__= MBSV.GetNextFree(numchld);
				sizemem+=numchld;
				for(i=0;i<numchld;i++) { // copy to the actual nodes
					newnode=Node->children__ + i;
					tempnode= MBSV.GetAt(m_TempChildren + i,&NodeObjTemp);
					Newnode = MBSV.GetAddress(newnode, &NodeObjTemp1);
					CopyNode(tempnode, Newnode);

					InsertNodeIntoHeap(newnode,Newnode);

					MBSV.SetAt(newnode,Newnode);
				}
			
L1:
			MBSV.SetAt(node,Node);

		} else { // this was a branch, recursively process children
			NC=Node->GetNumChildren();
			CHLD=Node->children__;

			for(i=0;i<NC;i++) 
			{
				if( ProcessNode(CHLD + i, SV, initvec)==3) // everything is done here
				// need to erase this branch, Dval too small
				{				}
			}
		}
		// on exit undo IDX
		break;
	case 3: // cut the branch, the min is >= frecord
		return 3;
		break;
	case 0: //test (2) passed
		// do nothing
		// this is critical piece. If I proceed as in case 1, no branches are cut and
		// the algorithm becomes exponential
		break;
	}

/*------------------- for indices --*/

	return 0;
}

void Forest::ProcessAll(support_vector* SV)
{
	deque<HeadStruc>::iterator iter;
//	SVSetNode *node;

	SVSetNode NodeObj;

	AtLeastOneFound=0;

	for(iter=Heads.begin(); iter!=Heads.end(); iter++) {

// the root will code its initial vector in the possible killers list
// in predefined order. One exception is ROOT at iteration 1. This case is handled in
// GenerateInitVector. for all other roots ensure this is the case (cannot be virtual node)
//		node=MBSV.GetAt( (*iter).Head, &NodeObj );
		m_initvec = *((*iter).p_index);
		if( ProcessNode( (*iter).Head, SV, &m_initvec) == 3 )
		{ // permanently remove, above best known value

			iter=Heads.erase(iter);
			if(Heads.empty()) break;
		}
	}
}


int Forest::ProcessNodeQ(SVSetNodePtr node, support_vector* SV, siVec* initvec)
{
	// recursive calls
	SVSetNode NodeObj;
	SVSetNode* Node=MBSV.GetAt(node, &NodeObj);

	if(!Node->IsValid()) return 0;

	int P,i;
	// here I can keep indexvector (starting from the top), so no computevectors is necessary.	
	P=Node->TestVectorQ(SV->vec, initvec);

	switch(P) {
	case 1: // dominance, split this node
		if(Node->GetNumChildren()<0/*Node->children__ == MB_BADINDEX*/) { // means leaf
			m_ListAffectedNodes.push_back(node);
		} else { // this was a branch, recursively process children
			for(i=0;i<Node->GetNumChildren();i++) 
			{
				 ProcessNodeQ(Node->children__ +i, SV, initvec);
			}
		}
		break;
	case 2: // special split
		break;
	case 3: // cut the branch, the min is >= frecord
		return 3;
		break;
	case 0: //test (2) passed
		// do nothing
		// this is critical piece. If I proceed as in case 1, no branches are cut and
		// the algorithm becomes exponential
		break;
	}
	return 0;
}
	
void Forest::ProcessAllQ(support_vector* SV)
{
	m_ListAffectedNodes.clear();
	deque<HeadStruc>::iterator iter;
	SVSetNode NodeObj;
	SVSetNode *node;

	for(iter=Heads.begin(); iter!=Heads.end(); iter++) {
		node=MBSV.GetAt( (*iter).Head, &NodeObj );
		m_initvec = *((*iter).p_index);
		ProcessNodeQ( (*iter).Head, SV, &m_initvec);
	}
}




double	Forest::GetConst()
{ 	return Parent->GetConst(); }


void	Forest::ComputeHeapWorst(int n, int extra)
{
	// will resort the heap and exclude all bu extra first elements
	int heapsize=HeapP->nodeCount;
	DATA_TYPE el;
	DATA_TYPE* array=(DATA_TYPE*)malloc(sizeof(DATA_TYPE)*extra);
	int total=MMmin(n,heapsize);

	SVSetNodePtr node;
	SVSetNode* Node;
	SVSetNode  NodeSpace; 

	int i=0;
	while(i<total && HeapP->nodeCount>0) {
		el=bh_delete_min(HeapP);
		node=el & 0xFFFFFFFF;
		if(node !=BADINDEX) { 
			Node=MBSV.GetAt(node, &NodeSpace);
			if(Node->GetHeapNode() == 0) { //only if it was not deleted earlier
				array[i]=el;
				i++;
			}
		}
	}
	m_heapworst=- _getKey(array[i-1]);

	// continue with extra
	while(i<extra && el !=BADINDEX) {
		el=bh_delete_min(HeapP);
		node=el & 0xFFFFFFFF;
		if(node !=BADINDEX) { 
			Node=MBSV.GetAt(node, &NodeSpace);
			if(Node->GetHeapNode() == 0) { //only if it was not deleted earlier
				array[i]=el;
				i++;
			}
		}
	}
// remove the rest
// here we can probably just destroy the heap, no need to re-sorting
	bh_free(HeapP);
	HeapP=bh_alloc();
//	while(HeapP->nodeCount>0 && el !=BADINDEX) {
//		el=bh_delete_min(HeapP);
//	}

// repopulate the heap
	total=i;
	for(i=0;i<total;i++)
		 bh_insert(HeapP, array[i]); 

	delete(array);
}




ChildrenList::ChildrenList() 
{
	FirstChild=MBCL.GetNextFree();

	First=Last=0;
	LastChild=FirstChild;
	Children* block=MBCL.GetAt(FirstChild);
	block->NextBlock=MB_BADINDEX;
}
ChildrenList::~ChildrenList() 
{
	if(FirstChild==MB_BADINDEX) return;
	Children* block=MBCL.GetAt(FirstChild);
	EUINT Next=block->NextBlock;
	MBCL.FreeBlock(FirstChild);
	FirstChild=Next;

	while(FirstChild!=MB_BADINDEX) {
		block=MBCL.GetAt(FirstChild);
		Next=block->NextBlock;
		MBCL.FreeBlock(FirstChild);
		FirstChild=Next;	
	}
}
EUINT ChildrenList::AddBlock()
{
	return MBCL.GetNextFree();
}

void ChildrenList::push_back(EUINT v)
{
	Children* block;
	if(Last==0) {
		block=MBCL.GetAt(FirstChild);
		block->m_children[0]=v; Last++; return;}

	div_t b=div(Last, ChildrenBlockSize);
	b.rem;
	block=MBCL.GetAt(LastChild);
	if(b.rem!=0) {
		block->m_children[b.rem]=v; Last++;
	}
	else {
		block->NextBlock=AddBlock();
		LastChild=block->NextBlock;
		block=MBCL.GetAt(LastChild);
		block->m_children[0]=v;
		block->NextBlock=MB_BADINDEX;
		Last++;
	}	
}

EUINT ChildrenList::pop_back()
{
	Children* block, *block1;
	div_t b=div(Last-1, ChildrenBlockSize);
	block=MBCL.GetAt(LastChild);
	EUINT r=block->m_children[b.rem]; // to return
	// adjust size
	Last--;

	b=div(Last, ChildrenBlockSize);
	if(Last>0 && b.rem == 0) //need to remove last block
	{
		block1=block=MBCL.GetAt(FirstChild);
		if(block->NextBlock==LastChild) {
			LastChild=FirstChild;
			MBCL.FreeBlock(block->NextBlock);
		} else {
			block=GetMemBlockAt(block->NextBlock);
			while(block->NextBlock!=LastChild) {
				block1=block;
				block=GetMemBlockAt(block->NextBlock);
			}
			LastChild=block1->NextBlock;
			MBCL.FreeBlock(block->NextBlock);
		}
	}
	return r;
}

EUINT ChildrenList::pop_front()
{
	Children* block=MBCL.GetAt(FirstChild);
	EUINT r=block->m_children[First];
	First++;
	if(First==ChildrenBlockSize) {
		// adjust size: remove first block
		EUINT B=FirstChild;
		FirstChild=block->NextBlock;
		MBCL.FreeBlock(B);
		First=0; Last -= ChildrenBlockSize;
	}
	return r;
}

EUINT ChildrenList::erase(int i)
{
	// swap with the last one. does not preserve order. do not use for the first Dim entries
	// returns next i (that is, this same i)
//	EUINT r=pop_back();
//	SetAt(i,r);
	pop_front();
	return i;
}
inline EUINT& ChildrenList::operator[](int i) { Children* block=MBCL.GetAt(FirstChild); return (*block)[i+First]; }
void ChildrenList::SetAt(int i, EUINT val) {Children* block=MBCL.GetAt(FirstChild); (*block)[i+First]=val; } 

void ChildrenList::copy(ChildrenList* source, int start)
{
	for(int i=start;i<source->size();i++) push_back( (*source)[i] );
}
