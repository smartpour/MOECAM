
#include "ecam.h"
#include "../ganso.h"



 double Meps=1.0e-5;

 // just to time the algorithm
//time_t	start, stop;

SawToothCover*	Parent;

#ifdef MPI
extern MemoryBlockMPI<SVSetNode>	MBSV;
#else
extern MemoryBlock<SVSetNode>	MBSV;
#endif

extern MemoryBlockE<ConstrMin>   MBCM; // blocks for ConstMin array

double TVal=0;

/* implementation --------------------------------------- */

SawToothCover::SawToothCover()
{	//A.newsize(0); B.newsize(0); 
	Iters=MaxIter;
	Parent=this;
	frecord=Infinity;
	GVA=-Infinity;
	C_minus=1;
	m_Constants.newsize(1);
	m_Constants=1;
	m_bounds.newsize(1);
	m_SimplexBound=1; // standard simplex by default
	logfile=NULL;
	FinalList=NULL;
	ListSize=0;
	m_myobj = 0;
	Cut = _Cut;
	Const = _Const;

}; // constructor

SawToothCover::~SawToothCover()
{
	if(FinalList!=NULL) free(FinalList);
}


double SawToothCover::RequestGap()
{
	return (frecord - GVA)/fabs(GVA);
}
// Class Support Vector -------------------------------------------------------

support_vector* support_vector::This()   {	return this;}

void support_vector::PrintMe(FILE *f)
{
	if(f==NULL ) return;
	int i;
	for(i=0;i<vec.size();i++)
		fprintf(f," %.14le",funvalue/Parent->m_Constants[i] - vec[i]);
//	fprintf(f," %.14le %u\n",funvalue,label);
	fprintf(f,"\n");
	fflush(f);
}


void support_vector::SVForm(dVec& x, real val)
// forms the support vector from x and value
{
	vec.newsize(x.size());
	funvalue=val;
	for(int i=0;i<vec.size();i++) 
		vec[i]=funvalue/Parent->m_Constants[i] - x[i];
}

void support_vector::Increment()
{
	for(int i=0;i<vec.size();i++) 
		vec[i]+= (max_(funvalue,1)*Meps)/Parent->m_Constants[i];
	funvalue += max_(funvalue,1)*Meps;
}

short int	support_vector::ChangeF(real newval)
{
	short int i;
	for( i=0;i<vec.size();i++) 
		vec[i]+= (-funvalue + newval)/Parent->m_Constants[i];
	i= (newval>=funvalue);
	funvalue=newval;
	return i;
}

void support_vector::ReturnX(dVec& x)
{
	for(int i=0;i<vec.size();i++)
		x[i]=funvalue/Parent->m_Constants[i] - vec[i];
}

void SawToothCover::InitialSteps(int d, int iters)
{
	Parent = this;
	Dim = d + 1; // we are on simplex
	Iters = iters;
	radius = 0;

	if (m_bounds.size() < Dim - 1)
	{
		m_bounds.newsize(Dim - 1); m_bounds = 1;
	}

	//	div_t dt;

	// assume constants are given
	ComputeCminus();

	HeapPossibleMin.Init();
	HeapPossibleMin.SetHeapLimit(Iters + 100);

	m_lastindex.newsize(Dim);
	m_TempSiVec.newsize(Dim);

	// this is the main loop
	GVA = -Infinity;

	//	ResetTime();
	iteration = 0;

	if (logfile != NULL) {
		fprintf(logfile, "Format: #vectors, #auxMin, #nodes #virtual #total  duration, frecord, GVA\n");
		fprintf(logfile, "Const: %f\n", Const);
	}

	SVectors.clear(); InitPopulate();

	// if need to shift the Aux function up (only for the heap!!)
	ComputeConstant();

	AddVertices();

}

int SawToothCover::MakeStepbefore()
{
	SVSetNode NodeObj;
	SVSetNode *Node;
	SVSetNodePtr node;
	int res;


L3:		Node = &NodeObj;
	// wait till we eliminate virtual branches and have acceptable node
	while (HeapPossibleMin.ProcessFront(node, &Node)) //WILL EXIT itself
	{
		Node = &NodeObj;// reset the address	
	}

	if (node == MB_BADINDEX) return 1;

	res = AssessCandidate(Node, node);
#ifdef PERIODIC
	res = 1;
#endif // periodic
#ifdef ANTIPERIODIC
	res = 1;
#endif // antiperiodic

	if (!res) { // if we need to modify this node as not satisfying constraints
		bh_delete_min(HeapPossibleMin.HeapP);  // OK, remove it from heap and reinsert

		HeapPossibleMin.ProcessNodeConstraintsSimple(node, Node, m_lastindex, temp_x);
		//				HeapPossibleMin.RemoveNodeFromHeap(node);
		MBSV.SetAt(node, Node);
		goto L3; // start all over
	}
LC:

	if (CheckExistingSV(temp_x)) // perturb ?
	{
		HeapPossibleMin.RemoveNodeFromHeap(node);	goto L3;
	}


	return 0;
}

int SawToothCover::MakeStep(dVec&  X) 
{
	temp_x = X;
	SVSetNodePtr node;


	// OK, evaluate f() and form the SV
	FormNewSV();


	HeapPossibleMin.ProcessAll(&(SVectors.back()));

	if (HeapPossibleMin.HeapP->nodeCount > Iters + 100 && HeapPossibleMin.m_heapworst == limitworst) {

		HeapPossibleMin.m_heapworst = HeapPossibleMin.GetGlobalWorstHeap(); // start screening nodes

	}
	else if (HeapPossibleMin.HeapP->nodeCount > HeapPossibleMin.m_heapoverlimit) {
		HeapPossibleMin.ComputeHeapWorst(HeapPossibleMin.m_heaplimit - LastLabel, HeapPossibleMin.m_heaplimit - LastLabel);
		//HeapPossibleMin.m_heaplimit*1.1);
	}

	return 0;
}

void SawToothCover::Minimise(int d, int iters)
{
	m_lasterror=0;
#ifdef MPI
	int myid;
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	if(myid!=0) {
		MBSV.SlaveMain();
		return;
	} else
	MBSV.Init();
#endif

//		scanf("%d",&myid);


	int res;
	int sz1,sz2,sz3,sz4;
	SVSetNode NodeObj;
	SVSetNode *Node;
	SVSetNodePtr node;

	int objchosen;

	while(StopCriteria()==0) {

		MyParent->MakeStepbefore();

		objchosen = MyParent->DecideObjectives();

	//	printf("%d \n", objchosen);

		MyParent->MakeStep(objchosen);
		 


		/*
L3:		Node=&NodeObj;
// wait till we eliminate virtual branches and have acceptable node
		while(HeapPossibleMin.ProcessFront(node,&Node)) //WILL EXIT itself
		{ Node=&NodeObj;// reset the address	
		}

		if(node == MB_BADINDEX) break;
		res=AssessCandidate( Node,node ); 
#ifdef PERIODIC
		res=1;
#endif // periodic
#ifdef ANTIPERIODIC
		res=1;
#endif // antiperiodic

		if( !res ) { // if we need to modify this node as not satisfying constraints
			bh_delete_min(HeapPossibleMin.HeapP);  // OK, remove it from heap and reinsert

			HeapPossibleMin.ProcessNodeConstraintsSimple(node,Node,m_lastindex,temp_x); 
//				HeapPossibleMin.RemoveNodeFromHeap(node);
			MBSV.SetAt(node, Node);
			goto L3; // start all over
		}
	LC:             

		

		if(CheckExistingSV(temp_x)) // perturb ?
		{ 	HeapPossibleMin.RemoveNodeFromHeap(node);	goto L3;	}

// OK, evaluate f() and form the SV
		FormNewSV();

		*/

// we have the index, so just process the top node: will be faster
//		HeapPossibleMin.ProcessTopNodeWithSV(node,&m_lastindex,&(SVectors.back()));


/*------These are all aux reporting functions  ---------*/
		sz1=HeapPossibleMin.Size() & 0x00FFFFFF;
		sz2=HeapPossibleMin.GetSize() ;
		sz4=HeapPossibleMin.GetSizeMem() ;

/*
// now process the rest
		HeapPossibleMin.ProcessAll(&(SVectors.back()));

		if(HeapPossibleMin.HeapP->nodeCount > Iters+100 && HeapPossibleMin.m_heapworst==limitworst) {

			HeapPossibleMin.m_heapworst=HeapPossibleMin.GetGlobalWorstHeap(); // start screening nodes

		} else if(HeapPossibleMin.HeapP->nodeCount > HeapPossibleMin.m_heapoverlimit) {
			HeapPossibleMin.ComputeHeapWorst(HeapPossibleMin.m_heaplimit-LastLabel, HeapPossibleMin.m_heaplimit-LastLabel);
				//HeapPossibleMin.m_heaplimit*1.1);
		}

		*/
		sz3=HeapPossibleMin.sizevirtual;

		if(logfile!=NULL)
		{fprintf(logfile," %d %d %d %d %d | %f %f %f %f\n",SVectors.size(),sz1,sz2,sz3,sz4,
				duration,frecord,GVA,flast); //#endif
			fflush(logfile);
		}


		iteration++;

/*------These are all aux reporting functions  ---------*/


	}

#ifdef PERIODIC
			ExtendPeriodic(xglob_x);
#endif // periodic
#ifdef ANTIPERIODIC
			ExtendAntiPeriodic(xglob_x);
#endif // antiperiodic

			if(ListSize>0) {
				FinalList=(double*)malloc(sizeof(double)*(Dim-1)*ListSize);
				PrepareFinalList(FinalList,&ListSize);
			}


	HeapPossibleMin.EraseAll();
return;
}



void SawToothCover::ExtendPeriodic(dVec& x)
{
	int i;
	double slack=0;
	for(i=0;i<Dim-1;i++) {
		if(x[i]>m_bounds[i]) x[i] -= m_bounds[i];
		else if(x[i]<0) x[i]+= m_bounds[i];

		slack += x[i];
	}
	x[Dim-1]=1.0-slack;
}

void SawToothCover::ExtendAntiPeriodic(dVec& x)
{
	int i;
	double slack=0;
	for(i=0;i<Dim-1;i++) {
		if(x[i]>m_bounds[i]) x[i] = 2*m_bounds[i] - x[i];
		else if(x[i]<0) x[i] = -x[i];

		slack += x[i];
	}
	x[Dim-1]=1.0-slack;
}


int SawToothCover::PopNextCandidate(dVec &x)
{
	if(HeapPossibleMin.Empty()) return 0;

	SVSetNode NodeObj;
	SVSetNode *Node;
	SVSetNodePtr node;
	int res;
L3:		if(!HeapPossibleMin.Empty()) {
			Node=MBSV.GetAt(node=HeapPossibleMin.Front(),&NodeObj);
			res=AssessCandidate( Node, node) ; 
		}
		else return 0;

		if( !res ) { 
#ifdef BOXCONSTRAINTS
			HeapPossibleMin.ProcessNodeConstraintsSimple(node,Node,m_lastindex,temp_x); 
#endif
			MBSV.SetAt(node,Node);
			goto L3;
		}

		MBSV.Discard(node,Node);

	if(HeapPossibleMin.Empty()) return 0;
	HeapPossibleMin.RemoveTop(); // forget this one
	x=temp_x;
	return 1;
}



// extended version of the above
int SawToothCover::PopNextCandidate(dVec& x, double& d, dVec& RHS)
{
	if(HeapPossibleMin.Empty()) return 0;

	SVSetNode NodeObj;
	SVSetNode *Node;

	SVSetNodePtr node;
	int res;
L3:		if(!HeapPossibleMin.Empty()) {
			Node=MBSV.GetAt(node=HeapPossibleMin.Front(),&NodeObj);
			res=AssessCandidate( Node, node) ; 
		}
		else return 0;

		if( !res ) { 
//#ifdef BOXCONSTRAINTS
			HeapPossibleMin.ProcessNodeConstraintsSimple(node,Node,m_lastindex,temp_x); 
//#else
//			HeapPossibleMin.ProcessNodeConstraints(node,Node,m_lastindex); 
//#endif
			MBSV.SetAt(node,Node);

			goto L3;
		}
		MBSV.Discard(node,Node);

	if(HeapPossibleMin.Empty()) return 0;
	HeapPossibleMin.RemoveTop(); // forget this one
	x=temp_x;

	int i;
	for(i=0;i<Dim;i++) { 
		RHS[i]=  SVectors[m_lastindex[i]].vec[i] * m_Constants[i];  // diag * constants
	}

	d = fabs(Node->Dval) - GetConst();

	return 1;
}



int		SawToothCover::StopCriteria()
{
	if (LastLabel>Iters) return 1;
	if(GVA>frecord ) { // need to stop: the constant was too small
			m_lasterror=-2;
			return 2;
	}
	return 0;
}

int		SawToothCover::LastError() {return m_lasterror;}

void	SawToothCover::ComputeCminus()
{
	Parent=this; // just in case for the interpolant
	int i,j;
// just to ensure right consts
	if(m_Constants.size()!=Dim) {
		dVec temp=m_Constants;
		m_Constants.newsize(Dim);
		j=min_(temp.size(),Dim);
		for( i=0;i<j;i++)
		if(temp[i]>0)
			m_Constants[i]=temp[i];
		else m_Constants[i]=1;
		for(i=j;i<Dim;i++) m_Constants[i]=m_Constants[j-1];
	}

	C_minus=0;
	for( i=0;i<Dim;i++) C_minus += 1.0/m_Constants[i];
}

void	SawToothCover::ComputeConstant()
{
	if(GVA - 2*(fabs(GVA)) + GetConst()  > 0) return; // constant is satisfactory, all local minima >0

	double t = -GVA + 1; // the thing we add to all f
	Const += t; 
	// update the first n vectors
	int i;
	;
	SVSetNode NodeObj;
	SVSetNode *S;

	S=MBSV.GetAt(HeapPossibleMin.Front(),&NodeObj);


	S->Dval=0;
	for(i=0;i<Dim;i++)	S->Dval += SVectors[i].vec[i];
	S->Dval = (S->Dval+1)/C_minus;

	S->Dval += GetConst();

	MBSV.SetAt(HeapPossibleMin.Front(),S);
}

void SawToothCover::AddVertices()
{
// add all vertices of the hypercube

	TVal=100; // allow virtual
//	return;

	support_vector	sv;
	int j;
	double val;
	unsigned long int mask;


	unsigned long int limit=1<<(Dim-1);
	double slack=0;

	int fl=0;
	int sz1,sz2,sz3;

// temporary fix just one vertex
	slack=0;

		for(j=0;j<Dim-1;j++) {
			temp_x[j]=m_bounds[j];
			slack += temp_x[j];
		}
		temp_x[Dim-1]=1-slack;

		val=Value(temp_x);
		if(frecord>val) { frecord=val; xglob_x=temp_x;} // copy vector
		sv.SVForm(temp_x,val);
		sv.label=LastLabel;
		SVectors.push_back(sv);
		LastLabel++;
		HeapPossibleMin.ProcessAll(&(SVectors.back()));

		return;

	for(mask=1;mask < limit;mask++)
	{
		slack=0;
		fl=0;
		for(j=0;j<Dim-1;j++) {
			if((mask>>j) &0x1) {temp_x[j]=m_bounds[j];fl++;} else temp_x[j]=0;
			slack += temp_x[j];
		}
		temp_x[Dim-1]=1-slack;
		
		if(fl==1) goto Cont; // already this vertex there in the basis
		val=Value(temp_x);
		if(frecord>val) { frecord=val; xglob_x=temp_x;} // copy vector
		sv.SVForm(temp_x,val);
		sv.label=LastLabel;
		SVectors.push_back(sv);
		LastLabel++;
		HeapPossibleMin.ProcessAll(&(SVectors.back()));


		sz1=HeapPossibleMin.Size() & 0x00FFFFFF;
		sz2=HeapPossibleMin.GetSize() ;
		sz3=HeapPossibleMin.sizevirtual;
		if(logfile!=NULL)
		{fprintf(logfile," %d %d %d %d  %f\n",SVectors.size(),sz1,sz2,sz3,
				frecord); //#endif
			fflush(logfile);
		}
Cont:;
	}

		TVal=100; // allow virtual
}
int		SawToothCover::AllowVirtual() {return TVal>1; }



dVec TempVector;
double TempA;

double SawToothCover::Value(dVec& x)
{
	 double FunValue;

	 FunValue = radius;
#ifdef PERIODIC
			savedtemp_x=x;
			ExtendPeriodic(savedtemp_x);
			MyParent->CurrentFunction(Dim-1,&(savedtemp_x[0]),&FunValue);
				//fv(&(savedtemp_x[0]),FunValue);
			goto LC;
#endif // periodic
#ifdef ANTIPERIODIC
			savedtemp_x=x;
			ExtendAntiPeriodic(savedtemp_x);
			MyParent->CurrentFunction(Dim-1,&(savedtemp_x[0]),&FunValue, &m_myobj);
//			MyParent->fv(&(savedtemp_x[0]),FunValue);
			goto LC;
#endif // antiperiodic

			MyParent->CurrentFunction(Dim-1,&(savedtemp_x[0]),&FunValue, &m_myobj);
//	MyParent->fv(&(x[0]),FunValue);

LC:
	double	fboxbasic = MyParent->m_box1;
	double	fboxnonbasic =MyParent->m_box2;
	double  fbox =  fboxbasic+fboxnonbasic;

	return ValueNormalise( FunValue);
}


void	SawToothCover::ReturnSol(dVec &x)
{
#ifdef SIMPLEXCONSTRAINTS
	for(int i=0;i<Dim;i++)
		x[i]= (-xglob_x[i] + 1)*TempA;
#else
		x=xglob_x;
#endif

}

int SawToothCover::PopNextCandidateExternal(dVec &x)
{
	PopNextCandidate(x);
#ifdef SIMPLEXCONSTRAINTS
	// here we need to translate the vector, assuming we have a simplex Sum x_i <= A
	for(int i=0;i<Dim;i++)
		x[i]= (-x[i] + 1)*TempA;
#endif
	return 1;
}

void SawToothCover::PrepareFinalList(double* List, int* size)
{
	int m,i,j,k;
	k=1;
	m=0;
	SVSetNode NodeObj;
	SVSetNode *Node;
	SVSetNodePtr node;

	while(k) {
// get next minimum

		Node=&NodeObj;
// wait till we have acceptable node
		while(HeapPossibleMin.ProcessFront(node,&Node)) //WILL EXIT itself
			 Node=&NodeObj;// reset the address	
		
		if(node != MB_BADINDEX) 
			AssessCandidate( Node,node ); 

		if(m>=ListSize || node == MB_BADINDEX)
		{k=0; goto L1;}  // break the loop
			bh_delete_min(HeapPossibleMin.HeapP); 

#ifdef PERIODIC
			savedtemp_x=temp_x;
			ExtendPeriodic(savedtemp_x);
			for(j=0;j<Dim-1;j++) List[m*(Dim-1) + j] = savedtemp_x[j];
#endif // periodic
#ifdef ANTIPERIODIC
			savedtemp_x=temp_x;
			ExtendAntiPeriodic(savedtemp_x);
			for(j=0;j<Dim-1;j++) List[m*(Dim-1) + j] = savedtemp_x[j];
#endif // antiperiodic
		m++;
L1:;
	}
	*size=m;
}


double	SawToothCover::ValueNormalise(double val)
{
	if(val>Cut) val = Cut;
//	val += Const;
	return val;
}

//extern int LU_factor( dMat &A, iVec &indx);

void	SawToothCover::InitPopulate() // vertices of simplex
//?????????? other points?
{

#ifdef SIMPLEXCONSTRAINTS
	TempVector.newsize(Dim);
	TempA=1.0/(Dim-1.0);
	m_bounds=1; // we need this for unit simplex
#endif

	LastLabel=0;
	FunCallsWhenFound=Dim;
	temp_x.newsize(Dim);

	temp_x=0.0;
	int i,j;
	double val;
	support_vector	sv;
	frecord=Infinity;

// not vertices, prepare correct locations
	TNT::Array2D<real> Atemp(Dim-1,Dim-1);
	TNT::Array1D<real> RHStemp(Dim-1);
	iVec ipiv(Dim-1);

	Atemp=m_Constants[Dim-1];
	for(i=0;i<Dim-1;i++)
			Atemp[i][i]=Atemp[i][i] + m_Constants[i];

	JAMA::LU<real> MyLU(Atemp);

//	if(LU_factor(Atemp,ipiv)==1) exit(1); // error

	double rt=0;
	for(j=0;j<Dim-1;j++)
		rt += m_bounds[j];
	
	rt *= m_Constants[Dim-1];

	for(j=0;j<Dim-1;j++)
	{
		RHStemp=rt;
		RHStemp[j] += m_Constants[j]*m_bounds[j];

		RHStemp=MyLU.solve (RHStemp);
	//	LU_solve(Atemp,ipiv,RHStemp);

		val=0;	

RHStemp[j]=m_bounds[j]*(1.0+0.2/m_Constants[j]); //?????????this is a hack!!!!

		for(i=0;i<Dim-1;i++){
			temp_x[i]=RHStemp[i];
			val += RHStemp[i];
		}
		temp_x[Dim-1] = 1-val;

		xglob_x=temp_x;
		val=Value(temp_x);
		xglob_val=val;
		if(frecord>val) { frecord=val; xglob_x=temp_x;} // copy vector
		sv.SVForm(temp_x,val);
		sv.label=LastLabel;
		SVectors.push_back(sv);
		LastLabel++;
	}

	// last one
	temp_x=0;
	temp_x[Dim-1]=1;
	val=Value(temp_x);
	if(frecord>val) { frecord=val; xglob_x=temp_x;} // copy vector
	sv.SVForm(temp_x,val);
	sv.label=LastLabel;
	SVectors.push_back(sv);
	LastLabel++;


	for(i=0;i<Dim-1;i++)
	for(int j=0;j<Dim;j++)
	{
		if(i!=j)
		 SVectors[i].vec[j] =Infinity;
	}

#ifdef SIMPLEXCONSTRAINTS
	for(i=0;i<Dim-1;i++)
	SVectors[Dim-1].vec[i]=Infinity;
#endif
	// same as above

// initial step is different, but similar to
//	AddGoodCandidates();

	// now generate new SVset, with every support vector we just constructed
	SVSetNodePtr n=MBSV.GetNextFree();

	SVSetNode NodeObj;
	SVSetNode*	S=MBSV.GetAt(n,&NodeObj);

	S->Init();
		//new SVSetNode;
	S->Dval=0;
	for(i=0;i<Dim;i++)	S->Dval += SVectors[i].vec[i];
	S->Dval = (S->Dval+1)/C_minus;

	GVA=S->Dval;

	S->Dval += GetConst();

	MBSV.SetAt(n,S);

	for(i=0;i<Dim;i++)	HeapPossibleMin.temp_index[i]=i;
	HeapPossibleMin.AddRootNode(n);

	LastLabel=Dim;
}



int SawToothCover::AssessCandidate(SVSetNode* SH, SVSetNodePtr thisnode)
{

// Here the node is assessed. Firstly, the index is generated. Then it is checked
// for the domain (not too close to the boundary, perhaps constraints?
// then temp_val, temp_x and tempFval are assigned. They are used in FormSVector.

	int i,j;
	double d;

	int tr=0;
	ConstrMin* cm=NULL;

	//GB check this, why it is here?

// m_lastindex remains valid on exit. Use it!
	SH->ComputeVectors(m_lastindex, HeapPossibleMin.GetVecAddress(), thisnode);

// Here process special case of constrained minumum
	if(SH->Dval == Special1) {
		// we can't be here! report error
		//cout<<"error"<<endl;
		// Yes we can, if we have a short list of minima in the slave processor
		return 2; //

	} else if(SH->Dval <0) {
		// we locate the minimum saved with this node 
		// retrieve the constrained local min
		//cm=MBCM.GetAt((EUINT) (SH->killers));
		if(cm!=NULL) cm->GetVec(temp_x);
			else {SH->Dval=fabs(SH->Dval); return 0;}
		temp_val = fabs(SH->Dval) - GetConst();
		return 1; // OK
	} else {

	j=0;
	for(i=0;i<Dim;i++) { 
		temp_x[j]=  SVectors[m_lastindex[i]].vec[j]; 
		j++;
	}
	// normalise to simplex

	d = fabs(SH->Dval) - GetConst();

	radius=0;

	for(i=0;i<Dim;i++) {
		temp_x[i] = d/m_Constants[i] - temp_x[i]; 
		if(radius< (frecord-d)/m_Constants[i]) radius = (frecord-d)/m_Constants[i];
	}

	temp_val=d;

#ifdef BOXCONSTRAINTS
// cube constraints
	for(i=0;i<Dim-1;i++) { 
		if(temp_x[i]<0) return 0; 
		if(temp_x[i]>m_bounds[i]) return 0;
	}

#endif

	} // else (normal minimum)
	return 1;
}

void SawToothCover::FormNewSV()
{
	double GVAt;
	GVAt=temp_val; 

// now calculate the actual value of the function
	double val=Value(temp_x);
	GVA=max_(GVA,GVAt);
	flast=val;
	if(frecord > val) { // keep this value if it is the best
		frecord=val;
		xglob_x=temp_x; // copy vector
		xglob_val=val;

		FunCallsWhenFound = MyParent->FunctionCounter;
	}
	
	// form new support vector
	support_vector	sv;

	sv.SVForm(temp_x,val);
	sv.label=LastLabel;
	LastLabel++;

	SVectors.push_back(sv);
}

int		SawToothCover::CheckExistingSV(dVec&  X)
{
	int i;
	dVec Y=X;
	for(i=Dim;i<LastLabel;i++)
	{	SVectors[i].ReturnX(Y);
		if(EquivalentSV(X,Y)) 
		{ 
//			cout << i<< " "<< X;
			return 1;
		}
	}
	return 0;
}

int SawToothCover::EquivalentSV(dVec& X, dVec& Y)
{
	float x,y;
	int i;
	for(i=0;i<Dim;i++)
	{
		x=(float)X[i];
		y=(float)Y[i];
		if(x!=y) return 0;
	}
	return 1;
}






/*------------------------------------------------*/


double SawToothCover::GetConst() {return Const;}
void SawToothCover::SetConst(double newconst) { Const=newconst; }
void SawToothCover::SetIter(int iter) {Iters=iter;}
void SawToothCover::SetCut(double newcut) {Cut=newcut;}
double SawToothCover::GetCut() {return Cut;}
//void	SawToothCover::SetBoxConstraints(dVec& a, dVec& b) {	A=a; B=b; }

void	SawToothCover::SetConstants(dVec& newconst) {m_Constants=newconst;}
void	SawToothCover::SetConstants(double newconst, int dim) {m_Constants.newsize(dim); m_Constants=newconst; m_Constants[dim-1] *=sqrt(dim-1.0);}
// create a vector, the last one *= sqrt(n-1)

void  SawToothCover::SetListSize(int m) {ListSize=m;}	





/* coordinate transformation block -----------* /
dVec	ScaleV;
dVec	ShiftV;
double  TransAlpha;

// linear transformation
dVec y_t, t_t;
real v_t;
dMat Aprod,AprodT;
double VolFactor;
double LipFactor;


void	FormTransformMatricesScale(int Dim, dVec& a, dVec& b)
{
	// Dim is the  dimension of x1..xn, and sum(x_n)<=1 (ie like in cube case, extra variable not there)

	VolFactor=LipFactor=1;

	ScaleV.newsize(Dim);
	ShiftV.newsize(Dim);
	int i;
	for( i=0;i<Dim;i++) { 
			double Uglobal= (2.0-(Dim+1))/ (Dim+1.0*sqrt((double)Dim)) +1;
			ScaleV[i]=(b[i]-a[i])/Uglobal;
			ShiftV[i]=(a[i]);
			VolFactor *=ScaleV[i];
			LipFactor = max_(LipFactor,ScaleV[i]);
		}
	
	y_t.newsize(Dim);  // in cube
	t_t.newsize(Dim+1); // on simplex	
}


void	TransformSystemOfConstraintsScale(dMat& Constr, dVec& RHS, dMat& ModConstr, dVec& ModRHS)
{
	// use the formulas: ModC=Constr * Scale * AprodT ; ModRHS=RHS+ModC*v_t-Con*shift
	int i,j;

	ModConstr = Constr;
	for(i=0;i<ModConstr.dim(1); i++)
		for(j=0;j<ModConstr.dim(2); j++) ModConstr[i][j] *= ScaleV[j];


	ModRHS = RHS;
	for(i=0;i<Constr.dim(1); i++)
		for(j=0;j<Constr.dim(2); j++) ModRHS[i] -= ShiftV[i] * Constr[i][j];

}
*/
