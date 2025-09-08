#include "ecam.h"
#ifdef MPI
MemoryBlockMPI<SVSetNode>	MBSV;
#else
MemoryBlock<SVSetNode>	MBSV;
#endif
MemoryBlockE<Children>   MBCL; // blocks for SVSetNode and children array 

MemoryBlockE<SavedIndex>   MBSI; // saved index for queries
MemoryBlockE<ConstrMin>   MBCM; // blocks for ConstMin array

Children* GetMemBlockAt(EUINT b)
{
	return MBCL.GetAt(b);
}

EUINT CreateChildrenArray(EUINT size)
{
	EUINT r;
	Children* block;
	if(size<=ChildrenBlockSize) {
		r= MBCL.GetNextFree(); 
		block=MBCL.GetAt(r);
		block->NextBlock=MB_BADINDEX;
		return r;
	}

	div_t b=div(size-1,ChildrenBlockSize);
	b.quot++; // that's how many blocks are needed
	r=MBCL.GetNextFree();
	block=MBCL.GetAt(r);
	for(int i=1; i<b.quot;i++) {
		block->NextBlock=MBCL.GetNextFree(); 
		block = MBCL.GetAt(block->NextBlock);
	}
	block->NextBlock=MB_BADINDEX;

	return r;
}

void DeleteChildrenArray(EUINT& first)
{
	EUINT r=first;
	Children* block;
	while(r!=MB_BADINDEX) {
		block=MBCL.GetAt(r);
		first=block->NextBlock;	
		MBCL.FreeBlock(r);
		r=first;
	}
	first=MB_BADINDEX;
}
