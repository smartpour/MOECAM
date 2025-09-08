#include <stdlib.h>
#include <malloc.h>
#include <memory.h>

#ifndef MEMORYBLOCK
#define MEMORYBLOCK
/*typedef struct MyStuct{
	double v;
	int i;
} MyStruct_t;
*/

//#define _UNIX64
#ifdef _UNIX_64
#define EUINT unsigned long long int
#else
#define EUINT unsigned int
#endif

#define MB_IDX_SHF 20  // that's how we split the index into 2 parts: block and index within the block
#define MB_IDX_MASK ((1 << (MB_IDX_SHF)) - 1)
#define MB_BLK_SHF (32 - (MB_IDX_SHF))

/***  Macros for calculating the correct location of the node  ***/
#define MB_BLOCK(A) ((A) >> (MB_IDX_SHF))
#define MB_INDEX(A) ((A) & (MB_IDX_MASK))
#define MB_INDEXB(A,B) (((A) << (MB_IDX_SHF)) + B)

/***  The upper limits - 4GB  ***/
#define MB_MAX_NODES 0xFFFFFFFF

/***  Define the ranges  ***/
#define MB_MAX_INDEXES (1 << (MB_IDX_SHF))
#define MB_MAX_BLOCKS (((MB_MAX_NODES) / (MB_MAX_INDEXES) + 1))

#define MB_BADINDEX		0xFFFFFFFF
#define MB_SPECIALINDEX 0xEFFFFFFF
//#define MB_BLOCKSIZE	0x7FFFF
#define MB_BLOCKSIZE	0xFFFFF  //??????? why not
//#define MB_BLOCKSIZE	0xFF

#define	HasFree(B) ( !((B) & 0x1) )
#define	HasAnyFree(B) ( (B!=0xFFFFFFFE) )
#define SetFree(B) ( ((B) &= 0xFFFFFFFE) )

#define SetFreeI(B,r)   { (B) &= (~(0x1 << r)) ;(B) &= 0xFFFFFFFE; } 


inline int  SetOccupied(EUINT &B, short i) { 
	B |= (0x1 << i); 
	if(HasAnyFree(B)) {SetFree(B); return 0; } 
	else  {B |= 0x00000001; return 1;}
};

inline short WhichFree(EUINT B) {
	for(short register i=1;i<32;i++)
		if(!((B>>i) & 0x1) ) return i;
	return 0;
};

template <class T>
class MemBlock {
public:

    typedef         T   value_type;
    typedef         T*  pointer;
    typedef         T&  reference;
    typedef const   T&  const_reference;


	value_type * m_data;
//	EUINT m_index[1024];
	EUINT	m_NextAvail,m_temp;

	short i,j,k;
	//0 means free block, 1 means occupied

	MemBlock() {       // 7FFF= 31*32*32 +31*32 +31
		m_data=(value_type*) calloc(MB_BLOCKSIZE,sizeof(T) ); // 31000 blocks of size T, >64kb
//		for(short  i=0;i<1024;i++) m_index[i]=0;
		memset(m_data,0xFF,MB_BLOCKSIZE*sizeof(T));
		m_NextAvail=0;
	};

	~MemBlock() {
		free(m_data);
	};

	EUINT	GetNextFree() {
		if(m_NextAvail>=MB_BLOCKSIZE) return MB_BADINDEX;
		m_temp=m_NextAvail;
		m_NextAvail++;
		return m_temp;
	};

	EUINT	GetNextFree(int M) {
		if(m_NextAvail+M>=MB_BLOCKSIZE) return MB_BADINDEX;
		m_temp=m_NextAvail;
		m_NextAvail+=M;
		return m_temp;
	};


	inline int	IsFree() { return (m_NextAvail<MB_BLOCKSIZE); };
	inline int	IsFreeM(int M) { return (m_NextAvail+M < MB_BLOCKSIZE); };
	inline void	FreeBlock(EUINT B) { memset(m_data+B,0xFF,sizeof(T)); };

	T* GetAt(EUINT B) { return (T*) (m_data+B); };

	void SetAt(EUINT B, T* Value) { memcpy(m_data + B, Value, sizeof(T)); };

	inline reference operator()(EUINT B) { return *((T*)(GetAt(B)));	};
    inline const_reference operator() (EUINT B) const { *((T*)(GetAt(B))); };

	void ClearAll() {	m_NextAvail=0;};
};




template <class T>
class MemoryBlock 
{
public:
    typedef         T   value_type;
    typedef         T*  pointer;
    typedef         T&  reference;
    typedef const   T&  const_reference;

  MemBlock<T>** block;

  EUINT nodeCount, emptyBlocks, currentBlock;
  int valid;
  int LimitBlocks;

	MemoryBlock(void)
	{
		block =  (MemBlock<T>** ) calloc(MB_MAX_BLOCKS, sizeof(MemBlock<T>*));
//		assert(block != NULL);
		nodeCount = emptyBlocks = 0;
/***  The starting point is 0 but set to -1 because the
      _createNextBlock will increament the value
      before using it
***/
		currentBlock = (EUINT) -1;
		_createNextBlock();
		valid=1;
		LimitBlocks=MB_MAX_BLOCKS;
	};

	~MemoryBlock(void)
	{
		for(EUINT  loop = currentBlock + emptyBlocks; loop > 0; loop--) delete (block[loop]);
/***  To free the first block!  ***/
	  delete(block[0]);
	  free(block);
	  valid=0;
	};

	EUINT	GetNextFree() {
		nodeCount++;
		EUINT loop;
		for(loop=0; loop <= currentBlock; loop++)
			if(block[loop]->IsFree()) {
				loop = MB_INDEXB(loop, block[loop]->GetNextFree());
				return loop;
			}
		// no space left
		if(currentBlock < (EUINT)(LimitBlocks-1) /*MB_MAX_BLOCKS-2*/) {
			_createNextBlock();
			loop = MB_INDEXB(loop, block[currentBlock]->GetNextFree());
			return loop;
		}
		nodeCount--;
//		exit(20);
		return MB_BADINDEX;
	};

	EUINT	GetNextFree(int M) {
		nodeCount+=M;
		EUINT loop;
		for(loop=0; loop <= currentBlock; loop++)
			if(block[loop]->IsFreeM(M)) {
				loop = MB_INDEXB(loop, block[loop]->GetNextFree(M));
				return loop;
			}
		// no space left
		if(currentBlock < (EUINT)(LimitBlocks-1 )/*MB_MAX_BLOCKS-2*/) {
			_createNextBlock();
			loop = MB_INDEXB(loop, block[currentBlock]->GetNextFree(M));
			return loop;
		}
		nodeCount-=M;
//		exit(20);
		return MB_BADINDEX;
	};

	inline void	FreeBlock(EUINT B)	{
		nodeCount--;
		block[MB_BLOCK(B)]->FreeBlock(MB_INDEX(B));
	//	B=MB_BADINDEX;
	};

	inline int	IsFree() { 
		if(currentBlock < MB_MAX_BLOCKS-1) return 1;
		for(EUINT  loop=0; loop <= currentBlock; loop++)
			if(block[loop]->IsFree()) return 1;
		return 0; 
	};

	inline T* GetAt(EUINT B, T* Value=NULL) { return block[MB_BLOCK(B)]->GetAt(MB_INDEX(B)); };

	inline T* GetAddress(EUINT B, T* Value=NULL) {return block[MB_BLOCK(B)]->GetAt(MB_INDEX(B)); }

	inline void SetAt(EUINT B, T* Value) { block[MB_BLOCK(B)]->SetAt(MB_INDEX(B),Value); };

	inline reference operator()(EUINT B) { return *((T*)(GetAt(B)));	};
    inline const_reference operator() (EUINT B) const { *((T*)(GetAt(B))); };

	inline void _createNextBlock()
	{
	  currentBlock++;

	  if(emptyBlocks == 0)
	  {
		block[currentBlock] = new MemBlock<T>;
//		assert(block[currentBlock] != NULL);
	  }
	  else
		emptyBlocks--;
	}

	void Discard(EUINT B, T* Value) {}  // does nothing

	void ClearAll() 
	{
//		cout << "commited blocks " <<currentBlock + emptyBlocks<<" of size "<<block[0]->m_NextAvail <<endl;

		for(EUINT  loop = currentBlock + emptyBlocks; loop > 0; loop--) block[loop]->ClearAll();
/***  To free the first block!  ***/
	  block[0]->ClearAll();
	};

	inline int IsValid() {return valid;}
};




/*-------------old version for lists------------------------*/
template <class T>
class MemBlockE {
public:

    typedef         T   value_type;
    typedef         T*  pointer;
    typedef         T&  reference;
    typedef const   T&  const_reference;


	value_type * m_data;
	EUINT m_index[1024];

	short i,j,k;
	//0 means free block, 1 means occupied

	MemBlockE() {       // 7FFF= 31*32*32 +31*32 +31
		m_data=(value_type*) calloc(0x7FFF,sizeof(T) ); // 31000 blocks of size T, >64kb
		for(short register i=0;i<1024;i++) m_index[i]=0;
		memset(m_data,0xFF,0x7FFF*sizeof(T));
	};

	~MemBlockE() {
		free(m_data);
	};

	EUINT	GetNextFree() {
		//short register i,j,k;
		i=WhichFree(m_index[0]);
		j=WhichFree(m_index[i]);
		k=WhichFree(m_index[i*32+j]);
		if(SetOccupied(m_index[i*32+j],k))
			if(SetOccupied(m_index[i],j))
				SetOccupied(m_index[0],i);
		return GetAddress();// i,j,k);
	};

	void	FreeBlock(EUINT B)	{
		//short i,j,k;
		GetIJK(B);//,i,j,k);
		SetFreeI(m_index[i*32+j],k); 
		SetFreeI(m_index[i],j); 
		SetFreeI(m_index[0],i);
	};

	EUINT	GetAddress(){ //short i, short j, short k)	{
		EUINT r=(i-1);
		r *= 1024;
		r = r+ (j-1)*32 + k-1;
		//return ((i-1)*32*32+(j-1)*32 +k-1 ); //sizeof(MyStruct_t)*
		return r;
	};

	void	GetIJK(EUINT B)//, short& i, short& j, short &k)
	{
		div_t t=div(B,32); ///sizeof(MyStruct_t)
		j=t.quot; k=t.rem+1;
		t=div(j,32);
		j=t.rem+1; i=t.quot+1;
	};

	inline int	IsFree() { return HasFree(m_index[0]); };

	T* GetAt(EUINT B) { return (T*) (m_data+B); };

	void SetAt(EUINT B, T* Value) { memcpy(m_data + B, Value, sizeof(T)); };

	inline reference operator()(EUINT B) { return *((T*)(GetAt(B)));	};
    inline const_reference operator() (EUINT B) const { *((T*)(GetAt(B))); };

	void ClearAll() {
		for(short register i=0;i<1024;i++) m_index[i]=0;
	};
};




template <class T>
class MemoryBlockE 
{
public:
    typedef         T   value_type;
    typedef         T*  pointer;
    typedef         T&  reference;
    typedef const   T&  const_reference;

  MemBlockE<T>** block;

  EUINT nodeCount, emptyBlocks, currentBlock;
  int valid;

	MemoryBlockE(void)
	{
		block =  (MemBlockE<T>** ) calloc(MB_MAX_BLOCKS, sizeof(MemBlockE<T>*));
//		assert(block != NULL);
		nodeCount = emptyBlocks = 0;
/***  The starting point is 0 but set to -1 because the
      _createNextBlock will increament the value
      before using it
***/
		currentBlock = (EUINT) -1;
		_createNextBlock();
		valid=1;
	};

	~MemoryBlockE(void)
	{
		for(EUINT register loop = currentBlock + emptyBlocks; loop > 0; loop--) delete (block[loop]);
/***  To free the first block!  ***/
	  delete(block[0]);
	  free(block);
	  valid=0;
	};

	EUINT	GetNextFree() {
		nodeCount++;
		EUINT loop;
		for(loop=0; loop <= currentBlock; loop++)
			if(block[loop]->IsFree()) {
				loop = MB_INDEXB(loop, block[loop]->GetNextFree());
				return loop;
			}
		// no space left
		if(currentBlock < MB_MAX_BLOCKS-2) {
			_createNextBlock();
			loop = MB_INDEXB(loop, block[currentBlock]->GetNextFree());
			return loop;
		}
		nodeCount--;
		exit(20);
		return MB_BADINDEX;
	};

	inline void	FreeBlock(EUINT& B)	{
		nodeCount--;
		block[MB_BLOCK(B)]->FreeBlock(MB_INDEX(B));
		B=MB_BADINDEX;
		//if(BLOCK(B) == currentBlock  && 
	};
	inline void	FreeBlockC(EUINT B)	{
		nodeCount--;
		block[MB_BLOCK(B)]->FreeBlock(MB_INDEX(B));
	};


	inline int	IsFree() { 
		if(currentBlock < MB_MAX_BLOCKS-1) return 1;
		for(EUINT register loop=0; loop <= currentBlock; loop++)
			if(block[loop]->IsFree()) return 1;
		return 0; 
	};

	inline T* GetAt(EUINT B) { return block[MB_BLOCK(B)]->GetAt(MB_INDEX(B)); };

	inline void SetAt(EUINT B, T* Value) { block[MB_BLOCK(B)]->SetAt(MB_INDEX(B),Value); };

	inline reference operator()(EUINT B) { return *((T*)(GetAt(B)));	};
    inline const_reference operator() (EUINT B) const { *((T*)(GetAt(B))); };

	inline void _createNextBlock()
	{
	  currentBlock++;

	  if(emptyBlocks == 0)
	  {
		block[currentBlock] = new MemBlockE<T>;
//		assert(block[currentBlock] != NULL);
	  }
	  else
		emptyBlocks--;
	}

	void ClearAll() 
	{
	  for(EUINT register loop = currentBlock + emptyBlocks; loop > 0; loop--) block[loop]->ClearAll();
/***  To free the first block!  ***/
	  block[0]->ClearAll();
	};

	inline int IsValid() {return valid;}
};



#endif

