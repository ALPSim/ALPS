// NAME
//   HeapWrapper
// PURPOSE
//   debugging versions of ::new and ::delete with necessary support
//   routines.  
// NOTES
//   Author: Stephen R. Davis, "C++ Programmers Companion", pp 343-351
//

#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;

#include "heapwrapper.h"

// Wrapper - the heap wrapper structure which we will
//           use to hold the prologue and epilogue information
struct Wrapper {
  // prologue information
  Wrapper* pNext;
  int length;
  char     prologue;
  // user data
  double data[1];		// type double to make sure we
				// are on a word boundary
  // epilogue follows
};

static Wrapper *pFirst = 0;	// pointer to linked list of allocated blocks

static size_t heap_size = 0;	// to keep track of the size of the heap

// prototype declarations
static void werror (Wrapper*, char*, int = 0);
void* allocateFn (size_t ,  char *, int );
void deleteFn (void*);

// werror - display heap wrapper error message and quit
static void werror (Wrapper *,
		    char *pErrorString, int )
{
  cerr << "Heap error: " << pErrorString << std::endl;
  exit (1);
}

// new - allocate a block from the heap and then tag on
//       the debugging type, data information; return
//       to the user the address of his data

void* operator new (size_t sizeBlock) throw (std::bad_alloc)
{
  return allocateFn (sizeBlock, "Unknown", 0);
}

void* operator new [] (size_t sizeBlock) throw(std::bad_alloc)
{
  return allocateFn (sizeBlock, "Unknown", 0);
}

void* allocateFn (size_t sizeBlock,
		  char *, int )
{
  Wrapper *pBlock;
  char    *pUserData;

  if (sizeBlock == 0)
    return (void*)0;
  else
    {
      pBlock = (Wrapper*)malloc (sizeBlock + sizeof(Wrapper));
      if (!pBlock)
	{
	  return (void*)0;
	}

      // now fill in the data
      pBlock->length = sizeBlock;
      heap_size += sizeBlock;

      // fill in the prologue and epilogue sections
      pBlock->prologue = 0x12;
      pUserData = (char*)pBlock->data;
      char *pEpilogue = pUserData + sizeBlock;
      *pEpilogue = 0x21;

      // and add it to the list
      pBlock->pNext = pFirst;
      pFirst = pBlock;
      // now return the user a pointer to his data
      return (void*) pUserData;
    }
}

// delete - accept either form of delete function;
//          size information unnecessary
void operator delete (void *pUserData) throw()
{
  if (pUserData)
    deleteFn(pUserData);
}

void operator delete (void *pUserData, size_t) 
{
   deleteFn(pUserData);
}

void operator delete [] (void *pUserData) throw()
{
  if (pUserData)
    deleteFn(pUserData);
}

void operator delete [] (void *pUserData, size_t)
{
  deleteFn(pUserData);
}

// deleteFn - check out the heap block to make sure all
//            is kosher then remove it from the list
void deleteFn (void *pUserData)
{
  // calculate the address of the original Wrapper block
  long long offset = (long long)&((Wrapper*)0)->data;
  // (offset is now the number of bytes the field
  // data is from the beginning of a Wrapper block)
  Wrapper *pBlock = (Wrapper*)((char*)pUserData - offset);

  // check the prologue and epilogue
  if (pBlock->prologue != 0x12)
    werror (pBlock, "Prologue overwritten", 1);
  char *pEpilogue = (char*)pUserData + pBlock->length;
  if (*pEpilogue != 0x21)
    werror (pBlock, "Epilogue overwritten", 1);

  // now unlink it from the list
  Wrapper *pWrapper = pFirst;
  int foundIt = 0;
  if (pWrapper == pBlock)
    {
      pFirst = pBlock->pNext;
      foundIt = 1;
    }
  else
    {
      while (pWrapper)
	{
	  if (pWrapper->pNext == pBlock)
	    {
	      pWrapper->pNext = pBlock->pNext;
	      foundIt = 1;
	      break;
	    }
	  pWrapper = pWrapper->pNext;
	}
    }
  if (!foundIt)
    {
	werror (pBlock, "Block not in list; released twice?");
    }

  // now free the block
  heap_size -= pBlock->length;
  if (foundIt)
    free ((void*)pBlock);
}


size_t size_of_heap ()
{
  return heap_size;
}
