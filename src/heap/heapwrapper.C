// NAME
//   HeapWrapper
// PURPOSE
//   debugging versions of ::new and ::delete with necessary support
//   routines.  
// NOTES
//   Author: Stephen R. Davis, "C++ Programmers Companion", pp 343-351
//   
//   10-Nov-1993  <heeb@phys.ethz.ch>
//   Created.
//


static const char HeapWrapper_cc_RCS_ID[] =
"$Id$";

//  $Log$
//  Revision 1.1  2004/09/03 14:24:32  troyer
//  *** empty log message ***
//
// Revision 1.1  1993/12/23  18:54:31  heeb
// Initial revision
//

#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;

#include "heapwrapper.h"

#ifdef USE_HEAP_WRAPPER

#define FILENAME_LENGTH 32

// Wrapper - the heap wrapper structure which we will
//           use to hold the prologue and epilogue information
struct Wrapper {
  // prologue information
  Wrapper* pNext;
  int length;
  char fileName[FILENAME_LENGTH+1];
  unsigned lineNumber;
  char     prologue;
  // user data
  double data[1];		// type double to make sure we
				// are on a word boundary
  // epilogue follows
};

static Wrapper *pFirst = 0;	// pointer to linked list of
				// allocated blocks

static size_t heap_size = 0;	// to keep track of the size of the heap

// prototype declarations
static void werror (Wrapper*, char*, int = 0);
// both new operators
void* operator new (size_t) throw (std::bad_alloc);
void* operator new (size_t, char*, int);
void* allocateFn (size_t sizeBlock,
		  char *pFile, int lineNo);

// and both delete operators
void operator delete (void*) throw();
//void operator delete (void*, size_t);
void deleteFn (void*);
void displayAllocated ();

// werror - display heap wrapper error message and quit
static void werror (Wrapper *pBlock,
		    char *pErrorString, int fatal)
{
  cerr << "Heap error: " << pErrorString << std::endl;
//    << , Allocated: "
//       << pBlock->fileName << ", " << pBlock->lineNumber << endl;
//  if (fatal)
    exit (1);
};

// new - allocate a block from the heap and then tag on
//       the debugging type, data information; return
//       to the user the address of his data
void* operator new (size_t sizeBlock,
		      char *pFile, int lineNo)
{
  return allocateFn (sizeBlock, pFile, lineNo);
};

void* operator new (size_t sizeBlock) throw (std::bad_alloc)
{
  return allocateFn (sizeBlock, "Unknown", 0);
};

#ifndef NO_VECTOR_NEW
void* operator new [] (size_t sizeBlock,
		      char *pFile, int lineNo)
{
  return allocateFn (sizeBlock, pFile, lineNo);
};

void* operator new [] (size_t sizeBlock) throw(std::bad_alloc)
{
  return allocateFn (sizeBlock, "Unknown", 0);
};
#endif

void* allocateFn (size_t sizeBlock,
		  char *pFile, int lineNo)
{
  Wrapper *pBlock;
  char    *pUserData;

  if (sizeBlock == 0)
    return (void*)0;
  else
    {
      if (sizeBlock == 8)
	{
	  //int waitfordebugger = 1;
//	  while (waitfordebugger);
	}
      // get a block from the heap
      pBlock = (Wrapper*)malloc (sizeBlock + sizeof(Wrapper));
      if (!pBlock)
	{
	  return (void*)0;
	}

      // now fill in the data
      strncpy (pBlock->fileName, pFile, FILENAME_LENGTH);
      pBlock->fileName[FILENAME_LENGTH] = '\0';
      pBlock->lineNumber = lineNo;
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
};

// delete - accept either form of delete function;
//          size information unnecessary
void operator delete (void *pUserData) throw()
{
  if (pUserData)
    deleteFn(pUserData);
};

void operator delete (void *pUserData, size_t size) 
{
   deleteFn(pUserData);
}

#ifndef NO_VECTOR_NEW
void operator delete [] (void *pUserData) throw()
{
  if (pUserData)
    deleteFn(pUserData);
};

void operator delete [] (void *pUserData, size_t size)
{
  deleteFn(pUserData);
}
#endif

// deleteFn - check out the heap block to make sure all
//            is kosher then remove it from the list
void deleteFn (void *pUserData)
{
  // calculate the address of the original Wrapper block
  int offset = (int)&((Wrapper*)0)->data;
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
};

// displayAllocated - this function simply displays all
//                    of the heap blocks that are still
//                    "check out"
void displayAllocated()
{
  Wrapper *pBlock;
  pBlock = pFirst;
  while (pBlock)
    {
      cerr << pBlock->length << " bytes allocated at "
	   << pBlock->fileName << ";" << pBlock->lineNumber << endl;
      pBlock = pBlock->pNext;
    }
};

size_t size_of_heap ()
{
  return heap_size;
};

#endif
