// NAME
//   HeapWrapper
// PURPOSE
//   header stuff for the Heap Wrapper by Stephen R. Davis
// NOTES
//   Author: Stephen R. Davis, "C++ Programmers Companion"
//   10-Nov-1993  <heeb@phys.ethz.ch>
//   Created.
//

// define the following for a compiler which does not support
// ``operator new []'' and ``operator delete []''
// GCC started supporting this from version 2.6.x
//#define NO_VECTOR_NEW

#ifndef HEAPWRAPPER_H__
#define HEAPWRAPPER_H__


/*
 *  $Log$
 *  Revision 1.1  2004/09/03 14:24:32  troyer
 *  *** empty log message ***
 *
 * Revision 1.1  1993/12/23  18:54:31  heeb
 * Initial revision
 *
 */

#define USE_HEAP_WRAPPER

#include <cstdlib>

#ifdef USE_HEAP_WRAPPER

#define NEW new(__FILE__, __LINE__)
void* operator new (size_t, char*, int);
void* operator new (size_t) throw(std::bad_alloc);
void operator delete (void*) throw();
//void operator delete (void*, size_t);
#ifndef NO_VECTOR_NEW
void* operator new [] (size_t, char*, int);
void* operator new [] (size_t) throw(std::bad_alloc);
void operator delete [] (void*) throw();
//void operator delete [] (void*, size_t);
#endif
void displayAllocated ();
size_t size_of_heap ();

#else

#define NEW new

#endif

#define DELETE(x)  ((delete (x)), (x) = 0)


#endif
