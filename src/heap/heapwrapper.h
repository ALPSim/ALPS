// NAME
//   HeapWrapper
// PURPOSE
//   header stuff for the Heap Wrapper by Stephen R. Davis
// NOTES
//   Author: Stephen R. Davis, "C++ Programmers Companion"
//   10-Nov-1993  <heeb@phys.ethz.ch>
//   Created.
//

#ifndef HEAPWRAPPER_H__
#define HEAPWRAPPER_H__

#include <cstdlib>

void* operator new (size_t) throw(std::bad_alloc);
void operator delete (void*) throw();
void* operator new [] (size_t) throw(std::bad_alloc);
void operator delete [] (void*) throw();
size_t size_of_heap ();

#endif
