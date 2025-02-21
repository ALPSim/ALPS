/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2004 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

/* $Id$ */

#ifndef OSIRIS_BUFFER_H
#define OSIRIS_BUFFER_H

#include <alps/config.h>
#include <vector>

namespace alps {
namespace detail {

class Buffer;

/** a simple Buffer class. values can be written into it or read from it. */
class Buffer : public std::vector<char>
{
  public:
  /** create a buffer. */
  Buffer()
    : std::vector<char>(),
      read_pos(0) {}

 /// erase the Buffer
 void clear() {(*this)=Buffer();}

  /** get a pointer to the Buffer. 
      This pointer might be invalidated by writing to the Buffer.
  */
  operator char* () { return (size() ? &(this->operator[](0)) : 0); }

  // write basic data types and arrays of them

template <class T>
  void write(const T* p,size_type n)
  {
    write_buffer(p, n*sizeof(T));
  }

  template <class T>
  void write(const T x) 
  {
    write_buffer(&x, sizeof(T));
  }

  // read basic data types and arrays of them
  
  template <class T>
  void read(T* p,size_type n=1) { read_buffer(p, n*sizeof(T));}

  template <class T>
  void read(T& x) { read_buffer(&x, sizeof(T)); }

private:
  // the position at which reading will take place
  uint32_t read_pos; 

  void write_buffer(const void*, size_type);
  void read_buffer(void*, size_type);
};

} // end namespace detail
} // end namespace alps

#endif // OSIRIS_BUFFER_H
