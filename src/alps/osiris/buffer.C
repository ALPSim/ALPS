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

#include <alps/osiris/buffer.h>
#include <alps/osiris/dump.h>
#include <alps/osiris/std/vector.h>

#include <boost/throw_exception.hpp>
#include <stdexcept>

namespace alps {
namespace detail {

// write a few bytes
void Buffer::write_buffer(const void* p, size_type n)
{
  size_type write_pos=size();
  resize(write_pos+n);
  memcpy(&((*this)[write_pos]),p,n);
}


// read a few bytes and update position
void Buffer::read_buffer(void* p, size_type n)
{
  if(read_pos+n>size())
    boost::throw_exception(std::runtime_error("read past buffer"));
  memcpy(p,&((*this)[read_pos]),n);
  read_pos+=n;
}

} // end namespace detail
} // end namespace alps
