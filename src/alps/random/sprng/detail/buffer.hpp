/* 
 * Copyright Matthias Troyer 2005
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
 */

#ifndef ALPS_RANDOM_SPRNG_DETAIL_BUFFER_HPP
#define ALPS_RANDOM_SPRNG_DETAIL_BUFFER_HPP

#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include <vector>
#include <boost/throw_exception.hpp>
#include <boost/assert.hpp>

namespace alps { namespace random { namespace sprng { namespace detail {

/// A buffer class to pack and unpack SPRNG generators
///
/// The SPRNG libraries can pack and unpack random number generators into and from
/// an array of char. This class is wrapper around these packing and unpacking
/// functions of the SPRNG libraries and has convenience functions to write
/// and read the buffer contents as a sequence of integers

class buffer
{
public:

  typedef std::vector<char> buffer_type;

  /// creates an empty buffer  
  buffer()
  {}

  /// packs a SPRNG generator into a new buffer
  /// @param genptr the SPRNG pointer to the generator
  /// @param pack_fun the SPRNG function used to pack the generator
  buffer(int* genptr, int(*pack_fun)(int*,char**) )
  {
    pack(genptr,pack_fun);
  }
  
  /// packs a SPRNG generator into an existing buffer.
  ///
  /// The previous contents of the buffer is erased.
  ///
  /// @param genptr the SPRNG pointer to the generator
  /// @param pack_fun the SPRNG function used to pack the generator
  
  void pack(int* genptr, int(*pack_fun)(int*,char**) )
  {
    buffer_.clear();
    char* buf;
    std::size_t len = (*pack_fun)(genptr,&buf);
    if (buf==0)
      boost::throw_exception(std::bad_alloc());
    std::copy(buf,buf+len,std::back_inserter(buffer_));
  }

  /// unpacks a SPRNG generator from a buffer.
  /// 
  /// @param unpack_fun the SPRNG function used to unpack the generator
  /// @returns a newly created SPRNG generator. The return value is 0 if the creation of the generator failed
  int* unpack(int*(*unpack_fun)(const char*) )
  {
    BOOST_ASSERT(!buffer_.empty());
    int* genptr = (*unpack_fun)(&buffer_[0]);
    if (genptr==0)
      boost::throw_exception(std::bad_alloc());
    return genptr;
  }
  
  /// writes the length and contents of the buffer into an output stream
  template<class OutputStream>
  void write(OutputStream& os) const
  {
    os << buffer_.size();
    for (buffer_type::const_iterator it=buffer_.begin();it!=buffer_.end();++it)
      os << " " << static_cast<int>(*it);
  }

  /// read the length and contents of the buffer from an input stream
  template<class InputStream>
  void read(InputStream& is)
  {
    std::size_t len;
    is >> len;
    buffer_.resize(len);
    for (buffer_type::iterator it=buffer_.begin();it!=buffer_.end();++it) 
    {
      int x;
      is >> x;
      *it =x;
    }
  }

  bool operator==(const buffer& b)
  {
    return b.buffer_ == buffer_;
  }  
private:
  buffer_type buffer_;
};



} } } } // namespace alps::random::sprng::detail

#endif // ALPS_RANDOM_SPRNG_DETAIL_BUFFER_HPP
