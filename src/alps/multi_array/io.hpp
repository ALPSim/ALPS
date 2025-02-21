/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2012 by Ilia Zintchenko <iliazin@gmail.com>                       *
 *                       Jan Gukelberger                                           *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_MULTI_ARRAY_IO_HPP
#define ALPS_MULTI_ARRAY_IO_HPP

#include <alps/multi_array/multi_array.hpp>

namespace alps{

  template <class T, std::size_t D, class Allocator>
  std::ostream& operator<<(std::ostream& out, const multi_array<T,D,Allocator>& a)
  {
    boost::array<typename boost::multi_array<T,D,Allocator>::index, D> index;
    T* rE = const_cast< T* >(a.data());

    for(std::size_t dir = 0; dir < D; dir++) out << "{";
      
    for(int i = 0; i < a.num_elements(); i++)
      {
  	for(std::size_t dir = 0; dir < D; dir++ )
  	  index[dir] = (rE - a.origin()) / a.strides()[dir] % a.shape()[dir] +  a.index_bases()[dir];

	if(index[D-1] == a.shape()[D-1]-1){

	  int M = 1;
	  for(std::size_t dir = D-1; dir > 0; --dir){
	    if(index[dir-1] == a.shape()[dir-1]-1) M++;
	    else break;
	  }

	  out << a(index);

	  for(int m = 0; m < M; ++m)
	    out << "}";
	
	  if(M < D){
	    
	    if(M > 1 || M == D - 1)
	      out << ",\n";
	    else
	      out << ", ";  

	    for(int m = 0; m < M; ++m)
	      out << "{";  

	  }

	}
	else 
	  out << a(index) << ", ";

  	++rE;
      }

    out << ";";

    return out;
  }

}//namespace alps

#endif // ALPS_MULTI_ARRAY_IO_HPP
