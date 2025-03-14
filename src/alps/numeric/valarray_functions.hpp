/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

/* $Id: valarray_functions.h 3520 2009-12-11 16:49:53Z tamama $ */

#ifndef ALPS_NUMERIC_VALARRAY_FUNCTIONS_HPP
#define ALPS_NUMERIC_VALARRAY_FUNCTIONS_HPP



namespace alps { 
  namespace numeric {

    template <class T>
    std::ostream& operator<< (std::ostream &out, std::valarray<T> const & val)
    {
      std::copy(&const_cast<std::valarray<T>&>(val)[0],&const_cast<std::valarray<T>&>(val)[0]+val.size(),std::ostream_iterator<T>(out,"\t"));
      return out;
    }


  }
}

#endif // ALPS_NUMERIC_VALARRAY_FUNCTIONS_HPP




