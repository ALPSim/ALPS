/*****************************************************************************
 *
 * ALPS DMFT Project - BLAS Compatibility headers
 *  BLAS headers for accessing BLAS from C++.
 *
 * Copyright (C) 2010 Matthias Troyer <gtroyer@ethz.ch>
 *
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

#include <boost/numeric/bindings/blas.hpp>
#include<vector>

// provide overloads for types where blas can be used        

namespace blas{

#define IMPLEMENT_FOR_REAL_BLAS_TYPES(F) F(float) F(double)

#define IMPLEMENT_FOR_COMPLEX_BLAS_TYPES(F) \
F(std::complex<float>) \
F(std::complex<double>)

#define IMPLEMENT_FOR_ALL_BLAS_TYPES(F) \
IMPLEMENT_FOR_REAL_BLAS_TYPES(F) \
IMPLEMENT_FOR_COMPLEX_BLAS_TYPES(F) 
} // namespace
