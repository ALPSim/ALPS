/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
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

#ifndef __GLOBALS_H__
#define __GLOBALS_H__

#ifndef MIN_VECTOR_SIZE
#define MIN_VECTOR_SIZE 1000
#endif
#ifndef MIN_MATRIX_SIZE
#define MIN_MATRIX_SIZE 100
#endif

#include "vector.h"
#include "matrix.h"

namespace dmtk
{

template<class T>
class DMTKglobals
{
  public:
    Matrix<T> m1;
    Matrix<T> m2;
    Matrix<T> m3;
    Matrix<T> m4;
    Vector<T> v1;
    Vector<T> v2;
    Vector<T> v3;
    Vector<T> v4;

    DMTKglobals() { m1 = m2 = m3 = m4 = Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE); v1 = v2 = v3 = v4 = Vector<T>(MIN_VECTOR_SIZE); }
};

static DMTKglobals<double> globals_double;
static DMTKglobals<complex<double> > globals_complex;

DMTKglobals<double> &
get_globals(const double&)
{
  return globals_double;
}

DMTKglobals<complex<double> >&
get_globals(const complex<double>&)
{
  return globals_complex;
}

} //namespace dmtk

#endif // __GLOBALS_H__
