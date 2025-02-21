/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2002 by Rene Villiger <rvilliger@smile.ch>,
*                            Prakash Dayal <prakash@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
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

/* $Id: fmatrix.h,v 1.8 2003/04/24 13:49:40 renev Exp $ */

#ifndef IETL_FMATRIX_H
#define IETL_FMATRIX_H

#include <cstddef>

#undef minor
#undef data

namespace ietl
{

    template <class T>
    class FortranMatrix 
    {
      private:
         T* p;
         std::size_t n_;
         std::size_t m_;
      public:
         typedef std::size_t size_type;
         FortranMatrix(size_type n, size_type m) : n_(n), m_(m) { p = new T[m*n]; };
         ~FortranMatrix() { delete[] p; };
         T* data() { return p; };
         const T* data() const { return p; };
         T operator()(size_type i, size_type j) const { return p[i+j*n_]; };
         T& operator()(size_type i, size_type j) { return p[i+j*n_]; };
         void resize(size_type n, size_type m) { m_=m; n_=n; delete[] p; p = new T[m*n]; };
         size_type nrows() { return n_; };
         size_type ncols() { return m_; };
         size_type minor() { return n_; };
    };  
}
#endif

