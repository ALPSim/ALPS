/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2001-2011 by Prakash Dayal <prakash@comp-phys.org>,
 *                            Matthias Troyer <troyer@comp-phys.org>
 *                            Bela Bauer <bauerb@phys.ethz.ch>
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

/* $Id: ublas.h,v 1.4 2003/09/05 09:27:53 prakash Exp $ */

#ifndef IETL_UBLAS_H
#define IETL_UBLAS_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <ietl/traits.h>

namespace ietl {
    
    template<class T, class Gen> 
    inline void generate(boost::numeric::ublas::vector<T>& c, const Gen& gen) {
        std::generate(c.begin(),c.end(),gen);
    }  
    
    template<class T, class S>
    inline void clear(boost::numeric::ublas::vector<T,S>& c) {
        c.clear();
    }
    
    template<class V>
    typename V::value_type dot(boost::numeric::ublas::vector_expression<V> const & x,
                               boost::numeric::ublas::vector_expression<V> const & y)
    {
        return inner_prod(conj(x), y);
    }
    
    template<class V>
    typename number_traits<typename V::value_type>::magnitude_type
    two_norm(boost::numeric::ublas::vector_expression<V> const & x)
    {
        return norm_2(x);
    }
    
    template < class T>
    void copy(boost::numeric::ublas::vector<T> const & x,
              boost::numeric::ublas::vector<T>& y) {
        y.assign(x);
    }
    
    template<class M, class V>
    void mult(boost::numeric::ublas::matrix_expression<M> const & m,
              boost::numeric::ublas::vector_expression<V> const & x,
              V & y)
    {
        y = prod(m, x);
    }
}

#endif
