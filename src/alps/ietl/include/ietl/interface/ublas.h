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
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
