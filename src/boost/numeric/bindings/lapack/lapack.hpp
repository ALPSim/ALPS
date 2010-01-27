/*
 * 
 * Copyright (c) Toon Knapen & Kresimir Fresl 2003
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * KF acknowledges the support of the Faculty of Civil Engineering, 
 * University of Zagreb, Croatia.
 *
 */

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_LAPACK_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_LAPACK_HPP

// linear systems

#include <boost/numeric/bindings/lapack/gesv.hpp>
#include <boost/numeric/bindings/lapack/posv.hpp>
#include <boost/numeric/bindings/lapack/ppsv.hpp>
#include <boost/numeric/bindings/lapack/sysv.hpp>
#include <boost/numeric/bindings/lapack/spsv.hpp>
#include <boost/numeric/bindings/lapack/hesv.hpp>
#include <boost/numeric/bindings/lapack/hpsv.hpp>

// eigenproblems 

#include <boost/numeric/bindings/lapack/gees.hpp>
#include <boost/numeric/bindings/lapack/trevc.hpp>
#include <boost/numeric/bindings/lapack/trexc.hpp>
#include <boost/numeric/bindings/lapack/hbev.hpp>

// SVD

#include <boost/numeric/bindings/lapack/gesvd.hpp>
#include <boost/numeric/bindings/lapack/gesdd.hpp>

// Miscellaneous

#include <boost/numeric/bindings/lapack/geqrf.hpp>
#include <boost/numeric/bindings/lapack/ormqr.hpp>



#endif // BOOST_NUMERIC_BINDINGS_LAPACK_LAPACK_HPP

