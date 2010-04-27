//
//  Copyright (C) 2002, 2003 Si-Lab b.v.b.a., Toon Knapen and Kresimir Fresl
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_DETAIL_CONFIG_FORTRAN_HPP
#define BOOST_NUMERIC_BINDINGS_DETAIL_CONFIG_FORTRAN_HPP

#if defined(BIND_FORTRAN_LOWERCASE_UNDERSCORE) || defined(BIND_FORTRAN_LOWERCASE)
// Allow manual override of the defaults, e.g. if you want to use a fortran
// lib compiled with gcc from MSVC
#else

// First we need to know what the conventions for linking
// C with Fortran is on this platform/toolset
#if defined(__GNUC__) || defined(__ICC) || defined(__sgi) || defined(__COMO__) || defined(__KCC)
#define BIND_FORTRAN_LOWERCASE_UNDERSCORE
#elif defined(__IBMCPP__) || defined(_MSC_VER)
#define BIND_FORTRAN_LOWERCASE
#else
#error do not know how to link with fortran for the given platform
#endif

#endif

// Next we define macro's to convert our symbols to
// the current convention
#if defined(BIND_FORTRAN_LOWERCASE_UNDERSCORE)
#define FORTRAN_ID( id ) id##_
#elif defined(BIND_FORTRAN_LOWERCASE)
#define FORTRAN_ID( id ) id
#else
#error do not know how to bind to fortran calling convention
#endif


// "g77" or clapack or "gfortran -ff2c"
//#define BIND_FORTRAN_F2C_RETURN_CONVENTIONS
// "g77 -fno-f2c" or "gfortran"
//#define BIND_FORTRAN_NO_F2C_RETURN_CONVENTIONS

// As long as f2c return conventions are the common case,
// we turn them on unless requested otherwise
#ifndef BIND_FORTRAN_NO_F2C_RETURN_CONVENTIONS
#define BIND_FORTRAN_F2C_RETURN_CONVENTIONS
#endif

#ifndef BIND_FORTRAN_INTEGER_8
    typedef int fortran_int_t;
#else
    typedef std::ptrdiff_t fortran_int_t;
#endif

#endif // BOOST_NUMERIC_BINDINGS_TRAITS_FORTRAN_H
