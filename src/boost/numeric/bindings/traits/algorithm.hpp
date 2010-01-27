/*
 * 
 * Copyright (c) 2002, 2003 Kresimir Fresl, Toon Knapen and Karl Meerbergen
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * KF acknowledges the support of the Faculty of Civil Engineering, 
 * University of Zagreb, Croatia.
 *
 */

#ifndef boost_numeric_bindings_traits_algorithm_hpp
#define boost_numeric_bindings_traits_algorithm_hpp

#include <boost/numeric/bindings/traits/type_traits.hpp>

namespace boost { namespace numeric { namespace bindings { namespace traits {

  ///
  /// To be used instead of operator== for numeric types
  /// Implemented as functor instead of free function because of specialisation
  /// rationale: operator== on builtin types can not be overloaded.

  template < typename T >
  struct is_equal
  {
    is_equal(typename type_traits< T >::real_type tolerance) : tolerance_( tolerance ) {} 

    bool operator()(const T& a, const T& b)
    { return std::abs( a - b ) < tolerance_ ; }

    // bool operator()(const T& a, const T& b, typename value_traits< T >::value_type tolerance) 
    // { return std::abs( a - b ) < tolerance ; }

    typename type_traits< T >::real_type tolerance_ ;
  } ;

}}}}

#endif // boost_numeric_bindings_traits_algorithm_hpp
