/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1999-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_NUMERIC_SCALAR_PRODUCT_HPP
#define ALPS_NUMERIC_SCALAR_PRODUCT_HPP

#include <alps/numeric/matrix/scalar_product.hpp>
#include <alps/functional.h>
#include <alps/type_traits/element_type.hpp>

#include <algorithm>
#include <functional>
#include <numeric>
#include <valarray>

namespace alps { namespace numeric {

// The generic implementation of the scalar_product moved to alps/numeric/matrix/scalar_product.hpp

/// \overload
template <class T>
inline T scalar_product(const std::valarray<T>& c1, const std::valarray<T>& c2) 
{
  return std::inner_product(data(c1),data(c1)+c1.size(),data(c2),T(), std::plus<T>(),conj_mult<T,T>());
}

} } // namespace alps::numeric

#endif // ALPS_NUMERIC_SCALAR_PRODUCT_HPP
