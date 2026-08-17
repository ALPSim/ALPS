/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Andreas Hehn <hehn@phys.ethz.ch>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_VECTOR_INTERFACE_HPP
#define ALPS_VECTOR_INTERFACE_HPP

#include <vector>

namespace alps {
namespace numeric {
    // This macro creates free functions that call member functions with the same
    // name, e.g. swap_cols(A,i,j) -> A.swap_cols(i,j)
#define IMPLEMENT_FORWARDING(RET,NAME,ARGS,VARS) \
template <typename Vector> \
RET NAME ARGS \
{ \
return m.NAME VARS; \
} 
    IMPLEMENT_FORWARDING(void, resize, (Vector& m, typename Vector::size_type i1), (i1) )
    
#undef IMPLEMENT_FORWARDING
} // end namespace numeric
} // end namespace alps


#endif //ALPS_VECTOR_INTERFACE_HPP

