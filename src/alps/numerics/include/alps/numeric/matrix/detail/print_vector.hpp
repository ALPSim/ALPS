/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Andreas Hehn <hehn@phys.ethz.ch>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef ALPS_NUMERIC_MATRIX_DETAIL_PRINT_VECTOR_HPP
#define ALPS_NUMERIC_MATRIX_DETAIL_PRINT_VECTOR_HPP

#include <ostream>

namespace alps {
namespace numeric {
namespace detail {

template <typename Vector>
void print_vector(std::ostream& os, Vector const& v)
{
    os<<"[";
    if(v.size() > 0)
    {
        for(unsigned int i=0;i<v.size()-1;++i)
          os<<v(i)<<", ";
        os<< v(v.size()-1);
    }
    os << "]"<<std::endl;
}

} // end namespace detail
} // end namespace numeric
} // end namespace alps

#endif // ALPS_NUMERIC_MATRIX_DETAIL_PRINT_VECTOR_HPP
