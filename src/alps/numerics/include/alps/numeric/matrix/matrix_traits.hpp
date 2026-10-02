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


#ifndef ALPS_MATRIX_TRAITS_HPP
#define ALPS_MATRIX_TRAITS_HPP

namespace alps {
namespace numeric {

    template <typename Matrix>
    struct associated_diagonal_matrix
    {
    };

    template <typename Matrix>
    struct associated_real_diagonal_matrix
    {
    };

    template <typename Matrix>
    struct associated_vector
    {
    };

    template <typename Matrix>
    struct associated_real_vector
    {
    };

} // end namespace numeric
} // end namespace alps
#endif //ALPS_MATRIX_TRAITS_HPP
