/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2013 by Andreas Hehn <hehn@phys.ethz.ch>                          *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef ALPS_NUMERIC_MATRIX_ENTITY_HPP
#define ALPS_NUMERIC_MATRIX_ENTITY_HPP

namespace alps {
namespace numeric {

namespace tag {
    struct scalar {};
    struct vector {};
    struct matrix {};
}


template <typename T>
struct entity
{
    typedef tag::scalar type;
};

template <typename T>
struct get_entity : entity<typename boost::remove_const<T>::type>
{
};

} // end namespace numeric
} // end namespace alps

#endif // ALPS_NUMERIC_MATRIX_ENTITY_HPP
