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
#ifndef ALPS_NUMERIC_MATRIX_DETAIL_DEBUG_OUTPUT_HPP
#define ALPS_NUMERIC_MATRIX_DETAIL_DEBUG_OUTPUT_HPP

#ifdef ALPS_NUMERIC_MATRIX_DEBUG
#include <typeinfo>
#include <iostream>
#define ALPS_NUMERIC_MATRIX_DEBUG_OUTPUT(T) \
    std::cerr << __FILE__ <<" " << __LINE__ << ":" << T << std::endl;
#else
#define ALPS_NUMERIC_MATRIX_DEBUG_OUTPUT(T)
#endif

#endif //ALPS_NUMERIC_MATRIX_DETAIL_DEBUG_OUTPUT_HPP
