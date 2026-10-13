/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_DETAIL_REMOVE_CVR_HPP
#define ALPS_NGS_DETAIL_REMOVE_CVR_HPP

namespace alps {
    namespace detail {

        template<typename T> struct remove_cvr {
            typedef T type;
        };
    
        template<typename T> struct remove_cvr<T const> {
            typedef typename remove_cvr<T>::type type;
        };
    
        template<typename T> struct remove_cvr<T volatile> {
            typedef typename remove_cvr<T>::type type;
        };
    
        template<typename T> struct remove_cvr<T &> {
            typedef typename remove_cvr<T>::type type;
        };

    }
}
#endif
