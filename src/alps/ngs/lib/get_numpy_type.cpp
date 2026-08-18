/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/detail/get_numpy_type.hpp>

namespace alps {
    namespace detail {

        int get_numpy_type(bool) { return NPY_BOOL; }
        int get_numpy_type(char) { return NPY_CHAR; }
        int get_numpy_type(unsigned char) { return NPY_UBYTE; }
        int get_numpy_type(signed char) { return NPY_BYTE; }
        int get_numpy_type(short) { return NPY_SHORT; }
        int get_numpy_type(unsigned short) { return NPY_USHORT; }
        int get_numpy_type(int) { return NPY_INT; }
        int get_numpy_type(unsigned int) { return NPY_UINT; }
        int get_numpy_type(long) { return NPY_LONG; }
        int get_numpy_type(unsigned long) { return NPY_ULONG; }
        int get_numpy_type(long long) { return NPY_LONGLONG; }
        int get_numpy_type(unsigned long long) { return NPY_ULONGLONG; }
        int get_numpy_type(float) { return NPY_FLOAT; }
        int get_numpy_type(double) { return NPY_DOUBLE; }
        int get_numpy_type(long double) { return NPY_LONGDOUBLE; }
        int get_numpy_type(std::complex<float>) { return NPY_CFLOAT; }
        int get_numpy_type(std::complex<double>) { return NPY_CDOUBLE; }
        int get_numpy_type(std::complex<long double>) { return NPY_CLONGDOUBLE; }
    }
}
