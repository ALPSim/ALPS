/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * This software is part of the ALPS libraries, published under the ALPS           *
 * Library License; you can use, redistribute it and/or modify it under            *
 * the terms of the license, either version 1 or (at your option) any later        *
 * version.                                                                        *
 *                                                                                 *
 * You should have received a copy of the ALPS Library License along with          *
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also       *
 * available from http://alps.comp-phys.org/.                                      *
 *                                                                                 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        *
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT       *
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE       *
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,     *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER     *
 * DEALINGS IN THE SOFTWARE.                                                       *
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
