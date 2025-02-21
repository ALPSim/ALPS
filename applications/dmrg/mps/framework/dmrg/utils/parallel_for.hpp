/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Alexandr Kosenkov <alex.kosenkov@gmail.com>
 *                         by Michele Dolfi <dolfim@phys.ethz.ch>
 * 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef PARALLEL_FOR_HPP
#define PARALLEL_FOR_HPP

#ifdef MAQUIS_OPENMP
    typedef std::size_t locale;
    typedef std::size_t locale_shared;
    #define parallel_pragma(a) _Pragma( #a )
    #define parallel_for(constraint, ...) parallel_pragma(omp parallel for schedule(dynamic, 1)) for(__VA_ARGS__)
    #define semi_parallel_for(constraint, ...) for(__VA_ARGS__)
#else
    typedef std::size_t locale;
    typedef std::size_t locale_shared;
    #define parallel_for(constraint, ...) for(__VA_ARGS__)
    #define semi_parallel_for(constraint, ...) for(__VA_ARGS__)
#endif

#endif
