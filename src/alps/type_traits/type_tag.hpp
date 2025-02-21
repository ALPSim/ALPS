/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1999-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
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

/* $Id$ */

#ifndef ALPS_TYPE_TRAITS_TYPE_TAG_H
#define ALPS_TYPE_TRAITS_TYPE_TAG_H

#include <alps/config.h>
#include <boost/mpl/int.hpp>
#include <complex>
#include <string>

namespace alps {

template <class T>
struct type_tag /* : public boost::mpl::int_<-1> */ {};

#define DEFINE_TYPE_TAG(TYPE,TAG) \
template<> struct type_tag< TYPE > : public boost::mpl::int_<TAG> {};

DEFINE_TYPE_TAG(float,0)
DEFINE_TYPE_TAG(double,1)
DEFINE_TYPE_TAG(long double,2)
DEFINE_TYPE_TAG(std::complex<float>,3)
DEFINE_TYPE_TAG(std::complex<double>,4)
DEFINE_TYPE_TAG(std::complex<long double>,5)
DEFINE_TYPE_TAG(int16_t,6)
DEFINE_TYPE_TAG(int32_t,7)
#ifndef BOOST_NO_INT64_T
DEFINE_TYPE_TAG(int64_t,8)
#endif
DEFINE_TYPE_TAG(uint16_t,9)
DEFINE_TYPE_TAG(uint32_t,10)
#ifndef BOOST_NO_INT64_T
DEFINE_TYPE_TAG(uint64_t,11)
#endif
DEFINE_TYPE_TAG(int8_t,12)
DEFINE_TYPE_TAG(uint8_t,13)
DEFINE_TYPE_TAG(std::string,14)
DEFINE_TYPE_TAG(bool,15)
#undef DEFINE_TYPE_TAG

} // namespace alps

#endif // ALPS_TYPE_TRAITS_TYPE_TAG_H
