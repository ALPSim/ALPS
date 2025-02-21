/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_TYPE_H
#define LOOPER_TYPE_H

#include <boost/mpl/bool.hpp>

namespace looper {

//
// QMC types
//

struct classical {};
struct path_integral {};
struct sse {};

//
// meta functions
//

template<typename QMC>
struct is_path_integral
{ typedef boost::mpl::false_ type; };

template<>
struct is_path_integral<path_integral>
{ typedef boost::mpl::true_ type; };

template<typename QMC>
struct is_sse
{ typedef boost::mpl::false_ type; };

template<>
struct is_sse<sse>
{ typedef boost::mpl::true_ type; };

} // end namepspace looper

#endif // LOOPER_TYPE_H
