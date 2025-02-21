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

#ifndef ALPS_NUMERIC_DOUBLE2INT_HPP
#define ALPS_NUMERIC_DOUBLE2INT_HPP

#include <boost/numeric/conversion/converter.hpp>

namespace alps { namespace numeric {

//
// double2int
//

/// \brief rounds a floating point value to the nearest integer
/// ex) double2int(3.6) -> 3
///     double2int(1.2) -> 1
///     duoble2int(-0.7) -> -1 (!= int(-0.7 + 0.5))
///
/// \return nearest integer of the input
inline int double2int(double in) {
  typedef boost::numeric::converter<int, double, boost::numeric::conversion_traits<int, double>,
    boost::numeric::def_overflow_handler, boost::numeric::RoundEven<double> > converter;
  return converter::convert(in);
}

} } // end namespace alps::numeric

#endif // ALPS_NUMERIC_DOUBLE2INT_HPP
