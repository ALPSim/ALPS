/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1999-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_MATH_HPP
#define ALPS_MATH_HPP

#if defined(_MSC_VER) || defined(__BORLANDC__) || defined(__DMC__)
#  pragma message ("This header is deprecated. Please use the new headers in alps/numeric")
#elif defined(__GNUC__) || defined(__HP_aCC) || defined(__SUNPRO_CC) || defined(__IBMCPP__)
#  warning "This header is deprecated. Please use the new headers in alps/numeric"
#endif

#include <alps/numeric/real.hpp>
#include <alps/numeric/abs2.hpp>
#include <alps/numeric/binomial.hpp>
#include <alps/numeric/is_equal.hpp>
#include <alps/numeric/is_zero.hpp>
#include <alps/numeric/is_nonzero.hpp>
#include <alps/numeric/round.hpp>
#include <alps/numeric/is_negative.hpp>
#include <alps/numeric/is_positive.hpp>
#include <alps/numeric/double2int.hpp>

namespace alps {

using numeric::real;
using numeric::binomial;
using numeric::abs2;
using numeric::is_equal;
using numeric::is_zero;
using numeric::is_nonzero;
using numeric::round;
using numeric::is_negative;
using numeric::is_positive;
using numeric::double2int;

} // end namespace

#endif // ALPS_MATH_HPP
