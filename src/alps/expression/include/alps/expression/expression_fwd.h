/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2006 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_EXPRESSION_EXPRESSION_FWD_H
#define ALPS_EXPRESSION_EXPRESSION_FWD_H

#include <alps/config.h>

#include <alps/cctype.h>
#include <alps/parameter.h>
#include <alps/random.h>
#include <alps/parser/parser.h>
#include <alps/utility/vectorio.hpp>
#include <alps/type_traits/is_symbolic.hpp>

#include <boost/call_traits.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/throw_exception.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/bool.hpp>

#include <cmath>
#include <complex>
#include <string>
#include <vector>
#include <stdexcept>

namespace alps {
namespace expression {

template<class T = std::complex<double> > class Expression;
template<class T = std::complex<double> > class Term;
template<class T = std::complex<double> > class Factor;
template<class T = std::complex<double> > class Evaluator;
template<class T = std::complex<double> > class ParameterEvaluator;
template <class T> class Block;
template <class T> class Function;
template <class T> class Number;
template <class T> class Symbol;

}

typedef expression::Expression<> Expression;
typedef expression::Term<> Term;
typedef expression::Factor<> Factor;
typedef expression::Evaluator<> Evaluator;
typedef expression::ParameterEvaluator<> ParameterEvaluator;

template <class T>
struct is_symbolic<expression::Expression<T> > : public boost::mpl::true_ {};

template <class T>
struct is_symbolic<expression::Term<T> > : public boost::mpl::true_ {};

template <class T>
struct is_symbolic<expression::Factor<T> > : public boost::mpl::true_ {};

template <class T>
struct is_symbolic<expression::Block<T> > : public boost::mpl::true_ {};

template <class T>
struct is_symbolic<expression::Function<T> > : public boost::mpl::true_ {};

template <class T>
struct is_symbolic<expression::Symbol<T> > : public boost::mpl::true_ {};



}

#endif
