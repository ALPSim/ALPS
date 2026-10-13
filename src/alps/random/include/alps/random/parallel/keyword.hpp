/* 
 * Copyright Matthias Troyer 2006
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
 */

/// \file keyword.hpp
/// \brief This file defines the named parameter keywords for the random library

#include <boost/parameter/keyword.hpp>
#include <boost/parameter/parameters.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/mpl/placeholders.hpp>

#ifndef ALPS_RANDOM_PARALLEL_KEYWORD_HPP
#define ALPS_RANDOM_PARALLEL_KEYWORD_HPP

#ifndef ALPS_RANDOM_MAXARITY
/// The number of named parameter arguments for the random library. The default value can be changed by defining the macro.
#define ALPS_RANDOM_MAXARITY 5
#endif


namespace alps { namespace random { namespace parallel {
  using boost::mpl::placeholders::_;
  
  BOOST_PARAMETER_KEYWORD(random_tag,global_seed)
  BOOST_PARAMETER_KEYWORD(random_tag,stream_number)
  BOOST_PARAMETER_KEYWORD(random_tag,total_streams)
  BOOST_PARAMETER_KEYWORD(random_tag,first)
  BOOST_PARAMETER_KEYWORD(random_tag,last)

  /// INTERNAL ONLY
  typedef boost::parameter::parameters<
      boost::parameter::optional<random_tag::global_seed>
    , boost::parameter::optional<random_tag::stream_number, boost::is_convertible<_,unsigned int> >
    , boost::parameter::optional<random_tag::total_streams, boost::is_convertible<_,unsigned int> >
  > seed_params;

  /// INTERNAL ONLY
  typedef boost::parameter::parameters<
      boost::parameter::required<random_tag::first>
    , boost::parameter::required<random_tag::last>
    , boost::parameter::optional<random_tag::stream_number, boost::is_convertible<_,unsigned int> >
    , boost::parameter::optional<random_tag::total_streams, boost::is_convertible<_,unsigned int> >
  > iterator_seed_params;

} } } // namespace alps::random::parallel

#endif // ALPS_RANDOM_PARALLEL_KEYWORD_HPP
