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

#include <alps/random/parallel/keyword.hpp>

#ifndef ALPS_RANDOM_SPRNG_KEYWORD_HPP
#define ALPS_RANDOM_SPRNG_KEYWORD_HPP

namespace alps { namespace random { namespace parallel {
  using boost::mpl::placeholders::_;
  
  BOOST_PARAMETER_KEYWORD(random_tag,parameter)

  /// INTERNAL ONLY
  typedef boost::parameter::parameters<
      boost::parameter::optional<random_tag::stream_number, boost::is_convertible<_,unsigned int> >
    , boost::parameter::optional<random_tag::total_streams, boost::is_convertible<_,unsigned int> >
    , boost::parameter::optional<random_tag::global_seed, boost::is_convertible<_,int> >
    , boost::parameter::optional<random_tag::parameter, boost::is_convertible<_,unsigned int> >
  > sprng_seed_params;

} } } // namespace alps::random::parallel

#endif // ALPS_RANDOM_SPRNG_KEYWORD_HPP
