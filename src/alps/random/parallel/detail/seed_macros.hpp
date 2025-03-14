/* boost alps/parallel/lcg64.hpp header file
 *
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

#ifndef ALPS_RANDOM_PARALLEL_DETAIL_SEED_MACROS_HPP
#define ALPS_RANDOM_PARALLEL_DETAIL_SEED_MACROS_HPP

#include <boost/parameter/macros.hpp>
#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/repeat_from_to.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/array/elem.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>

#include <alps/random/parallel/keyword.hpp>

#define ALPS_RANDOM_PARALLEL_CONSTRUCTOR(z, n, P)                           \
  BOOST_PP_IF(n,template <, BOOST_PP_EMPTY())                               \
  BOOST_PP_ENUM_PARAMS(n,class T) BOOST_PP_IF(n,>, BOOST_PP_EMPTY() )       \
  BOOST_PP_ARRAY_ELEM(0,P) (BOOST_PP_ENUM_BINARY_PARAMS(n,T,const& x))      \
  BOOST_PP_ARRAY_ELEM(1,P)                                                  \
  {                                                                         \
    seed(BOOST_PP_ENUM_PARAMS(n,x) );                                       \
  }                                                                         \
                                                                            \
  template <class It BOOST_PP_COMMA_IF(n) BOOST_PP_ENUM_PARAMS(n,class T)>  \
  BOOST_PP_ARRAY_ELEM(0,P) (It& f, It const& l BOOST_PP_COMMA_IF(n)         \
          BOOST_PP_ENUM_BINARY_PARAMS(n,T,const& x))                        \
  BOOST_PP_ARRAY_ELEM(1,P)                                                  \
  {                                                                         \
    seed(f,l BOOST_PP_COMMA_IF(n) BOOST_PP_ENUM_PARAMS(n,x) );              \
  }                                                                         \
                                                                            \
  BOOST_PP_IF(n,template <, BOOST_PP_EMPTY())                               \
  BOOST_PP_ENUM_PARAMS(n,class T) BOOST_PP_IF(n,>, BOOST_PP_EMPTY())        \
  void seed(BOOST_PP_ENUM_BINARY_PARAMS(n,T,const& x))                      \
  {                                                                         \
    seed_named(BOOST_PP_ENUM_PARAMS(n,x) );                                 \
  }

#define ALPS_RANDOM_PARALLEL_ITERATOR_SEED_DEFAULT_IMPL(z, n, unused)       \
  template <class It BOOST_PP_COMMA_IF(n) BOOST_PP_ENUM_PARAMS(n,class T)>  \
  void seed(It& f, It const& l BOOST_PP_COMMA_IF(n)                         \
               BOOST_PP_ENUM_BINARY_PARAMS(n,T,const& x))                   \
  {                                                                         \
    if(f == l)                                                              \
      throw std::invalid_argument("invalid seeding argument");              \
    seed(alps::random::parallel::global_seed=*f++                           \
         BOOST_PP_COMMA_IF(n) BOOST_PP_ENUM_PARAMS(n,x) );                  \
  }
   
#define ALPS_RANDOM_PARALLEL_ITERATOR_SEED_IMPL(z, n, unused)               \
  template <class It BOOST_PP_COMMA_IF(n) BOOST_PP_ENUM_PARAMS(n,class T)>  \
  void seed(It& f, It const& l BOOST_PP_COMMA_IF(n)                         \
               BOOST_PP_ENUM_BINARY_PARAMS(n,T,const& x))    {              \
     iterator_seed_named(alps::random::parallel::first=f                    \
                        , alps::random::parallel::last=l                    \
                        BOOST_PP_COMMA_IF(n) BOOST_PP_ENUM_PARAMS(n,x) );   \
  }                                                                         \
       

#define ALPS_RANDOM_PARALLEL_SEED_PARAMS(RNG,PARAMS,INIT)                                           \
  BOOST_PP_REPEAT_FROM_TO(0, ALPS_RANDOM_MAXARITY, ALPS_RANDOM_PARALLEL_CONSTRUCTOR,(2,(RNG,INIT))) \
  BOOST_PARAMETER_MEMFUN(void,seed_named,0,ALPS_RANDOM_MAXARITY,PARAMS)     

#define ALPS_RANDOM_PARALLEL_SEED(RNG) \
  ALPS_RANDOM_PARALLEL_SEED_PARAMS(RNG, alps::random::parallel::seed_params,BOOST_PP_EMPTY()) 

#define ALPS_RANDOM_PARALLEL_ITERATOR_SEED_PARAMS(PARAMS)                                      \
  BOOST_PP_REPEAT_FROM_TO(0, BOOST_PP_SUB(ALPS_RANDOM_MAXARITY,2), ALPS_RANDOM_PARALLEL_ITERATOR_SEED_IMPL,~) \
  BOOST_PARAMETER_MEMFUN(void,iterator_seed_named,2,ALPS_RANDOM_MAXARITY,PARAMS)     

#define ALPS_RANDOM_PARALLEL_ITERATOR_SEED() \
  ALPS_RANDOM_PARALLEL_ITERATOR_SEED_PARAMS(alps::random::parallel::iterator_seed_params) 

#define ALPS_RANDOM_PARALLEL_ITERATOR_SEED_DEFAULT() \
  BOOST_PP_REPEAT_FROM_TO(0, BOOST_PP_SUB(ALPS_RANDOM_MAXARITY,1), ALPS_RANDOM_PARALLEL_ITERATOR_SEED_DEFAULT_IMPL,~)

#endif // ALPS_RANDOM_PARALLEL_DETAIL_SEED_MACROS_HPP
