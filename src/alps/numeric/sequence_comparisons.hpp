/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* Copyright (C) 2011-2012 by Lukas Gamper <gamperl@gmail.com>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Maximilian Poprawe <poprawem@ethz.ch>
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

#ifndef ALPS_UTILITY_SEQUENCE_COMPARISONS_HPP
#define ALPS_UTILITY_SEQUENCE_COMPARISONS_HPP

#include <alps/type_traits/is_sequence.hpp>
#include <alps/type_traits/element_type.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/bind.hpp>

#include <algorithm>

namespace alps {
namespace numeric {

  template <class X, class BinaryPredicate>                                                                                                                     
  inline typename boost::disable_if< is_sequence<X>, bool>::type                                                                         
        at_least_one (const X& value1, const X& value2, const BinaryPredicate& pred) {                                                                               
    return pred(value1, value2);                                                                                                  
  }                                                                                                                                      
                                                                                                                                         
  template <class X, class BinaryPredicate>                                                                                                                   
  inline typename boost::enable_if< is_sequence<X>, bool>::type                                                                          
        at_least_one (const X& sequence, const typename element_type<X>::type & value, const BinaryPredicate& pred) {                                                
    return sequence.end() != std::find_if(sequence.begin(), sequence.end(), boost::bind<bool>(pred, boost::lambda::_1, value) );                      
  }                                                                                                                                      
                                                                                                                                         
  template <class X, class BinaryPredicate>                                                                                                                   
  inline typename boost::enable_if< is_sequence<X>, bool>::type                                                                          
        at_least_one (const typename element_type<X>::type & value, const X& sequence, const BinaryPredicate& pred) {                                                
    return sequence.end() != std::find_if(sequence.begin(), sequence.end(), boost::bind<bool>(pred, value, boost::lambda::_1 ) );                      
  }                                                                                                                                      
                                                                                                                                         
  template <class X, class BinaryPredicate>                                                                                                                    
  inline typename boost::enable_if< is_sequence<X>, bool>::type                                                                          
        at_least_one (const X& sequence1, const X& sequence2, const BinaryPredicate& pred) {                                                                         
    return !(std::equal(sequence1.begin(), sequence1.end(), sequence2.begin(), !boost::bind<bool>(pred, boost::lambda::_1, boost::lambda::_2)));     
  }


} // ending namespace numeric
} // ending namespace alps

#endif

