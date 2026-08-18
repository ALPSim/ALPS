/****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: nobinning.h 3520 2009-12-11 16:49:53Z gamperl $ */

#include <alps/numeric/special_functions.hpp>
#include <alps/numeric/accumulate_if.hpp>

#include <iostream>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>
#include <functional>

#include <boost/bind.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/lambda/lambda.hpp>


int main(int argc, char** argv)
{
  std::vector<double> A;
  std::transform(boost::counting_iterator<int>(0), boost::counting_iterator<int>(30), std::back_inserter(A), 0.1 * boost::lambda::_1);

  std::cout << "\nA: \t";
  std::copy(A.begin(),A.end(),std::ostream_iterator<double>(std::cout,"\t"));
  std::cout << "\n";

  double conditional_sum = alps::numeric::accumulate_if
                                            ( A.begin()
                                            , A.end()
                                            , double()
                                            , boost::bind
                                               ( std::less_equal<double>()
                                               , boost::lambda::_1
                                               , 2.
                                               )
                                            );
  std::cout << "\nSum of all elements in A if <= 2. :\t " << conditional_sum << "\n";


  double conditional_sum_sq = alps::numeric::accumulate_if
                                                ( A.begin()
                                                , A.end()
                                                , double()
                                                , boost::bind
                                                   ( std::plus<double>()
                                                   , boost::lambda::_1
                                                   , boost::bind
                                                       ( static_cast<double (*)(double)>(&alps::numeric::sq)
                                                       , boost::lambda::_2
                                                       )
                                                   )
                                                , boost::bind
                                                   ( std::greater<double>()
                                                   , boost::lambda::_1
                                                   , 1.
                                                   )
                                                );

  std::cout << "\nSum of the square of all elements in A if > 1. :\t" << conditional_sum_sq << "\n";

  return 0;
}

