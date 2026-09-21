/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2007 - 2010 by Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: lattice.h 1744 2005-08-26 07:55:30Z wistaria $ */

#include <cmath>
double factorial(unsigned int x)
{
  double res=1.;
  while (x)
    res *=x--;
  return res;
}


double probability(int L, int N)
{
  double probability=factorial(L)/factorial(N)/factorial(L-N);
  return probability * std::pow(0.5,L);
}
