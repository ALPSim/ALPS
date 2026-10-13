/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1999-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_NUMERIC_BINOMIAL_HPP
#define ALPS_NUMERIC_BINOMIAL_HPP

namespace alps { namespace numeric {

/// \brief calculate the binomial coefficient
/// \return the binomial coefficient l over n
inline std::size_t binomial(std::size_t l, std::size_t n)
{
  double nominator=1;
  double denominator=1;
  std::size_t n2=std::max BOOST_PREVENT_MACRO_SUBSTITUTION (n,l-n);
  std::size_t n1=std::min BOOST_PREVENT_MACRO_SUBSTITUTION (n,l-n);
  for (std::size_t i=n2+1;i<=l;i++)
    nominator*=i;
  for (std::size_t i=2;i<=n1;i++)
    denominator*=i;
  return std::size_t(nominator/denominator+0.1);
}
} } // end namespace alps::numeric

#endif // ALPS_NUMERIC_BINOMIAL_HPP
