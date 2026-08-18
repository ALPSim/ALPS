/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "util.h"
#include <alps/random/pseudo_des.h>
#include <boost/classic_spirit.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/throw_exception.hpp>
#include <stdexcept>

namespace alps {

int hash(int n, int s) {
  static const unsigned int hash_seed = 3777549;
  return (alps::pseudo_des::hash(s, n ^ hash_seed) ^ alps::pseudo_des::hash(s, hash_seed)) &
    ((1<<30) | ((1<<30)-1));
}

std::string id2string(int id, std::string const& pad) {
  int i = id;
  std::string str;
  while (i >= 10) {
    str += pad;
    i /= 10;
  }
  str += boost::lexical_cast<std::string>(id);
  return str;
}

double parse_percentage(std::string const& str) {
  using namespace boost::spirit;
  double r;
  if (!parse(str.c_str(), real_p[assign_a(r)] >> '%' >> end_p, space_p).full)
    boost::throw_exception(std::runtime_error("error in parsing \"" + str + '\"'));
  return 0.01 * r;
}

} // end namespace alps
