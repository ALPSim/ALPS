/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
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
