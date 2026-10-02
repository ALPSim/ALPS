/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/integer_range.h>
#include <cstdlib>
#include <iostream>
#include <string>

template<class T>
void check_range_boundaries() {
  const T lo = (std::numeric_limits<T>::min)();
  const T hi = (std::numeric_limits<T>::max)();
  const alps::integer_range<T> full(lo, hi), empty;
  if (full.empty() || !full.valid() || !empty.empty() || empty.size() != 0 ||
      alps::integer_range<T>(lo).size() != 1 || alps::integer_range<T>(hi).size() != 1 ||
      alps::integer_range<T>(0, hi - 1).size() != hi ||
      unify(full, full) != full)
    throw std::runtime_error("integer_range boundary check failed");
  try {
    full.size();
  } catch (std::overflow_error const&) {
    return;
  }
  throw std::runtime_error("integer_range failed to reject an unrepresentable size");
}

int main()
{
  check_range_boundaries<int>();
  check_range_boundaries<unsigned int>();
  check_range_boundaries<long long>();
  check_range_boundaries<unsigned long long>();
  std::string str;
  alps::Parameters params;
  while (std::getline(std::cin, str) && str.size())
    params << alps::Parameter(str);
  std::cout << "Parameters:\n" << params;

  std::cout << "Test for integer_range<int>:\n";
  while (std::getline(std::cin, str)) {
    if (str.size() == 0 || str[0] == '\0') break;
    std::cout << "parse " << str << ": ";
    try {
      alps::integer_range<int> r(str, params);
      std::cout << "result " << r << std::endl;
    }
    catch (std::exception& exp) {
      std::cout << exp.what() << std::endl;
    }
  }
  std::cout << "Test for integer_range<unsigned int>:\n";
  while (std::getline(std::cin, str)) {
    if (!std::cin || str[0] == '\0') break;
    std::cout << "parse " << str << ": ";
    try {
      alps::integer_range<unsigned int> r(str, params);
      std::cout << "result " << r << std::endl;
    }
    catch (std::exception& exp) {
      std::cout << exp.what() << std::endl;
    }
  }

  alps::integer_range<int> r(0, 5);
  std::cout << "initial: " << r << std::endl;
  std::cout << "7 is included? " << r.is_included(7) << std::endl;
  r = 3;
  std::cout << "3 is assigned: " << r << std::endl;
  r.include(8);
  std::cout << "8 is included: " << r << std::endl;
  std::cout << "7 is included? " << r.is_included(7) << std::endl;

  alps::integer_range<int> s("[3:]");
  std::cout << "initial: " << s << std::endl;
  std::cout << "multiplied by 3.5: " << 3.5 * s << std::endl;
  std::cout << "multiplied by 2000000000: " << s * 2000000000 << std::endl;
  s *= 0.1;
  std::cout << "multiplied by 0.1: " << s << std::endl;

  alps::integer_range<int> t(0, 5);
  alps::integer_range<int> u(-2, 3);
  alps::integer_range<int> v(-2, 10);
  alps::integer_range<int> w(7, 10);
  std::cout << "initial: " << t << std::endl;
  std::cout << "overlap with " << u << ": " << overlap(t, u) << std::endl;
  std::cout << "union with " << u << ": " << unify(t, u) << std::endl;
  std::cout << "overlap with " << v << ": " << overlap(t, v) << std::endl;
  std::cout << "union with " << v << ": " << unify(t, v) << std::endl;
  std::cout << "overlap with " << w << ": " << overlap(t, w) << std::endl;
  std::cout << "union with " << w << ": ";
  try { std::cout << unify(t, w) << std::endl; }
  catch (std::exception& exp) { std::cout << exp.what() << std::endl; }
  return 0;
}
