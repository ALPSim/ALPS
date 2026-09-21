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

#include <alps/parapack/util.h>
#include <iostream>

void test(std::string const& str) {
  std::cout << '\"' << str << "\" is parsed as " << alps::parse_percentage(str) << std::endl;
}

int main() {
  test("10%");
  test("10 %");
  test(" 10%");
  test(" 10 %  ");
  test("0.1%");
  return 0;
}
