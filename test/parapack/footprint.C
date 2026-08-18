/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/footprint.h>
#include <iostream>

class myclass
{
public:
  std::size_t footprint() const { (void)a; (void)b; return sizeof(*this); }
private:
  int a;
  double b;
};

int main() {
  int a = 0;
  std::cerr << "footprint of int is " << alps::footprint(a) << std::endl;
  std::vector<int> b(30);
  std::cerr << "footprint of std::vector<int> is " << alps::footprint(b) << std::endl;
  std::string c("my string");
  std::cerr << "footprint of std::string is " << alps::footprint(c) << std::endl;
  myclass d;
  std::cerr << "footprint of myclass is " << alps::footprint(d) << std::endl;
  double* ptr = nullptr;
  std::cerr << "footprint of double* is " << alps::footprint(ptr) << std::endl;
  return 0;
}
