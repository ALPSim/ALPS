/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
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

#include <alps/parapack/footprint.h>
#include <iostream>

class myclass
{
public:
  std::size_t footprint() const { return sizeof(*this); }
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
  double* ptr;
  std::cerr << "footprint of double* is " << alps::footprint(ptr) << std::endl;
  return 0;
}
