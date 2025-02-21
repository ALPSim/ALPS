/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003 by Synge Todo <wistaria@comp-phys.org>
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

/* $Id$ */

#include <alps/config.h>
#include <cstddef>
#include <iostream>

int main() {

#define DO_TYPE(T) \
  std::cout << "size of "#T" is " << sizeof(T) << std::endl;

  DO_TYPE(bool)
  DO_TYPE(char)
  DO_TYPE(short)
  DO_TYPE(int)
  DO_TYPE(long)
  DO_TYPE(long long)
  DO_TYPE(float)
  DO_TYPE(double)
  DO_TYPE(long double)

  DO_TYPE(alps::int8_t)
  DO_TYPE(alps::int16_t)
  DO_TYPE(alps::int32_t)
  DO_TYPE(alps::int64_t)

  DO_TYPE(std::size_t)
  DO_TYPE(std::ptrdiff_t)

  return 0;
}
