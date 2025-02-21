/*****************************************************************************
*
* ALPS Project Applications
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

#ifndef LOOPER_VERSION_H
#define LOOPER_VERSION_H

#include <iostream>

/* Define the version of ALPS/looper */
#define LOOPER_VERSION "3.2b12-20100128"

/* Define the published date of ALPS/looper */
#define LOOPER_DATE "2010/01/28"

#define LOOPER_VERSION_STRING "ALPS/looper version " LOOPER_VERSION " (" LOOPER_DATE ")"

#define LOOPER_COPYRIGHT LOOPER_VERSION_STRING "\n" \
  "  multi-cluster quantum Monte Carlo algorithms for spin systems\n" \
  "  available from http://wistaria.comp-phys.org/alps-looper/\n" \
  "  copyright (c) 1997-2010 by Synge Todo <wistaria@comp-phys.org>\n" \

#include <alps/utility/copyright.hpp>
#include <iostream>

namespace looper {

inline std::string version() {
  return LOOPER_VERSION_STRING;
}

inline std::ostream& print_copyright(std::ostream& os = std::cout) {
  os << LOOPER_COPYRIGHT << "\n";
  return os;
}

inline std::ostream& print_license(std::ostream& os = std::cout) {
  os << "Please look at the file LICENSE for the license conditions.\n";
  return os;
}

} // end namespace looper

#endif // LOOPER_VERSION_H
