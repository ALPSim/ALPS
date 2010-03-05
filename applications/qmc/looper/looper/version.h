/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
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
