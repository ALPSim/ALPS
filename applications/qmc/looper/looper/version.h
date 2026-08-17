/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
  alps::print_license(os);
  return os;
}

} // end namespace looper

#endif // LOOPER_VERSION_H
