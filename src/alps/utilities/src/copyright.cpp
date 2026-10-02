/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2011 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/config.h>
#include <alps/utility/copyright.hpp>
#include <alps/version.h>

void alps::print_copyright(std::ostream& out) {
  out << "based on the ALPS libraries version " << ALPS_VERSION << "\n";
  out << "  available from https://alps.comp-phys.org/\n";
  out << "  copyright (c) 1994-" << ALPS_YEAR
      << " by the ALPS collaboration.\n";
  out << "  Licensed under the MIT License.\n";
  out << "  License text: https://github.com/ALPSim/ALPS/blob/master/LICENSE.txt\n";
  out << "  For details see the publication: \n"
      << "  B. Bauer et al., J. Stat. Mech. (2011) P05001.\n\n";
}

void alps::print_license(std::ostream& out) {
  out << "Licensed under the MIT License.\n";
  out << "License text: https://github.com/ALPSim/ALPS/blob/master/LICENSE.txt\n";
}

std::string alps::version() { return ALPS_VERSION; }

std::string alps::version_string() { return ALPS_VERSION_STRING; }

std::string alps::year() { return ALPS_YEAR; }

std::string alps::config_host() { return ALPS_CONFIG_HOST; }

std::string alps::config_user() { return ALPS_CONFIG_USER; }

std::string alps::compile_date() {
#if defined(__DATE__) && defined(__TIME__)
  return __DATE__ " " __TIME__;
#else
  return "unknown";
#endif
}
