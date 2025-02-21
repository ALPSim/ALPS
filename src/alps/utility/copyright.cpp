/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2011 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
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
#include <alps/utility/copyright.hpp>
#include <alps/version.h>

void alps::print_copyright(std::ostream& out) {
  out << "based on the ALPS libraries version " << ALPS_VERSION << "\n";
  out << "  available from http://alps.comp-phys.org/\n";
  out << "  copyright (c) 1994-" << ALPS_YEAR
      << " by the ALPS collaboration.\n";
  out << "  Consult the web page for license details.\n";
  out << "  For details see the publication: \n"
      << "  B. Bauer et al., J. Stat. Mech. (2011) P05001.\n\n";
}

void alps::print_license(std::ostream& out) {
  out << "Please look at the file LICENSE.txt for the license conditions\n";
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
