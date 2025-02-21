/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2010 by Synge Todo <wistaria@comp-phys.org>
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

#include <alps/utility/copyright.hpp>
#include <alps/version.h>
#include <alps/utility/os.hpp>
#include <boost/version.hpp>
#include <iostream>

int main() {
  std::cout << "ALPS version:          " << alps::version() << std::endl
            << "Boost version:         " << BOOST_LIB_VERSION << std::endl
            << "source directory:      " << ALPS_SRCDIR << std::endl
            << "installed at:          " << ALPS_PREFIX << std::endl
            << "configured on:         " << alps::config_host() << std::endl
            << "configured by:         " << alps::config_user() << std::endl
            << "compiled on:           " << alps::compile_date() << std::endl
            << "current hostname:      " << alps::hostname() << std::endl
            << "current user:          " << alps::username() << std::endl;
}
