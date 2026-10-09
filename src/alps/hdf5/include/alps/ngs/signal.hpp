/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_SIGNAL_HPP
#define ALPS_NGS_SIGNAL_HPP

#include <alps/ngs/config.hpp>

#include <boost/array.hpp>

#include <vector>

namespace alps {
  namespace ngs {

    class ALPS_DECL signal{

    public:

      /*!
Listen to the following posix signals SIGINT, SIGTERM, SIGXCPU, SIGQUIT, SIGUSR1, SIGUSR2, SIGSTOP SIGKILL. Those signals can be check by empty, top, pop

\verbatim embed:rst
.. note::
   If a SIGSEGV (segfault) or SIGBUS (bus error) occures, a stacktrace
   is printed an all open hdf5 archives are closed before it exits.
\endverbatim
      */
      signal();

      /*!
Returns if a signal has been captured.
      */
      bool empty();

      /*!
Returns the last signal that has been captured .
      */
      int top();


      /*!
Pops a signal form the stack.
       */
      void pop();

      /*!
Listen to the signals SIGSEGV (segfault) and SIGBUS (bus error). If one
of these signals are captured, a stacktrace is printed an all open hdf5
archives are closed before it exits. 
      */
      static void listen();

      static void slot(int signal);

      static void segfault(int signal);

    private:

      static std::size_t begin_;
      static std::size_t end_;
      static boost::array<int, 0x20> signals_;
    };
  }
}

#endif
