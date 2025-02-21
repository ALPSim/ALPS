/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2009 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

#ifndef ALPS_SIGNAL_H
#define ALPS_SIGNAL_H

//=======================================================================
// This file defines a signal handler class
//=======================================================================

#include <alps/config.h>

namespace alps {
namespace scheduler {

/** implements a signal handler.
    signals are intercepted and can be checked for.
*/

class ALPS_DECL SignalHandler
{
public:
  /** symbolic names for signals.
    SIGINT, SIGQUIT and SIGTERM are mapped to TERMINATE
    SIGTSTP is mapped to STOP
    SIGUSR1 is mapped to USER1
    SIGUSR2 is mapped to USER2*/

    enum SignalInfo
    {
      NOSIGNAL=0,
      USER1,
      USER2,
      STOP,
      TERMINATE
    };

  /// a default constructor
  SignalHandler();

  /** ask for signals.
      If more than one signal has been received the signal
      with the highest priority will be returned. Priorities are:
      USER1 > USER2 > STOP > TERMINATE. */

  SignalInfo operator()();

  /// send myself a noncatchable stop signal
  static void stopprocess();

private:
  static unsigned int u1; // number of user1 signals received
  static unsigned int u2; // number of user2 signals received
  static unsigned int k; // number of terminate signals received
  static unsigned int s; // number of stop signals received
  static unsigned int count; // total number of signals received
  static bool initialized;

  // functions to be called by the signal handlers, to register them with
  // this object, only for internal use

  static void tstp(int);
  static void kill(int);
  static void usr1(int);
  static void usr2(int);
};

} // end namespace
} // end namespace

#endif // ALPS_SIGNAL_H
