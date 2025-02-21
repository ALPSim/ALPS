/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
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

#include <alps/parapack/clone_info.h>
#include <alps/osiris/comm.h>
#if defined(ALPS_HAVE_UNISTD_H)
# include <unistd.h>
#elif defined(ALPS_HAVE_WINDOWS_H)
# include <windows.h>
#endif

int main(int argc, char **argv) {
  alps::comm_init(argc, argv);
  alps::Parameters params;
  params["SEED"] = 29832;
  alps::clone_info info(0, params, "info_test");
  info.start("test 1");
  #if defined(ALPS_HAVE_UNISTD_H)
    sleep(1); // sleep 1 Sec
  #elif defined(ALPS_HAVE_WINDOWS_H)
    Sleep(1000); // sleep 1000 mSec
  #endif
  info.stop();
  #if defined(ALPS_HAVE_UNISTD_H)
    sleep(1); // sleep 1 Sec
  #elif defined(ALPS_HAVE_WINDOWS_H)
    Sleep(1000); // sleep 1000 mSec
  #endif
  info.start("test 2");
  #if defined(ALPS_HAVE_UNISTD_H)
    sleep(1); // sleep 1 Sec
  #elif defined(ALPS_HAVE_WINDOWS_H)
    Sleep(1000); // sleep 1000 mSec
  #endif
  info.stop();
  info.set_progress(0.593483);
  if (alps::is_master()) {
    alps::oxstream oxs;
    oxs << info;
  }
  info.set_progress(1);
  if (alps::is_master()) {
    alps::oxstream oxs;
    oxs << info;
  }
  alps::comm_exit();
  return 0;
}
