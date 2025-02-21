/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2005 by Matthias Troyer <troyer@comp-phys.org>,
*                       Andreas Streich <astreich@student.ethz.ch>
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

#include "fitting_scheduler.h"

#include <alps/scheduler/types.h>
#include <alps/osiris/comm.h>

namespace alps {
namespace scheduler {

int start_fitting(int argc, char** argv, const Factory& p) 
{
  if (argc < 2)
    return -1;

  comm_init(argc,argv);

  if ((is_master()) || (!runs_parallel())) {
    p.print_copyright(std::cout);
    alps::scheduler::print_copyright(std::cout);
    alps::print_copyright(std::cout);
  }

  NoJobfileOptions njfo(argc-1,argv);
  if (!(njfo.valid)) {
    std::cerr << "invalid options, returning \n";
    return -1;
  }

  if (!runs_parallel()) 
    theScheduler = new FittingScheduler<SingleScheduler>(njfo,p,argv[argc-1]);
  else if (is_master())
    theScheduler = new FittingScheduler<MPPScheduler>(njfo,p,argv[argc-1]);
  else
    theScheduler = new Scheduler(njfo,p);

  theScheduler->run();

  delete theScheduler;
  comm_exit();
  return 0;
}

} // namespace scheduler
} // namespace alps
