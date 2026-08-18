/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2005 by Matthias Troyer <troyer@comp-phys.org>,
*                       Andreas Streich <astreich@student.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
