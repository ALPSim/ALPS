/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/scheduler/scheduler.h>
#include <boost/assert.hpp>
#include <alps/utility/copyright.hpp>
#include <alps/osiris/comm.h>

namespace alps {
namespace scheduler {

//=======================================================================
// SingleScheduler
// 
// a scheduler for a single task
//-----------------------------------------------------------------------

using namespace boost::posix_time;

SingleScheduler::SingleScheduler(const NoJobfileOptions& opt,const Factory& f)
  : Scheduler(opt,f)
{
  end_time=second_clock::local_time()+seconds(long(time_limit));
}

void SingleScheduler::create_task(Parameters const& p)
{
  destroy_task();
  theTask = proc.make_task(processes,p);
  if(theTask->cpus()>processes.size()) 
    boost::throw_exception(std::runtime_error("Task needs more CPUs than available"));
}

void SingleScheduler::destroy_task()
{
  if (theTask)
    delete theTask;
  theTask=0;
}

SingleScheduler::~SingleScheduler()
{
  destroy_task();
}

int SingleScheduler::run()
{
  BOOST_ASSERT(theTask != 0);
  ptime task_time(second_clock::local_time());
  std::cout << "Starting task.\n";
  
  theTask->start();

  ptime next_check=second_clock::local_time();
  ptime last_checkpoint=second_clock::local_time();
       
  bool task_finished = false;
  do {
    if (check_signals() == SignalHandler::TERMINATE)
      break;
         
    theTask->run();
          
    if(time_limit >0. && second_clock::local_time()>end_time) {
      std::cout << "Time limit exceeded\n";
      break;
    }
          
    if(second_clock::local_time()>next_check) {
      std::cout  << "Checking if it is finished: " << std::flush;
      double more_time=0;
      double percentage=0.;
      task_finished=theTask->finished(more_time,percentage);
              
      // next check after at more_time, restrained to min. and max. times
      more_time= (more_time < min_check_time ? min_check_time :
                  (more_time > max_check_time ? max_check_time : more_time));
      next_check=second_clock::local_time()+seconds(int(more_time));
      if(!task_finished)
        std::cout  << "not yet, next check in " << int(more_time) << " seconds ( "
        << static_cast<int>(100.*percentage) << "% done).\n";
    }
    if((!task_finished)&&(second_clock::local_time()>last_checkpoint+seconds(int(checkpoint_time)))) {
      // make regular checkpoints if not yet finished
      checkpoint();
      last_checkpoint=second_clock::local_time();
    }
  } while (!task_finished);
            
  theTask->halt();
  std::cout  << "This task took " << (second_clock::local_time()-task_time).total_seconds() << " seconds.\n";
  return task_finished ? 0 : -1;
}

// initialize a scheduler for real work, parsing the command line
SingleScheduler* start_single(const Factory& p, int argc, char** argv)
{
  alps::comm_init(argc,argv,false);
  if (is_master()) {
    p.print_copyright(std::cout);
    alps::scheduler::print_copyright(std::cout);
    alps::print_copyright(std::cout);
  }
  
  NoJobfileOptions opt;
  if (argc)
    opt = NoJobfileOptions(argc,argv);
  
  if (!opt.valid)
    return 0; // do not actually run
  
  if (is_master()) {
    SingleScheduler* s=new SingleScheduler(opt,p);
    theScheduler = s;
    return s;
  }
  else {
    theScheduler = new Scheduler(opt,p);
    theScheduler->run();
    return 0;
  }
}


void stop_single( bool exit_)
{  
  if (theScheduler)
    delete theScheduler; 
  if (exit_) alps::comm_exit();
}


} // namespace scheduler
} // namespace alps
