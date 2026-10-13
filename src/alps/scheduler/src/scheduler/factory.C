/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/scheduler/factory.h>

namespace alps {
namespace scheduler {

Worker* Factory::make_worker(const ProcessList&,const Parameters&,int) const
{
  boost::throw_exception(std::logic_error("Factory::make_worker() needs to be implemented"));
  return 0;
}

Task* Factory::make_task(const ProcessList& w,const boost::filesystem::path& fn) const
{
  alps::Parameters parms;
  { // scope to close file
    boost::filesystem::ifstream infile(fn);
    parms.extract_from_xml(infile);
  }
  return make_task(w,fn,parms);
}

Task* Factory::make_task(const ProcessList&,const boost::filesystem::path&,const Parameters&) const
{
  boost::throw_exception(std::logic_error("Factory::make_task(const ProcessList&,const boost::filesystem::path&,const Parameters&) needs to be implemented"));
  return 0;
}

Task* Factory::make_task(const ProcessList&,const Parameters&) const
{
  boost::throw_exception(std::logic_error("Factory::make_task(const ProcessList&,const Parameters&) needs to be implemented"));
  return 0;
}

} // namespace scheduler
} // namespace alps
