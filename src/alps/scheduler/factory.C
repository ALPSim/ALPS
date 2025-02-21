/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@comp-phys.org>
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
