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

#ifndef ALPS_SCHEDULER_FACTORY_H
#define ALPS_SCHEDULER_FACTORY_H

#include <alps/scheduler/task.h>
#include <alps/config.h>
#include <boost/filesystem/path.hpp>
#include <iostream>

namespace alps {
namespace scheduler {

//=======================================================================
// Factory
//
// a factory for user defined task and subtask objects
//-----------------------------------------------------------------------

class ALPS_DECL Factory
{
public:
  Factory() {}
  virtual ~Factory() {}
  virtual Task* make_task(const ProcessList&,const boost::filesystem::path&) const;
  virtual Task* make_task(const ProcessList&,const boost::filesystem::path&,const Parameters&) const;
  virtual Task* make_task(const ProcessList&,const Parameters&) const;
  virtual Worker* make_worker(const ProcessList&,const Parameters&,int) const;
  virtual void print_copyright(std::ostream&) const=0;
};

template <class TASK>
class SimpleFactory : public Factory
{
public:
  SimpleFactory() {}
  
  Task* make_task(const ProcessList& w,const boost::filesystem::path& fn) const
  {
    return new TASK(w,fn);
  }

  Task* make_task(const ProcessList& w,const Parameters& p) const
  {
    return new TASK(w,p);
  }
  
  void print_copyright(std::ostream& out) const
  {
    TASK::print_copyright(out);
  }
};

template <class TASK, class WORKER>
class BasicFactory : public SimpleFactory<TASK>
{
public:
  BasicFactory() {}
  
  Worker* make_worker(const ProcessList& where ,const Parameters& parms,int node) const
  {
    return new WORKER(where,parms,node);
  }

  void print_copyright(std::ostream& out) const
  {
    SimpleFactory<TASK>::print_copyright(out);
    WORKER::print_copyright(out);
  }
};

} // namespace scheduler
} // namespace alps
 
#endif
