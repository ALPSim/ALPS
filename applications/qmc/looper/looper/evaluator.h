/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2010 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_EVALUATOR_H
#define LOOPER_EVALUATOR_H

#include <alps/scheduler.h>
#ifdef HAVE_PARAPACK_13
# include <alps/parapack/serial.h>
#else
# include <alps/parapack/worker.h>
#endif

namespace looper {

class abstract_evaluator : public alps::parapack::simple_evaluator {
public:
  virtual ~abstract_evaluator() {}
  virtual void pre_evaluate(alps::ObservableSet& m, alps::Parameters const&,
    alps::ObservableSet const& m_in) const = 0;
  virtual void evaluate(alps::scheduler::MCSimulation&, alps::Parameters const&,
    boost::filesystem::path const&) const = 0;
  virtual void evaluate(alps::ObservableSet& m, alps::Parameters const&,
    alps::ObservableSet const& m_in) const = 0;
  virtual void evaluate(alps::ObservableSet&) const {};
};

} // end namespace looper

#endif // LOOPER_EVALUATOR_H
