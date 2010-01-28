/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2010 by Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
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
