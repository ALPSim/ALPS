/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003 by Matthias Troyer <troyer@comp-phys.org>
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

/* $Id$ */

#ifndef ALPS_APPLICATIONS_MC_SPIN_FACTORY_H_
#define ALPS_APPLICATIONS_MC_SPIN_FACTORY_H_

#include <alps/scheduler/montecarlo.h>
#include "abstractspinsim.h"

class SpinFactory : public alps::scheduler::Factory
{
public:
  SpinFactory() {}
  
  alps::scheduler::MCSimulation* make_task(const alps::ProcessList& w,
          const boost::filesystem::path& fn) const;

  alps::scheduler::Worker* make_worker(const alps::ProcessList& where,
          const alps::Parameters& parms, int node) const;
  void print_copyright(std::ostream&) const;

private:
  int countElements(const std::string& str) const;
  int findDominantMatrixString(const alps::Parameters& parms) const;
  void produceError(const alps::Parameters& parms) const;
};

#endif
