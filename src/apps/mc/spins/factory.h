/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003 by Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
