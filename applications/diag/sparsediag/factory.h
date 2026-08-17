/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2004 by Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_APPLICATIONS_SPARSEDIAG_FACTORY_H
#define ALPS_APPLICATIONS_SPARSEDIAG_FACTORY_H

#include <alps/scheduler/factory.h>

class SparseDiagFactory : public alps::scheduler::Factory {
  alps::scheduler::Task* make_task(const alps::ProcessList&, 
     const boost::filesystem::path&, const alps::Parameters&) const;  
  void print_copyright(std::ostream& out) const;
};

#endif
