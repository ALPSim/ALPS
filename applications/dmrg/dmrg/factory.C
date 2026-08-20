/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@comp-phys.org>,
*                            Adrian Feiguin <afeiguin@uwyo.edu>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <cstddef>
#include "factory.h"
#include "dmrg.h"

alps::scheduler::Task* DMRGFactory::make_task(const alps::ProcessList& w, const boost::filesystem::path& fn, const alps::Parameters& parms) const
{
  return parms.value_or_default("COMPLEX",false)  ?
    static_cast<alps::scheduler::Task*>(new DMRGTask<std::complex<double> >(w,fn)) :
    static_cast<alps::scheduler::Task*>(new DMRGTask<double>(w,fn));
}
  
void DMRGFactory::print_copyright(std::ostream& out) const
{
   print_dmrg_copyright(out);
}
