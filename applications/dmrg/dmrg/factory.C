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
   out << "ALPS/dmrg version " DMRG_VERSION " (" DMRG_DATE ")\n"
       << "  Density Matrix Renormalization Group algorithm\n"
       << "  for low-dimensional interacting systems.\n"
       << "  available from http://alps.comp-phys.org/\n"
       << "  copyright (c) 2006-2013 by Adrian E. Feiguin\n"
       << "  for details see the publication: \n"
       << "  A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).\n\n";
}
