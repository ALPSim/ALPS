/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2006 by Matthias Troyer <troyer@comp-phys.org>,
*                            Andreas Honecker <ahoneck@uni-goettingen.de>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <cstddef>
#include "factory.h"
#include "fulldiag.h"

alps::scheduler::Task* FullDiagFactory::make_task(const alps::ProcessList& w, const boost::filesystem::path& fn, const alps::Parameters& parms) const
{
  alps::LatticeLibrary lib(parms);
  bool is_graph = false;
  if (parms.defined("GRAPH"))
    is_graph = lib.has_graph("GRAPH");
  if (parms.defined("LATTICE"))
    is_graph = is_graph || lib.has_graph("LATTICE");  
  alps::params const settings =
    FullDiagMatrix<double>::read_parameters(parms, &FullDiagMatrix<double>::define_parameters);
  return FullDiagMatrix<double>::use_complex(settings, is_graph)  ?
    static_cast<alps::scheduler::Task*>(new FullDiagMatrix<std::complex<double> >(w,fn,&settings)) :
    static_cast<alps::scheduler::Task*>(new FullDiagMatrix<double>(w,fn,&settings));
}
  
void FullDiagFactory::print_copyright(std::ostream& out) const
{
  out << "ALPS full diagonalization application version 1.1\n"
      << "  available from http://alps.comp-phys.org/\n"
      << "  copyright (c) 2003-2007 by Matthias Troyer <troyer@comp-phys.org>\n"
      << "                          and Andreas Honecker <ahoneck@uni-goettingen.de>\n"
      << "\n";
}
