/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2006 by Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include "factory.h"
#include "sparsediag.h"

alps::scheduler::Task* SparseDiagFactory::make_task(const alps::ProcessList& w, const boost::filesystem::path& fn, const alps::Parameters& parms) const
{
  alps::LatticeLibrary lib(parms);
  bool is_graph = false;
  if (parms.defined("GRAPH"))
    is_graph = is_graph || lib.has_graph("GRAPH");
  if (parms.defined("LATTICE"))
      is_graph = is_graph || lib.has_graph("LATTICE");
  alps::params const settings =
    SparseDiagMatrix<double>::read_parameters(parms, &SparseDiagMatrix<double>::define_parameters);
  return SparseDiagMatrix<double>::use_complex(settings, is_graph) ?
    static_cast<alps::scheduler::Task*>(new SparseDiagMatrix<std::complex<double> >(w,fn)) :
    static_cast<alps::scheduler::Task*>(new SparseDiagMatrix<double>(w,fn));
}
  
void SparseDiagFactory::print_copyright(std::ostream& out) const
{
  out << "ALPS sparse diagonalization application version 1.1\n"
      << "  available from http://alps.comp-phys.org/\n"
      << "  copyright (c) 2003-2007 by Matthias Troyer <troyer@comp-phys.org>\n"
      << "\n";
}
