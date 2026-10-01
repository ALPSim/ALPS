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
  bool transl = parms.value_or_default("TRANSLATION_SYMMETRY",true) || parms.defined("TOTAL_MOMENTUM");
  bool use_complex = (is_graph ? false : transl);
  return parms.value_or_default("COMPLEX",use_complex) ?
    static_cast<alps::scheduler::Task*>(new SparseDiagMatrix<std::complex<double> >(w,fn)) :
    static_cast<alps::scheduler::Task*>(new SparseDiagMatrix<double>(w,fn));
}
  
void SparseDiagFactory::print_copyright(std::ostream& out) const
{
  out << "ALPS sparse diagonalization application version 1.1\n"
      << "  available from http://alps.comp-phys.org/\n"
      << "  copyright (c) 2003-2007 by Matthias Troyer <troyer@comp-phys.org>\n"
      << " for details see the publication:\n"
      << "  A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).\n\n";
}
