/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2006 by Matthias Troyer <troyer@comp-phys.org>,
*                            Andreas Honecker <ahoneck@uni-goettingen.de>
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
  bool transl = parms.value_or_default("TRANSLATION_SYMMETRY",true) || parms.defined("TOTAL_MOMENTUM");
  bool use_complex = (is_graph ? false : transl);
  return parms.value_or_default("COMPLEX",use_complex)  ?
    static_cast<alps::scheduler::Task*>(new FullDiagMatrix<std::complex<double> >(w,fn)) :
    static_cast<alps::scheduler::Task*>(new FullDiagMatrix<double>(w,fn));
}
  
void FullDiagFactory::print_copyright(std::ostream& out) const
{
  out << "ALPS full diagonalization application version 1.1\n"
      << "  available from http://alps.comp-phys.org/\n"
      << "  copyright (c) 2003-2007 by Matthias Troyer <troyer@comp-phys.org>\n"
      << "                          and Andreas Honecker <ahoneck@uni-goettingen.de>\n"
      << " for details see the publication:\n"
      << "  A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).\n\n";
}
