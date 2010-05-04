/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@comp-phys.org>,
*                            Adrian Feiguin <afeiguin@uwyo.edu>
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
