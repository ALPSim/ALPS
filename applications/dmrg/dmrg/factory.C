/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@comp-phys.org>,
*                            Adrian Feiguin <afeiguin@uwyo.edu>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
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
