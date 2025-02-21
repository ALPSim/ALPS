/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MAQUIS_DMRG_MODELS_CODED_FACTORY_LATTICE_H
#define MAQUIS_DMRG_MODELS_CODED_FACTORY_LATTICE_H

#include "dmrg/models/coded/lattice.hpp"

inline boost::shared_ptr<lattice_impl>
coded_lattice_factory(BaseParameters & parms)
{
    typedef boost::shared_ptr<lattice_impl> impl_ptr;
    if (parms["LATTICE"] == std::string("periodic chain lattice"))
        return impl_ptr(new ChainLattice(parms, true));
    else if (parms["LATTICE"] == std::string("chain lattice"))
        return impl_ptr(new ChainLattice(parms, false));
    else if (parms["LATTICE"] == std::string("open chain lattice"))
        return impl_ptr(new ChainLattice(parms, false));
    else if (parms["LATTICE"] == std::string("square lattice"))
        return impl_ptr(new SquareLattice(parms));
    else if (parms["LATTICE"] == std::string("open square lattice"))
        return impl_ptr(new SquareLattice(parms));
    else {
        throw std::runtime_error("Don't know this lattice!");
    }
}

#endif
