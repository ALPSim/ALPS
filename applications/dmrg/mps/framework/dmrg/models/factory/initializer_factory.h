/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MAQUIS_DMRG_MODELS_INITIALIZER_FACTORY_H
#define MAQUIS_DMRG_MODELS_INITIALIZER_FACTORY_H

#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps_initializers.h"


template <class Matrix, class SymmGroup>
typename model_impl<Matrix,SymmGroup>::initializer_ptr
model_impl<Matrix,SymmGroup>::initializer(Lattice const& lat, BaseParameters & parms) const
{
    typename SymmGroup::charge initc = this->total_quantum_numbers(parms);
    
    int max_site_type = 0;
    std::vector<int> site_types(lat.size(), 0);
    for (int p = 0; p < lat.size(); ++p) {
        site_types[p] = lat.get_prop<int>("type", p);
        max_site_type = std::max(site_types[p], max_site_type);
    }
    
    std::cout << "site_types: ";
    std::copy(site_types.begin(), site_types.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;
    
    std::vector<Index<SymmGroup> > site_bases(max_site_type+1);
    for (int type = 0; type < site_bases.size(); ++type) {
        site_bases[type] = this->phys_dim(type);
        std::cout << "phys["<<type <<"]: " << site_bases[type] << std::endl;
    }
    
    if (parms["init_state"] == "default")
        return initializer_ptr(new default_mps_init<Matrix, SymmGroup>(parms, site_bases, initc, site_types));
    
    else if (parms["init_state"] == "const")
        return initializer_ptr(new const_mps_init<Matrix, SymmGroup>(parms, site_bases, initc, site_types));
    
    else if (parms["init_state"] == "orthogonal")
        return initializer_ptr(new qr_mps_init<Matrix, SymmGroup>(parms, site_bases, initc, site_types));

    else if (parms["init_state"] == "ortho_symm")
        return initializer_ptr(new qr_symm_mps_init<Matrix, SymmGroup>(parms, site_bases, initc, site_types));
    
    else if (parms["init_state"] == "thin")
        return initializer_ptr(new thin_mps_init<Matrix, SymmGroup>(parms, site_bases, initc, site_types));
    
    else if (parms["init_state"] == "thin_const")
        return initializer_ptr(new thin_const_mps_init<Matrix, SymmGroup>(parms, site_bases, initc, site_types));
    
    else if (parms["init_state"] == "basis_state")
        return initializer_ptr(new basis_mps_init<Matrix, SymmGroup>(parms, site_bases, site_types));
    
    else if (parms["init_state"] == "basis_state_generic")
        return initializer_ptr(new basis_mps_init_generic<Matrix, SymmGroup>(parms, site_bases, initc, site_types));
    
    else if (parms["init_state"] == "coherent")
        return initializer_ptr(new coherent_mps_init<Matrix, SymmGroup>(parms, site_bases, site_types));
    
    else if (parms["init_state"] == "basis_state_dm")
        return initializer_ptr(new basis_dm_mps_init<Matrix, SymmGroup>(parms, site_bases, site_types));
    
    else if (parms["init_state"] == "coherent_dm")
        return initializer_ptr(new coherent_dm_mps_init<Matrix, SymmGroup>(parms, site_bases, site_types));
    
    else {
        throw std::runtime_error("Don't know this initial state.");
        return initializer_ptr();
    }

}

#endif
