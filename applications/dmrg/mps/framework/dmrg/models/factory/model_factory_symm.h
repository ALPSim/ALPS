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

#include "matrices.h"

#include "dmrg/models/coded/factory.h"
#include "dmrg/models/factory/initializer_factory.h"

#include "dmrg/models/alps/model.hpp"


// init MACROS
#define impl_model_factory(MATRIX, SYMMGROUP)                                           \
template boost::shared_ptr<model_impl<MATRIX, SYMMGROUP> >                              \
model_factory<MATRIX,SYMMGROUP>(Lattice const&, BaseParameters &);                      


// Implementation
template <class Matrix, class SymmGroup>
boost::shared_ptr<model_impl<Matrix, SymmGroup> >
model_factory(Lattice const& lattice, BaseParameters & parms)
{
    typedef boost::shared_ptr<model_impl<Matrix, SymmGroup> > impl_ptr;
    if (parms["model_library"] == "alps") {
        if (parms["lattice_library"] != "alps")
            throw std::runtime_error("ALPS models require ALPS lattice.");
        return impl_ptr( new ALPSModel<Matrix, SymmGroup>(lattice, parms) );
    } else if (parms["model_library"] == "coded") {
        return coded_model_factory<Matrix, SymmGroup>::parse(lattice, parms);
    } else {
        throw std::runtime_error("Don't know this model_library!");
    }
    
}

