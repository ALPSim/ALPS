/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
 *               2012      by Jan Gukelberger <gukelberger@phys.ethz.ch>
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

#include "dmrg/models/coded/models_none.hpp"
#include "dmrg/models/coded/super_models_none.hpp"

template<class Matrix>
struct coded_model_factory<Matrix, TrivialGroup> {
    static boost::shared_ptr<model_impl<Matrix, TrivialGroup> > parse
    (Lattice const& lattice, BaseParameters & parms)
    {
        typedef boost::shared_ptr<model_impl<Matrix, TrivialGroup> > impl_ptr;
        if (parms["MODEL"] == std::string("boson Hubbard"))
            return impl_ptr( new BoseHubbardNone<Matrix>(lattice, parms) );
        else if (parms["MODEL"] == std::string("super boson Hubbard"))
            return impl_ptr( new SuperBoseHubbardNone<Matrix>(lattice, parms) );
        else {
            throw std::runtime_error("Don't know this model with NONE symmetry group!");
            return impl_ptr();
        }
    }
};
