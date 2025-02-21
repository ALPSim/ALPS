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

#include "dmrg/models/coded/models_u1.hpp"
//#include "dmrg/models/coded/models_bela.hpp"

template<class Matrix>
struct coded_model_factory<Matrix, U1> {
    static boost::shared_ptr<model_impl<Matrix, U1> > parse
    (Lattice const& lattice, BaseParameters & parms)
    {
        typedef boost::shared_ptr<model_impl<Matrix, U1> > impl_ptr;
        if (parms["MODEL"] == std::string("heisenberg"))
            return impl_ptr( new Heisenberg<Matrix>(lattice, parms["Jxy"], parms["Jz"]) );
        else if (parms["MODEL"] == std::string("HCB"))
            return impl_ptr( new HCB<Matrix>(lattice) );
        else if (parms["MODEL"] == std::string("boson Hubbard"))
            return impl_ptr( new BoseHubbard<Matrix>(lattice, parms) );
//        else if (parms["MODEL"] == std::string("fermion Hubbard"))
//            return impl_ptr( new FermiHubbardU1<Matrix>(lattice, parms) );
        else if (parms["MODEL"] == std::string("FreeFermions"))
            return impl_ptr( new FreeFermions<Matrix>(lattice, parms["t"]) );
//        else if (parms["MODEL"] == std::string("bela_chiral"))
//            return impl_ptr( new Chiral<Matrix>(lattice, parms) );
        else {
            throw std::runtime_error("Don't know this model!");
            return impl_ptr();
        }
    }
};
