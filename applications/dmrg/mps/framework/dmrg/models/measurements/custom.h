/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MEASUREMENTS_CUSTOM_H
#define MEASUREMENTS_CUSTOM_H

#include "dmrg/models/measurement.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"

namespace measurements {
    
    template <class Matrix, class SymmGroup>
    class custom : public measurement<Matrix, SymmGroup> {
        typedef  measurement<Matrix, SymmGroup> base;
    public:
        custom(std::string const& name_,
               const Lattice & lat,
               std::vector<block_matrix<Matrix, SymmGroup> > const & identities,
               std::vector<block_matrix<Matrix, SymmGroup> > const & fillings,
               std::vector< std::vector< std::pair<int, block_matrix<Matrix, SymmGroup> > > > const & ops)
        : base(name_)
        {
            generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat, identities, fillings);
            
            for (int k = 0; k < ops.size(); ++k) {
                generate_mpo::Operator_Term<Matrix, SymmGroup> term;
                term.operators = ops[k];
                term.with_sign = true;
                mpom.add_term(term);
            }
            
            mpo = mpom.create_mpo();
            
            this->cast_to_real = false;
        }
        
        void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
        {
            this->result = expval(mps, mpo);
        }
    
    protected:
        measurement<Matrix, SymmGroup>* do_clone() const
        {
            return new custom(*this);
        }
        
    private:
        MPO<Matrix, SymmGroup> mpo;
    };
    
}

#endif
