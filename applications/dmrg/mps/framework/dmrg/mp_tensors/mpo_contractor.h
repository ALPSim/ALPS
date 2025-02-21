/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2016 Institute for Theoretical Physics, ETH Zurich
 *               2011-2016 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MP_TENSORS_MPO_CONTRACTOR_H
#define MP_TENSORS_MPO_CONTRACTOR_H

#include "dmrg/mp_tensors/mpo_contractor_ss.h"
#include "dmrg/mp_tensors/mpo_contractor_ts.h"


template<class Matrix, class SymmGroup, class Storage>
class mpo_contractor
{
public:
    mpo_contractor(MPS<Matrix, SymmGroup> const & mps,
                   MPO<Matrix, SymmGroup> const & mpo,
                   BaseParameters & parms)
    {
        /// Contractor factory
        if (parms["mpo_compression"] == "singlesite")
            impl_.reset(new mpo_contractor_ss<Matrix, SymmGroup, Storage>(mps, mpo, parms));
        else if (parms["mpo_compression"] == "twosite")
            impl_.reset(new mpo_contractor_ts<Matrix, SymmGroup, Storage>(mps, mpo, parms));
        else
            throw std::runtime_error("Do no know mpo_compression="+parms["mpo_compression"].str());
    }
    
    std::pair<double,double> sweep(int sweep)
    {
        return impl_->sweep(sweep);
    }
    void finalize()
    {
        return impl_->finalize();
    }
    MPS<Matrix, SymmGroup> const& get_original_mps() const
    {
        return impl_->get_original_mps();
    }
    MPS<Matrix, SymmGroup> const& get_current_mps() const
    {
        return impl_->get_current_mps();
    }
    
private:
    boost::shared_ptr<mpo_contractor_base<Matrix, SymmGroup, Storage> > impl_;
};

#endif

