/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2016 Institute for Theoretical Physics, ETH Zurich
 *               2011-2016 by Michele Dolfi <dolfim@phys.ethz.ch>
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

