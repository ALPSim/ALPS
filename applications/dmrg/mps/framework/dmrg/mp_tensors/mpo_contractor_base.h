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

#ifndef MP_TENSORS_MPO_CONTRACTOR_BASE_H
#define MP_TENSORS_MPO_CONTRACTOR_BASE_H

#include "dmrg/optimize/optimize.h"

#include "dmrg/utils/BaseParameters.h"
#include "dmrg/utils/placement.h"



template<class Matrix, class SymmGroup, class Storage>
class mpo_contractor_base
{
public:
    mpo_contractor_base(MPS<Matrix, SymmGroup> const & mps_,
                        MPO<Matrix, SymmGroup> const & mpo_,
                        BaseParameters & parms_)
    : mps(mps_)
    , mpsp(mps_)
    , mpo(mpo_)
    , parms(parms_)
    {
        init_left_right(mpo, 0);
    }
    
    virtual std::pair<double,double> sweep(int sweep) =0;
    
    virtual void finalize()
    {
        mpsp[0].normalize_right(DefaultSolver());
    }
    
    virtual MPS<Matrix, SymmGroup> const& get_original_mps() const { return mps; }
    virtual MPS<Matrix, SymmGroup> const& get_current_mps() const { return mpsp; }
    
    virtual ~mpo_contractor_base() {}

protected:
    inline void boundary_left_step(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        left_[site+1] = contraction::overlap_mpo_left_step(mpsp[site], mps[site], left_[site], mpo[site]);
        Storage::pin(left_[site+1]);
    }
    
    inline void boundary_right_step(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        right_[site] = contraction::overlap_mpo_right_step(mpsp[site], mps[site], right_[site+1], mpo[site]);
        Storage::pin(right_[site]);
    }
    
    void init_left_right(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        construct_placements(mpo);
        std::size_t L = mps.length();
        
        left_.resize(mpo.length()+1);
        right_.resize(mpo.length()+1);
        
        Storage::drop(left_[0]);
        left_[0] = mps.left_boundary();
        Storage::pin(left_[0]);
        
        for (int i = 0; i < site; ++i) {
            Storage::drop(left_[i+1]);
            boundary_left_step(mpo, i);
            Storage::evict(left_[i]);
        }
        Storage::evict(left_[site]);
        
#ifndef NDEBUG
        maquis::cout << "Boundaries are partially initialized...\n";
#endif
        
        Storage::drop(right_[L]);
        right_[L] = mps.right_boundary();
        Storage::pin(right_[L]);
        
        for (int i = L-1; i >= site; --i) {
            Storage::drop(right_[i]);
            boundary_right_step(mpo, i);
            Storage::evict(right_[i+1]);
        }
        Storage::evict(right_[site]);
        
#ifndef NDEBUG
        maquis::cout << "Boundaries are fully initialized...\n";
#endif
    }
    
    double get_cutoff(int sweep) const
    {
        double cutoff;
        if (sweep >= parms.template get<int>("ngrowsweeps"))
            cutoff = parms.template get<double>("truncation_final");
        else
            cutoff = log_interpolate(parms.template get<double>("truncation_initial"), parms.template get<double>("truncation_final"), parms.template get<int>("ngrowsweeps"), sweep);
        return cutoff;
    }
    
    std::size_t get_Mmax(int sweep) const
    {
        std::size_t Mmax;
        if (parms.is_set("sweep_bond_dimensions")) {
            std::vector<std::size_t> ssizes = parms.template get<std::vector<std::size_t> >("sweep_bond_dimensions");
            if (sweep >= ssizes.size())
                Mmax = *ssizes.rbegin();
            else
                Mmax = ssizes[sweep];
        } else
            Mmax = parms.template get<std::size_t>("max_bond_dimension");
        return Mmax;
    }
    
    MPS<Matrix, SymmGroup> const& mps;
    MPS<Matrix, SymmGroup> mpsp;
    MPO<Matrix, SymmGroup> const& mpo;
    
    BaseParameters & parms;
    std::vector<Boundary<typename storage::constrained<Matrix>::type, SymmGroup> > left_, right_;
};

#endif

