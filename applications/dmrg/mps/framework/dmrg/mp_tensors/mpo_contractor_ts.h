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

#ifndef MP_TENSORS_MPO_CONTRACTOR_TS_H
#define MP_TENSORS_MPO_CONTRACTOR_TS_H

#include "dmrg/mp_tensors/mpo_contractor_base.h"

#include "dmrg/mp_tensors/twositetensor.h"
#include "dmrg/mp_tensors/mpo_ops.h"


/// TODO: 1) implement zip-up compression. E. M. Stoudenmire and S. R. White, New Journal of Physics 12, 055026 (2010).

template<class Matrix, class SymmGroup, class Storage>
class mpo_contractor_ts : public mpo_contractor_base<Matrix, SymmGroup, Storage>
{
    typedef mpo_contractor_base<Matrix, SymmGroup, Storage> base;
    using base::mpo;
    using base::mps;
    using base::mpsp;
    using base::left_;
    using base::right_;
    using base::parms;

public:
    mpo_contractor_ts(MPS<Matrix, SymmGroup> & mps_,
                      MPO<Matrix, SymmGroup> const & mpo_,
                      BaseParameters & parms_)
    : base(mps_, mpo_, parms_)
    {
        make_ts_cache_mpo(mpo, ts_cache_mpo, mps);
    }
    
    inline int to_site(const int L, const int i) const
    {
        if (i < 0) return 0;
        /// i, or (L-1) - (i - (L-1))
        return (i < L-1) ? i : 2*L - 2 - i;
    }

    std::pair<double,double> sweep(int sweep)
    {
        boost::chrono::high_resolution_clock::time_point sweep_now = boost::chrono::high_resolution_clock::now();
        
        std::size_t L = mps.length();
        
        std::pair<double,double> eps;
        block_matrix<Matrix, SymmGroup> norm_boudary;
        norm_boudary.insert_block(Matrix(1, 1, 1), SymmGroup::IdentityCharge, SymmGroup::IdentityCharge);
    
        
        int _site = 0, site = 0;
        
        if (_site < L-1) {
            Storage::prefetch(left_[site]);
            Storage::prefetch(right_[site+2]);
        } else {
            Storage::prefetch(left_[site-1]);
            Storage::prefetch(right_[site+1]);
        }
        
        for (; _site < 2*L-2; ++_site) {
            /* (0,1), (1,2), ... , (L-1,L), (L-1,L), (L-2, L-1), ... , (0,1)
             | |                        |
             site 1                      |
             |         left to right  | right to left, lr = -1
             site 2                   |                               */
            
            int lr, site1, site2;
            if (_site < L-1) {
                site = to_site(L, _site);
                lr = 1;
                site1 = site;
                site2 = site+1;
                ts_cache_mpo[site1].placement_l = mpo[site1].placement_l;
                ts_cache_mpo[site1].placement_r = get_right_placement(ts_cache_mpo[site1], mpo[site1].placement_l, mpo[site2].placement_r);
            } else {
                site = to_site(L, _site);
                lr = -1;
                site1 = site-1;
                site2 = site;
                ts_cache_mpo[site1].placement_l = get_left_placement(ts_cache_mpo[site1], mpo[site1].placement_l, mpo[site2].placement_r);
                ts_cache_mpo[site1].placement_r = mpo[site2].placement_r;
            }
            
//            maquis::cout << std::endl;
//            maquis::cout << "Sweep " << sweep << ", optimizing sites " << site1 << " and " << site2 << std::endl;
            
            if (_site != L-1)
            {
                Storage::fetch(left_[site1]);
                Storage::fetch(right_[site2+1]);
            }
            
            if (lr == +1) {
                if (site2+2 < right_.size()){
                    Storage::prefetch(right_[site2+2]);
                }
            } else {
                if (site1 > 0){
                    Storage::prefetch(left_[site1-1]);
                }
            }
            
            
            boost::chrono::high_resolution_clock::time_point now, then;
            
            /// Create TwoSite objects
            /// input
            TwoSiteTensor<Matrix, SymmGroup> tst(mps[site1], mps[site2]);
            MPSTensor<Matrix, SymmGroup> twin_mps = tst.make_mps();
            /// output
            TwoSiteTensor<Matrix, SymmGroup> tstp(mpsp[site1], mpsp[site2]);
            MPSTensor<Matrix, SymmGroup> twin_mpsp = tstp.make_mps();
            
            
            /// Perform multiplication with MPO
            SiteProblem<Matrix, SymmGroup> sp(left_[site1], right_[site2+1], ts_cache_mpo[site1]);
            ietl::mult(sp, twin_mps, twin_mpsp);
            tstp << twin_mpsp;

            
            double cutoff = this->get_cutoff(sweep);
            std::size_t Mmax = this->get_Mmax(sweep);
            truncation_results trunc;
            
            if (lr == +1)
            {
                // Write back result from optimization
                boost::tie(mpsp[site1], mpsp[site2], trunc) = tstp.split_mps_l2r(Mmax, cutoff);
                
//                block_matrix<Matrix, SymmGroup> t;
                
//                if (site2 < L-1) mpsp[site2+1].multiply_from_left(t);
                
                if (site1 != L-2)
                    Storage::drop(right_[site2+1]);
                
                this->boundary_left_step(mpo, site1); // creating left_[site2]
                norm_boudary = contraction::overlap_left_step(mpsp[site1], MPSTensor<Matrix,SymmGroup>(mpsp[site1]), norm_boudary);

                
                if (site1 != L-2){
                    Storage::evict(mps[site1]);
                    Storage::evict(mpsp[site1]);
                    Storage::evict(left_[site1]);
                }
            }
            if (lr == -1){
                // Write back result from optimization
                boost::tie(mpsp[site1], mpsp[site2], trunc) = tstp.split_mps_r2l(Mmax, cutoff);
                
//                block_matrix<Matrix, SymmGroup> t;
                
//                if (site1 > 0) mpsp[site1-1].multiply_from_right(t);
                
                if(site1 != 0)
                    Storage::drop(left_[site1]);
                
                this->boundary_right_step(mpo, site2); // creating right_[site2]
                norm_boudary = contraction::overlap_right_step(mpsp[site2], MPSTensor<Matrix,SymmGroup>(mpsp[site2]), norm_boudary);

                
                if(site1 != 0){
                    Storage::evict(mps[site2]);
                    Storage::evict(mpsp[site2]);
                    Storage::evict(right_[site2+1]);
                }
            }

            
            /// Finalize forward sweep
            if (_site == L-2) {
                norm_boudary = contraction::overlap_left_step(mpsp[site2], MPSTensor<Matrix,SymmGroup>(mpsp[site2]), norm_boudary);
                this->boundary_left_step(mpo, site2); // creating left_[site2+1]

                double nn = maquis::real( norm_boudary.trace() );
                eps.first = nn - 2.*maquis::real(left_[L][0].trace());
                
                /// prepare backward sweep
                norm_boudary = block_matrix<Matrix, SymmGroup>();
                norm_boudary.insert_block(Matrix(1, 1, 1), mps[L-1].col_dim()[0].first, mps[L-1].col_dim()[0].first);
            }
            
            /// Finalize backward sweep
            if (_site == 2*L-3) {
                norm_boudary = contraction::overlap_right_step(mpsp[site1], MPSTensor<Matrix,SymmGroup>(mpsp[site1]), norm_boudary);
                this->boundary_right_step(mpo, site1); // creating right_[site1]

                double nn = maquis::real( norm_boudary.trace() );
                eps.second = nn - 2.*maquis::real(right_[0][0].trace());
            }
            
        }
        
        return eps; /// note: the actual eps contain a constant, which is not important here.
    }
    
private:
    MPO<Matrix, SymmGroup> ts_cache_mpo;
};

#endif

