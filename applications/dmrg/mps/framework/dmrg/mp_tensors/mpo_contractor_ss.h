/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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

#ifndef MP_TENSORS_MPO_CONTRACTOR_SS_H
#define MP_TENSORS_MPO_CONTRACTOR_SS_H

#include "dmrg/mp_tensors/mpo_contractor_base.h"


/// TODO: 1) implement two-site time evolution. (single-site is stuck in initial MPS structure)
///       2) implement zip-up compression. E. M. Stoudenmire and S. R. White, New Journal of Physics 12, 055026 (2010).

template<class Matrix, class SymmGroup, class Storage>
class mpo_contractor_ss : public mpo_contractor_base<Matrix, SymmGroup, Storage>
{
    typedef mpo_contractor_base<Matrix, SymmGroup, Storage> base;
    using base::mpo;
    using base::mps;
    using base::mpsp;
    using base::left_;
    using base::right_;
    using base::parms;

public:
    mpo_contractor_ss(MPS<Matrix, SymmGroup> & mps_,
                      MPO<Matrix, SymmGroup> const & mpo_,
                      BaseParameters & parms_)
    : base(mps_, mpo_, parms_)
    { }
    
    std::pair<double,double> sweep(int sweep)
    {
        boost::chrono::high_resolution_clock::time_point sweep_now = boost::chrono::high_resolution_clock::now();
        
        std::size_t L = mps.length();
        
        std::pair<double,double> eps;
        block_matrix<Matrix, SymmGroup> norm_boudary;
        norm_boudary.insert_block(Matrix(1, 1, 1), SymmGroup::IdentityCharge, SymmGroup::IdentityCharge);
    
        for (int _site = 0; _site < 2*L; ++_site)
        {
            int site, lr;
            if (_site < L) {
                site = _site;
                lr = 1;
            } else {
                site = 2*L-_site-1;
                lr = -1;
            }
            
            SiteProblem<Matrix, SymmGroup> sp(left_[site], right_[site+1], mpo[site]);
            ietl::mult(sp, mps[site], mpsp[site]);
            
            if (lr == +1) {
                if (site < L-1) {
                    block_matrix<Matrix, SymmGroup> t;
                    t = mpsp[site].normalize_left(DefaultSolver());
                    mpsp[site+1].multiply_from_left(t);
                }
                
                this->boundary_left_step(mpo, site); // creating left_[site+1]
                norm_boudary = contraction::overlap_left_step(mpsp[site], MPSTensor<Matrix,SymmGroup>(mpsp[site]), norm_boudary);
            } else if (lr == -1) {
                if (site > 0) {
                    block_matrix<Matrix, SymmGroup> t;
                    t = mpsp[site].normalize_right(DefaultSolver());
                    mpsp[site-1].multiply_from_right(t);
                }   
                
                this->boundary_right_step(mpo, site); // creating right_[site]
                norm_boudary = contraction::overlap_right_step(mpsp[site], MPSTensor<Matrix,SymmGroup>(mpsp[site]), norm_boudary);
            }
            
            if (_site == L-1) {
                double nn = maquis::real( norm_boudary.trace() );
                eps.first = nn - 2.*maquis::real(left_[L][0].trace());
                
                /// prepare backward sweep
                norm_boudary = block_matrix<Matrix, SymmGroup>();
                norm_boudary.insert_block(Matrix(1, 1, 1), mps[L-1].col_dim()[0].first, mps[L-1].col_dim()[0].first);
            }
            
            if (_site == 2*L-1) {
                double nn = maquis::real( norm_boudary.trace() );
                eps.second = nn - 2.*maquis::real(right_[0][0].trace());
            }
            
        }
        
        return eps; /// note: the actual eps contain a constant, which is not important here.
    }
    
};

#endif

