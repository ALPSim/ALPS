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

#include <boost/random.hpp>

#include "dmrg/optimize/ietl_lanczos_solver.h"
#include "dmrg/optimize/ietl_jacobi_davidson.h"
#ifdef HAVE_ARPACK
#include "dmrg/optimize/arpackpp_solver.h"
#endif

#include "dmrg/utils/BaseParameters.h"

#ifndef WIN32
#include <sys/time.h>
#define HAVE_GETTIMEOFDAY
#endif


template<class Matrix, class SymmGroup>
struct SiteProblem
{
    SiteProblem(MPSTensor<Matrix, SymmGroup> const & ket_tensor_,
                Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & left_,
                Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & right_,
                MPOTensor<Matrix, SymmGroup> const & mpo_)
    : ket_tensor(ket_tensor_)
    , left(left_)
    , right(right_)
    , mpo(mpo_) { }
    
    MPSTensor<Matrix, SymmGroup> const & ket_tensor;
    Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & left;
    Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & right;
    MPOTensor<Matrix, SymmGroup> const & mpo;
};

#ifdef HAVE_GETTIMEOFDAY
#define BEGIN_TIMING(name) \
gettimeofday(&now, NULL);
#define END_TIMING(name) \
gettimeofday(&then, NULL); \
maquis::cout << "Time elapsed in " << name << ": " << then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec) << std::endl;
#else
#define BEGIN_TIMING(name)
#define END_TIMING(name)
#endif

/// TODO: 1) implement two-site time evolution. (single-site is stuck in initial MPS structure)
///       2) implement zip-up compression. E. M. Stoudenmire and S. R. White, New Journal of Physics 12, 055026 (2010).

template<class Matrix, class SymmGroup, class Storage>
class mpo_contractor_ss
{
public:
    mpo_contractor_ss(MPS<Matrix, SymmGroup> const & mps_,
                      MPO<Matrix, SymmGroup> const & mpo_,
                      BaseParameters & parms_)
    : mps(mps_)
    , mpsp(mps_)
    , mpo(mpo_)
    , parms(parms_)
    {
        mps.canonize(0);
        init_left_right(mpo);
        
        mpsp = mps;
    }
    
    std::pair<double,double> sweep(int sweep)
    {
        timeval sweep_now, sweep_then;
        gettimeofday(&sweep_now, NULL);
        
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
            
            SiteProblem<Matrix, SymmGroup> sp(mps[site], left_[site], right_[site+1], mpo[site]);
            ietl::mult(sp, mps[site], mpsp[site]);
            
            if (lr == +1) {
                if (site < L-1) {
                    block_matrix<Matrix, SymmGroup> t;
                    t = mpsp[site].normalize_left(DefaultSolver());
                    mpsp[site+1].multiply_from_left(t);
                }
                
                left_[site+1] = contraction::overlap_mpo_left_step(mpsp[site], mps[site], left_[site], mpo[site]);
                norm_boudary = contraction::overlap_left_step(mpsp[site], MPSTensor<Matrix,SymmGroup>(mpsp[site]), norm_boudary);
            } else if (lr == -1) {
                if (site > 0) {
                    block_matrix<Matrix, SymmGroup> t;
                    t = mpsp[site].normalize_right(DefaultSolver());
                    mpsp[site-1].multiply_from_right(t);
                }   
                
                right_[site] = contraction::overlap_mpo_right_step(mpsp[site], mps[site], right_[site+1], mpo[site]);
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
    
    void finalize()
    {
        mpsp[0].normalize_right(DefaultSolver());
    }
    
    MPS<Matrix, SymmGroup> get_original_mps() const { return mps; }
    MPS<Matrix, SymmGroup> get_current_mps() const { return mpsp; }
    
private:
    void init_left_right(MPO<Matrix, SymmGroup> const & mpo)
    {
        std::size_t L = mps.length();
        
        left_.resize(mpo.length()+1);
        right_.resize(mpo.length()+1);
        
        Storage::drop(left_[0]);
        left_[0] = mps.left_boundary();
        Storage::evict(left_[0]);
        
        for (int i = 0; i < L; ++i) {
            Storage::drop(left_[i+1]);
            left_[i+1] = contraction::overlap_mpo_left_step(mpsp[i], mps[i], left_[i], mpo[i]);
            Storage::evict(left_[i+1]);
        }
        
        Storage::drop(right_[L]);
        right_[L] = mps.right_boundary();
        Storage::evict(right_[L]);
        
        for(int i = L-1; i >= 0; --i) {
            Storage::drop(right_[i]);
            right_[i] = contraction::overlap_mpo_right_step(mpsp[i], mps[i], right_[i+1], mpo[i]);
            Storage::evict(right_[i]);
        }
    }
    
    MPS<Matrix, SymmGroup> mps, mpsp;
    MPO<Matrix, SymmGroup> const& mpo;
    
    BaseParameters & parms;
    std::vector<Boundary<typename storage::constrained<Matrix>::type, SymmGroup> > left_, right_;
};

#endif
