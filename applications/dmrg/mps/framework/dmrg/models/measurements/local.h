/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MEASUREMENTS_LOCAL_H
#define MEASUREMENTS_LOCAL_H

#include "dmrg/models/measurement.h"
#include "dmrg/models/meas_prepare.hpp"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/super_mpo.h"

namespace measurements {
    
    template <class Matrix, class SymmGroup>
    class local : public measurement<Matrix, SymmGroup> {
        typedef measurement<Matrix, SymmGroup> base;
        typedef generate_mpo::MPOMaker<Matrix, SymmGroup> generator;
        typedef std::vector<block_matrix<Matrix, SymmGroup> > op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
    public:
        
        local(std::string const& name_, const Lattice & lat,
              op_vec const & identities_, op_vec const & fillings_,
              std::vector<bond_element> const& terms)
        : base(name_)
        , lattice(lat)
        , identities(identities_)
        , fillings(fillings_)
        , is_bond(terms.size() > 0 && terms[0].size() > 1)
        , mpo_terms(terms)
        {
            this->cast_to_real = all_true(mpo_terms.begin(), mpo_terms.end(), static_cast<bool (*)(bond_element const&)>(&is_hermitian_meas));
        }

        void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
        {
            this->vector_results.clear();
            this->labels.clear();
            
            typedef typename SymmGroup::subcharge subcharge;
            if (!rmps || this->is_super_meas || is_bond)
                return evaluate_with_mpo(mps);
            
            int L = mps.size();
            this->vector_results.reserve(this->vector_results.size() + L);
            this->labels.reserve(this->labels.size() + L);
            
            for (typename Lattice::pos_t p = 0; p < L; ++p) {
                subcharge type = lattice.get_prop<subcharge>("type", p);
                typename Matrix::value_type res = 0.;
                bool evaluated = false;
                for (std::size_t i=0; i<mpo_terms.size(); ++i) {
                    if (mpo_terms[i][0].first[type].n_blocks() > 0) {
                        MPOTensor<Matrix, SymmGroup> temp;
                        temp.set(0, 0, mpo_terms[i][0].first[type]);
                        
                        MPSTensor<Matrix, SymmGroup> vec2 = contraction::site_hamil2(mps[p], rmps.get().left(p), rmps.get().right(p), temp);
                        res += mps[p].scalar_overlap(vec2);
                        evaluated = true;
                    }
                }
                if (evaluated) {
                    this->vector_results.push_back(res);
                    this->labels.push_back( lattice.get_prop<std::string>("label", p) );
                }
            }
        }
        
    protected:
        measurement<Matrix, SymmGroup>* do_clone() const
        {
            return new local(*this);
        }
        
        void evaluate_with_mpo(MPS<Matrix, SymmGroup> const& mps)
        {
            
            typename MPS<Matrix, SymmGroup>::scalar_type nn;
            if (this->is_super_meas)
                nn = dm_trace(mps, this->phys_psi);
            
            typedef std::vector<std::pair<std::string, MPO<Matrix,SymmGroup> > > local_mpo_list_t;
            local_mpo_list_t mpos = meas_prepare::local(lattice, identities, fillings, mpo_terms);
            for (typename local_mpo_list_t::const_iterator it = mpos.begin(); it != mpos.end(); ++it) {
                if (!this->is_super_meas) {
                    typename Matrix::value_type val = expval(mps, it->second);
                    this->vector_results.push_back((this->cast_to_real) ? maquis::real(val) : val);
                    this->labels.push_back(it->first);
                } else {
                    MPS<Matrix, SymmGroup> super_mpo = mpo_to_smps(it->second, this->phys_psi);
                    // static_cast needed for icpc 12.x
                    typedef typename MPS<Matrix, SymmGroup>::scalar_type (*overlap_func)(MPS<Matrix, SymmGroup> const &, MPS<Matrix, SymmGroup> const &);
                    typename MPS<Matrix, SymmGroup>::scalar_type val = static_cast<overlap_func>(::overlap)(super_mpo, mps);
                    val = val/nn;
                    
                    this->vector_results.push_back((this->cast_to_real) ? maquis::real(val) : val);
                    this->labels.push_back(it->first);
                }
            }
        }
        
    private:
        Lattice lattice;
        op_vec identities, fillings;
        bool is_bond;
        std::vector<bond_element> mpo_terms;
    };
    
}

#endif
