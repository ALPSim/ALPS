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

#ifndef MEASUREMENTS_AVERAGE_H
#define MEASUREMENTS_AVERAGE_H

#include "dmrg/models/measurement.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/super_mpo.h"

namespace measurements {

    template <class Matrix, class SymmGroup>
    class average : public measurement<Matrix, SymmGroup> {
        typedef  measurement<Matrix, SymmGroup> base;
        typedef generate_mpo::MPOMaker<Matrix, SymmGroup> generator;
        typedef std::vector<block_matrix<Matrix, SymmGroup> > op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
    public:

        average(std::string const& name_,
                const Lattice & lattice,
                op_vec const & identities, op_vec const & fillings,
                std::vector<bond_element> const& ops)
        : base(name_)
        {
            generator mpom(lattice, identities, fillings);

            for (std::size_t i = 0; i < ops.size(); ++i) {
                for (std::size_t p = 0; p < lattice.size(); ++p) {
                    int type1 = lattice.get_prop<int>("type", p);

                    generate_mpo::Operator_Term<Matrix, SymmGroup> term;
                    term.operators.push_back( std::make_pair(p, ops[i][0].first[type1]) );
                    if (ops[i].size() == 1) { // site term
                        term.with_sign = false;
                        mpom.add_term(term);

                    } else {
                        term.with_sign = ops[i][0].second;
                        std::vector<Lattice::pos_t> neighs = lattice.forward(p);
                        for (typename std::vector<Lattice::pos_t>::const_iterator hopto = neighs.begin();
                             hopto != neighs.end(); ++hopto)
                        {
                            int type2 = lattice.get_prop<int>("type", *hopto);

                            generate_mpo::Operator_Term<Matrix, SymmGroup> term2(term);
                            term2.operators.push_back( std::make_pair(*hopto, ops[i][1].first[type2]) );

                            block_matrix<Matrix, SymmGroup> m;
                            if (ops[i][0].second && p <= *hopto) {
                                gemm(fillings[type1], term2.operators[0].second, m); // Note inverse notation because of notation in operator.
                                swap(term2.operators[0].second, m);
                            } else if (ops[i][0].second) {
                                gemm(fillings[type2], term2.operators[1].second, m); // Note inverse notation because of notation in operator.
                                swap(term2.operators[0].second, m);
                            }

                            mpom.add_term(term2);
                        }
                    }
                }
            }
            mpo = mpom.create_mpo();

            /// TODO: this doesn't really work
            this->cast_to_real = all_true(ops.begin(), ops.end(), static_cast<bool (*)(bond_element const&)>(&is_hermitian_meas));
        }

        void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
        {
            if (!this->is_super_meas){
                this->result = (this->cast_to_real) ? maquis::real(expval(mps, mpo)) : expval(mps, mpo);
            } else {
                typename MPS<Matrix, SymmGroup>::scalar_type nn = dm_trace(mps, this->phys_psi);
                MPS<Matrix, SymmGroup> super_mpo = mpo_to_smps(mpo, this->phys_psi);
                // static_cast needed for icpc 12.x
#ifdef __INTEL_COMPILER
                typedef typename MPS<Matrix, SymmGroup>::scalar_type (*overlap_func)(MPS<Matrix, SymmGroup> const &, MPS<Matrix, SymmGroup> const &);
                this->result = static_cast<overlap_func>(&overlap)(super_mpo, mps) / nn;
#else
                this->result =::overlap<Matrix,SymmGroup>(super_mpo, mps) / nn;
#endif
            }
        }

    protected:
        measurement<Matrix, SymmGroup>* do_clone() const
        {
            return new average(*this);
        }

    private:
        MPO<Matrix, SymmGroup> mpo;
    };

}

#endif
