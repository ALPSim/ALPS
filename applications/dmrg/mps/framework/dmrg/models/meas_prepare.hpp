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

#ifndef MAQUIS_DMRG_MODEL_MEAS_PREPARE_HPP
#define MAQUIS_DMRG_MODEL_MEAS_PREPARE_HPP

#include "dmrg/models/lattice.h"
#include "dmrg/models/generate_mpo.hpp"

namespace meas_prepare {
    
    template<class Matrix, class SymmGroup>
    std::vector<std::pair<std::string, MPO<Matrix,SymmGroup> > >
    local(const Lattice & lat,
          std::vector<block_matrix<Matrix, SymmGroup> > const & identities,
          std::vector<block_matrix<Matrix, SymmGroup> > const & fillings,
          std::vector<std::vector<std::pair<std::vector<block_matrix<Matrix, SymmGroup> >, bool> > > const & ops)
    {
        std::vector<std::pair<std::string, MPO<Matrix,SymmGroup> > > mpos;
        if (ops.size() == 0)
            return mpos;
        
        for (std::size_t p = 0; p < lat.size(); ++p)
        {
            if (ops[0].size() == 1) {
                int type = lat.get_prop<int>("type", p);
                bool is_empty = true;
                generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat, identities, fillings);
                for (std::size_t i=0; i<ops.size(); ++i) {
                    if (ops[i][0].first[type].n_blocks() > 0) {
                        generate_mpo::Operator_Term<Matrix, SymmGroup> term;
                        term.operators.push_back( std::make_pair(p, ops[i][0].first[type]) );
                        mpom.add_term(term);
                        is_empty = false;
                    }
                }
                if (!is_empty) mpos.push_back(std::make_pair( lat.get_prop<std::string>("label", p) , mpom.create_mpo() ));
            } else {
                std::vector<Lattice::pos_t> neighs = lat.forward(p);
                for (typename std::vector<Lattice::pos_t>::const_iterator hopto = neighs.begin();
                     hopto != neighs.end();
                     ++hopto)
                {
                    int type1 = lat.get_prop<int>("type", p);
                    int type2 = lat.get_prop<int>("type", *hopto);
                    bool is_empty = true;
                    generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat, identities, fillings);
                    for (std::size_t i=0; i<ops.size(); ++i) {
                        if (ops[i][0].first[type1].n_blocks() > 0 && ops[i][1].first[type2].n_blocks() > 0) {
                            generate_mpo::Operator_Term<Matrix, SymmGroup> term;
                            term.operators.push_back( std::make_pair(p, ops[i][0].first[type1]) );
                            term.operators.push_back( std::make_pair(*hopto, ops[i][1].first[type2]) );
                            
                            block_matrix<Matrix, SymmGroup> m;
                            if (ops[i][0].second && p <= *hopto) {
                                gemm(fillings[type1], term.operators[0].second, m); // Note inverse notation because of notation in operator.
                                swap(term.operators[0].second, m);
                            } else if (ops[i][0].second) {
                                gemm(fillings[type2], term.operators[1].second, m); // Note inverse notation because of notation in operator.
                                swap(term.operators[0].second, m);
                            }
                            
                            term.with_sign = ops[i][0].second;
                            mpom.add_term(term);
                            is_empty = false;
                        }
                    }
                    if (!is_empty) mpos.push_back(std::make_pair( lat.get_prop<std::string>("label", p, *hopto) , mpom.create_mpo() ));
                }
            }
        }
        
        return mpos;
    }
    

    
	template<class Matrix, class SymmGroup>
	MPO<Matrix, SymmGroup>
    average(const Lattice & lat,
            std::vector<block_matrix<Matrix, SymmGroup> > const & identities,
            std::vector<block_matrix<Matrix, SymmGroup> > const & fillings,
            std::vector<std::pair<std::vector<block_matrix<Matrix, SymmGroup> >, bool> > const & ops)
	{
        generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat, identities, fillings);
        
        for (std::size_t p = 0; p < lat.size(); ++p)
        {
            generate_mpo::Operator_Term<Matrix, SymmGroup> term;
            term.operators.push_back( std::make_pair(p, ops[0].first[lat.get_prop<int>("type", p)]) );
            if (ops.size() == 1) {
                mpom.add_term(term);
            } else {
                term.with_sign = ops[0].second;
            	std::vector<Lattice::pos_t> neighs = lat.forward(p);
            	for (typename std::vector<Lattice::pos_t>::const_iterator hopto = neighs.begin();
            		 hopto != neighs.end();
            		 ++hopto)
            	{
                    int type1 = lat.get_prop<int>("type", p);
                    int type2 = lat.get_prop<int>("type", *hopto);

                    generate_mpo::Operator_Term<Matrix, SymmGroup> term2(term);
                    term2.operators.push_back( std::make_pair(*hopto, ops[1].first[type1]) );
                    
                    block_matrix<Matrix, SymmGroup> m;
                    if (ops[0].second && p <= *hopto) {
                        gemm(fillings[type1], term.operators[0].second, m); // Note inverse notation because of notation in operator.
                        swap(term.operators[0].second, m);
                    } else if (ops[0].second) {
                        gemm(fillings[type2], term.operators[1].second, m); // Note inverse notation because of notation in operator.
                        swap(term.operators[0].second, m);
                    }

                    mpom.add_term(term2);
            	}
                
            }
        }
        
        return mpom.create_mpo();
    }
    
}

#endif
