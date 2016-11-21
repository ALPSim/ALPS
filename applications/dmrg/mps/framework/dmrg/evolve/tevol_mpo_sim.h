/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef APP_DMRG_TEVOL_MPO_SIM_H
#define APP_DMRG_TEVOL_MPO_SIM_H

#include <cmath>

#include "dmrg/utils/storage.h"
#include "dmrg/evolve/te_utils.hpp"
#include "dmrg/evolve/trotter_decomposer.h"
#include "dmrg/mp_tensors/mpo_contractor.h"
#include "dmrg/utils/results_collector.h"

// ******   SIMULATION CLASS   ******
template <class Matrix, class SymmGroup>
class mpo_evolver {
    typedef term_descriptor<typename Matrix::value_type> term_t;
    typedef trotter_decomposer::steps_iterator steps_iterator;
    typedef trotter_decomposer::sequence_type sequence_type;
public:
    mpo_evolver(DmrgParameters * parms_, MPS<Matrix, SymmGroup> * mps_,
                Lattice const& lattice_, Model<Matrix, SymmGroup> const& model_,
                int init_sweep=0)
    : parms(parms_)
    , mps(mps_)
    , lattice(lattice_) // shallow copy
    , model(model_) // shallow copy
    , hamils(separate_hamil_terms(model.hamiltonian_terms()))
    , decomp(hamils.size(),
             (*parms)["te_order"],
             (*parms)["te_optim"])
    {
        maquis::cout << "Using MPO time evolution." << std::endl;
        
        maquis::cout << "Found " << hamils.size() << " non overlapping Hamiltonians." << std::endl;
        maquis::cout << "Using " << decomp.description() << std::endl;

        /// compute the time evolution gates
        prepare_te_terms(init_sweep);
    }
    
    void prepare_te_terms(unsigned sweep)
    {
        double dt = (*parms)["dt"];
        typename Matrix::value_type I;
        if (sweep < (*parms)["nsweeps_img"])
            I = maquis::traits::real_identity<typename Matrix::value_type>::value;
        else
            I = maquis::traits::imag_identity<typename Matrix::value_type>::value;
        typename Matrix::value_type alpha = -I*dt;
        
        Uterms.resize(decomp.size());
        for (int i=0; i<decomp.size(); ++i)
            Uterms[i] = make_exp_mpo(lattice, model, hamils[decomp.trotter_term(i).first], decomp.trotter_term(i).second*alpha);
    }
    
    void operator()(unsigned sweep, unsigned nsteps)
    {
        iteration_results_.clear();
        for (steps_iterator it=decomp.steps_begin(nsteps); it != steps_iterator(); ++it)
            evolve_time_step(sweep+it.index(), *it);
    }
    
    results_collector const& iteration_results() const
    {
        return iteration_results_;
    }
    
private:
    void evolve_time_step(unsigned sweep, sequence_type const & terms_sequence)
    {
        for (int i=0; i<terms_sequence.size(); ++i) {
            int which = terms_sequence[i];
            
            int maxiter = 6; double tol = 1e-6;
            mpo_contractor<Matrix, SymmGroup, storage::nop> evolution(*mps, Uterms[which], (*parms));
            for (int k = 0; k < maxiter; ++k) {
                std::pair<double,double> eps = evolution.sweep(sweep);
                double rel_error = std::abs( (eps.first-eps.second) / eps.second );
                if (rel_error < tol)
                    break;
                if (k == maxiter-1)
                    std::cerr << "Accuracy of variational compression reached only " << rel_error << " (maxiter="<< maxiter << ", tol=" << tol << ")" << std::endl;
            }
            evolution.finalize();
            *mps = evolution.get_current_mps();
        }
    }

    
private:        
    DmrgParameters * parms;
    MPS<Matrix, SymmGroup> * mps;
    Lattice lattice;
    Model<Matrix,SymmGroup> model;
    
    results_collector iteration_results_;
    
    std::vector<std::vector<term_t> > hamils;
    trotter_decomposer decomp;
    std::vector<MPO<Matrix, SymmGroup> > Uterms;
};

#endif
