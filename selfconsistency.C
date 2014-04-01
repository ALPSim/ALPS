/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005-2008 by Philipp Werner <werner@comp-phys.org>,
 *                            Emanuel Gull <gullc@comp-phys.org>
 *                            Sebastian Fuchs <fuchs@comp-phys.org>
 *                            Matthias Troyer <troyer@comp-phys.org>
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

/* $Id$ */

/// @file selfconsistency.C
/// @brief implements the selfconsistency loop functions

#include "selfconsistency.h"
#include "green_function.h"
#include "fouriertransform.h"
#include "types.h"
#include <sys/types.h> //will this fail on Windows?
#include <boost/tuple/tuple.hpp> //will this fail on Windows?

void selfconsistency_loop_single_site(const alps::Parameters& parms, MatsubaraImpuritySolver& solver,
                                FrequencySpaceHilbertTransformer& hilbert) 
{
  unsigned int n_tau=boost::lexical_cast<unsigned int>(parms["N"]);
  unsigned int n_matsubara=boost::lexical_cast<unsigned int>(parms["NMATSUBARA"]);
  unsigned int n_orbital=parms.value_or_default("FLAVORS", 2);
  unsigned int n_site=parms.value_or_default("SITES", 1);
  
  double beta = static_cast<double>(parms["BETA"]);
  double h = static_cast<double>(parms["H"]);
  double converged = static_cast<double>(parms["CONVERGED"]);
  bool paramagnet = parms.value_or_default("PARAMAGNET", false);
  //bool degenerate = parms.value_or_default("DEGENERATE", false);
  double relax_rate=static_cast<double>(parms.value_or_default("RELAX_RATE", 1.));
  int max_it=static_cast<int>(parms.value_or_default("MAX_IT", 1000));
  std::string basename=parms["BASENAME"];
  
  
  matsubara_green_function_t G0_omega = hilbert.initial_G0(parms);
  matsubara_green_function_t G_omega = G0_omega;
  
  int iteration_ctr = 0;
  //define multiple vectors
  matsubara_green_function_t G_old_omega(G0_omega);
  matsubara_green_function_t G0_old_omega(G0_omega);
  itime_green_function_t G_tau(n_tau +1, n_site, n_orbital, diagonal);
  itime_green_function_t G0_tau(n_tau+1, n_site, n_orbital, diagonal);
  itime_green_function_t G_old_tau(n_tau+1, n_site, n_orbital, diagonal);
  itime_green_function_t G0_old_tau(n_tau+1, n_site, n_orbital, diagonal);
  //initialize space for multiple vectors
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  FourierTransformer::generate_transformer(parms, fourier_ptr);
  double max_diff;	
  do {
    iteration_ctr++;
    std::cout<<"starting iteration nr. "<<iteration_ctr<<" H is: "<<h<<std::endl;
    double mu = static_cast<double>(parms["MU"]);
    G_old_omega = G_omega;
    G0_old_omega = G0_omega;
    G_old_tau = G_tau;
    G0_old_tau = G0_tau;
    fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega);
    //fourier.setmu(mu);
    std::cout<<"running solver."<<std::endl<<std::flush;
    boost::tie(G_omega, G_tau) = solver.solve_omega(G0_omega,parms);
    G_tau = hilbert.symmetrize(G_tau, paramagnet);
    G_omega = hilbert.symmetrize(G_omega, paramagnet);
    /*if(parms.defined("FILLING")){ //adjust chemical potential
     bool converged=false;
     matsubara_green_function_t G_lattice_omega(G_omega);
     do{
     //standard hilbert trafo
     matsubara_green_function_t G0_omega_new= hilbert(G_omega, G0_omega, mu, h, beta);
     //re-obtain lattice Green's function from it
     
     //compute density
     double total_density=compute_density(G_omega);
     //adjust chemical potential
     converged=adjust_chemical_potential(parms, total_density, mu);
     }while(!converged);
     }*/
    std::cout<<"running Hilbert transform"<<std::endl<<std::flush;
    G0_omega = hilbert(G_omega, G0_omega, mu, h, beta);
    //relaxation to speed up /slow down convergence
    if(relax_rate !=1){
      std::cout<<"using over/underrelaxation with rate: "<<relax_rate<<std::endl;
      for(unsigned int orbital=0;orbital<n_orbital;++orbital){
        for(unsigned int site1=0;site1<n_site;++site1){
          for(unsigned int site2=0;site2<n_site;++site2){
            for(unsigned int i=0;i<n_matsubara;++i){
              G0_omega(i, site1, site2, orbital)=relax_rate* G0_omega(i, site1, site2, orbital)+
              (1-relax_rate)*G0_old_omega(i, site1, site2, orbital);
            }
          }
        }
      }
    }
    std::cout<<"comparing old and new result."<<std::endl<<std::flush;
    max_diff=0;
    if(iteration_ctr>1){
      //comparison for the dressed Green's function in Matsubara freq.
      for(unsigned int w=0; w<n_matsubara; w++) {
        for(unsigned int i1=0; i1<n_site; i1++) {
          for(unsigned int i2=0; i2<n_site; i2++) {
            for(unsigned int o=0; o<n_orbital; o++) {
              if (std::abs(G_omega(w,i1, i2, o)-G_old_omega(w, i1, i2, o)) > max_diff)
                max_diff = std::abs(G_omega(w,i1,i2,o)-G_old_omega(w,i1,i2,o));
            }
          }
        }
      }
      std::cout<<"convergence loop: max diff in dressed Green (Matsubara freq): "
      <<max_diff<<"\tconverged: "<<converged<<std::endl;
    }
    print_all_green_functions(basename, iteration_ctr, G0_omega, G_omega, G0_tau, G_tau, beta);
  }while ((max_diff > converged || iteration_ctr <= 1) && iteration_ctr < max_it);           
  // write G0 (to be read as an input Green function)
  G0_omega.write(parms.value_or_default("G0OMEGA_output", "G0omega_output").c_str());
}



//construct the Hilbert transformer and run a single site loop
void perform_single_site_loop(const alps::Parameters &parms){
  
  FrequencySpaceHilbertTransformer *transform_ptr;
  if(parms.defined("TWODBS"))
    transform_ptr=setup_2D_hilbert_transform(parms);
  else if (parms.defined("DOSFILE"))
    transform_ptr=setup_DOS_hilbert_transform(parms);
  else
    transform_ptr=setup_Bethe_hilbert_transform(parms);
  
  boost::shared_ptr<MatsubaraImpuritySolver> solver_ptr;
  solver_ptr.reset(new ExternalSolver(boost::filesystem::path(parms["SOLVER"])));
  selfconsistency_loop_single_site(parms, *solver_ptr, *transform_ptr);
}

