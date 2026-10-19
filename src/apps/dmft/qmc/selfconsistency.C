 /*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *               2012 - 2013 by Jakub Imriska <jimriska@phys.ethz.ch>
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: selfconsistency.C 379 2009-10-07 13:58:36Z haase $ */

/// @file selfconsistency.C
/// @brief implements the selfconsistency loop functions

#include "selfconsistency.h"
#include "green_function.h"
#include "fouriertransform.h"
#include "types.h"
#include <alps/input/dmft.hpp>
#include <sys/types.h> 
#include <boost/tuple/tuple.hpp>
                                 /// @brief Run the self consistency loop for G and G0 mainly in imaginary time. Perform Fourier transformations if needed.
                                 ///
                                 /// @param parms contains the ALPS parameters needed for the simulation.
                                 /// @param solver is the impurity solver (e.g. Hirsch Fye) that creates G out of G0.
                                 /// @param hilbert is the HilbertTransformer that solves the Dyson equation, i.e. generates G0 out of G
                                 /// @param G0 is the bare Green's function in imaginary time. It has to be provided as an initial guess
                                 /// @param G is the Green's function, it does not have to be initialized but reasonable values will be returned upon completion of the loop


void selfconsistency_loop(alps::Parameters& parms, ImpuritySolver& solver, HilbertTransformer& hilbert)
{
  const alps::input::dmft_loop_parameters input =
    alps::input::read_dmft_loop_parameters(parms, alps::input::dmft_loop::imaginary_time);
  int N = input.n_tau;
  int flavors = input.sizes.flavors;
  double beta = input.beta;
  double h = input.h;
  double converged = input.converged;
  bool symmetrization = input.symmetrization;
  int max_it=input.max_it;
  std::string basename=input.basename;
  if (input.h_init) parms["H"]=*input.h_init;
  itime_green_function_t G0_tau = hilbert.initial_G0(parms);
  itime_green_function_t G_tau = G0_tau;
  itime_green_function_t G0_tau_old(G0_tau);
  int iteration_ctr=0;
  double max_diff;	
  do {
    ++iteration_ctr;
    std::cout<<"starting iteration nr. "<<iteration_ctr<<std::endl;
    double mu = input.mu;
    G0_tau_old = G0_tau;
    std::cout<<"running solver"<<std::endl;
    G_tau = solver.solve(G0_tau, parms);
    //std::cout<<"G after solver: "<<G_tau.to_multiple_vector()<<std::endl;
    std::cout<<"running Hilbert transform"<<std::endl;
    G_tau = hilbert.symmetrize(G_tau, symmetrization);
    G0_tau= hilbert(G_tau, mu, h, beta);
    parms["H"]=h;
    std::cout<<"comparing old and new results"<<std::endl;
    max_diff=0;
    for(int f=0; f<flavors;++f){
      for(int i=0; i<N; i++) {
        if (fabs(G0_tau(i,f)-G0_tau_old(i,f)) > max_diff)
          max_diff = std::abs(G0_tau(i,f)-G0_tau_old(i,f));
      }
    }	
    std::cout<<"maximum difference in G0_tau is: "<<max_diff<<std::endl;
    print_tau_green_functions(basename, iteration_ctr, G0_tau_old, G_tau, beta);
  } while (max_diff > converged  && iteration_ctr < max_it);
  std::cout<<(max_diff > converged ? "NOT " : "")<<"converged!"<<std::endl;
  // write G0 (to be read in as an input for a new simulation)
  G0_tau.write(input.g0_tau_output.c_str());
}



void F_selfconsistency_loop(alps::Parameters& parms, ImpuritySolver& solver, HilbertTransformer& hilbert)
{
  const alps::input::dmft_loop_parameters input =
    alps::input::read_dmft_loop_parameters(parms, alps::input::dmft_loop::dressed_imaginary_time);
  int N = input.n_tau;
  int flavors = input.sizes.flavors;
  double beta = input.beta;
  double h = input.h;
  double converged = input.converged;
  bool symmetrization = input.symmetrization;
  int max_it=input.max_it;
  std::string basename=input.basename;
  if (input.h_init) parms["H"]=*input.h_init;
  itime_green_function_t G_tau = hilbert.initial_G0(parms);
  itime_green_function_t G_tau_old(G_tau.ntime(), G_tau.nsite(), G_tau.nflavor());
  int iteration_ctr=0;
  double max_diff;
  do {
    ++iteration_ctr;
    std::cout<<"starting iteration nr. "<<iteration_ctr<<std::endl;
    G_tau_old=G_tau;
    
    std::cout<<"running solver"<<std::endl;
    G_tau= solver.solve(G_tau, parms); //the Werner solver WANTS a G_tau as an input. It then makes an F function out of it.
    G_tau = hilbert.symmetrize(G_tau, symmetrization);

    std::cout<<"comparing old and new results"<<std::endl;
    max_diff=0;
    for(int f=0; f<flavors;++f){
      for(int i=0; i<N; i++) {
        if (fabs(G_tau(i,f)-G_tau_old(i,f)) > max_diff)
          max_diff = fabs(G_tau(i,f)-G_tau_old(i,f));
      }
    }
    std::cout<<"maximum difference in G_tau is: "<<max_diff<<std::endl;
    print_dressed_tau_green_functions(basename, iteration_ctr, G_tau, beta);
    parms["H"]=h;
  } while (max_diff > converged && iteration_ctr < max_it);
  std::cout<<(max_diff > converged ? "NOT " : "")<<"converged!"<<std::endl;
  // write G (to be read in as an input for a new simulation)
  G_tau.write(input.g0_tau_output.c_str());

}


/// @brief Run the self consistency loop for G and G0 mainly in Matsubara frequency space. Perform Fourier transformations where needed.
///
/// @param parms contains the ALPS parameters needed for the simulation.
/// @param solver is the impurity solver (e.g. Hirsch Fye) that creates G_omega out of G0_omega and G0_tau.
/// @param hilbert is the HilbertTransformer that solves the Dyson equation, i.e. generates G0_omega out of G_omega
/// @param G0_omega is the bare Green's function in Matsubara frequency. It has to be provided as an initial guess
/// @param G_omega is the dressed Green's function, it does not have to be initialized but reasonable values will be returned upon completion of the loop

void selfconsistency_loop_omega(alps::Parameters& parms, MatsubaraImpuritySolver& solver, 
                                FrequencySpaceHilbertTransformer& hilbert) 
{
  const alps::input::dmft_loop_parameters input =
    alps::input::read_dmft_loop_parameters(parms, alps::input::dmft_loop::matsubara);
  unsigned int n_tau=input.n_tau;
  unsigned int n_matsubara=input.n_matsubara;
  unsigned int n_orbital=input.sizes.flavors;
  unsigned int n_site=input.sizes.sites;
  
  double beta = input.beta;
  double h = input.h;
  double mu = input.mu;
  double converged = input.converged;
  bool symmetrization = input.symmetrization;
  double relax_rate=input.relax_rate;
  int max_it=input.max_it;
  std::string basename=input.basename;
  
  if (input.h_init) parms["H"]=*input.h_init;
  matsubara_green_function_t G0_omega = hilbert.initial_G0(parms);
  G0_omega = hilbert.symmetrize(G0_omega, symmetrization);
  
  //define multiple vectors
  matsubara_green_function_t G_omega(n_matsubara, n_site, n_orbital);
  matsubara_green_function_t G_omega_old(n_matsubara, n_site, n_orbital);
  matsubara_green_function_t G0_omega_old(n_matsubara, n_site, n_orbital);
  itime_green_function_t G_tau(n_tau +1, n_site, n_orbital);
  itime_green_function_t G0_tau(n_tau+1, n_site, n_orbital);
  itime_green_function_t G0_tau_old(n_tau+1, n_site, n_orbital);

  boost::shared_ptr<FourierTransformer> fourier_ptr;
  FourierTransformer::generate_transformer(parms, fourier_ptr);
  fourier_ptr->backward_ft(G0_tau, G0_omega);
  
  if (input.initial_g0_tau_output)           // it is not needed to store it by default, as it will be stored in the 1st iteration as G0_tau_1
    G0_tau.write(input.initial_g0_tau_output->c_str());
  
  double max_diff=0.;	
  int iteration_ctr = 0;
  do {
    iteration_ctr++;
    std::cout<<"starting iteration nr. "<<iteration_ctr<<std::endl;
    G_omega_old = G_omega;
    G0_omega_old = G0_omega;
    G0_tau_old = G0_tau;
    std::cout<<"running solver."<<std::endl<<std::flush;
    boost::tie(G_omega, G_tau) = solver.solve_omega(G0_omega,parms);
    G_tau = hilbert.symmetrize(G_tau, symmetrization);
    G_omega = hilbert.symmetrize(G_omega, symmetrization);
    std::cout<<"running Hilbert transform"<<std::endl<<std::flush;
    G0_omega = hilbert(G_omega, G0_omega, mu, h, beta);
    //relaxation to speed up/slow down convergence
    if(relax_rate !=1){
      std::cout<<"using over/underrelaxation with rate: "<<relax_rate<<std::endl;
      for(unsigned int o=0;o<n_orbital;++o){
        for(unsigned int i2=0;i2<n_site;++i2){
          for(unsigned int i1=0;i1<n_site;++i1){
            for(unsigned int w=0;w<n_matsubara;++w){
              G0_omega(w, i1, i2, o)=relax_rate*G0_omega(w, i1, i2, o) + (1.-relax_rate)*G0_omega_old(w, i1, i2, o);
            }
          }
        }
      }
    }
    if (iteration_ctr==1) {
      parms["H"]=h;
      FourierTransformer::generate_transformer(parms, fourier_ptr);
    }
    fourier_ptr->backward_ft(G0_tau, G0_omega);
    if(iteration_ctr>1){
      //comparison for the dressed Green's function in Matsubara freq.
      std::cout<<"comparing old and new result."<<std::endl;
      max_diff=0;
      for(unsigned int o=0;o<n_orbital;++o){
        for(unsigned int i2=0;i2<n_site;++i2){
          for(unsigned int i1=0;i1<n_site;++i1){
            for(unsigned int w=0;w<n_matsubara;++w){
              if (std::abs(G_omega(w,i1,i2,o)-G_omega_old(w,i1,i2,o)) > max_diff)
                max_diff = std::abs(G_omega(w,i1,i2,o)-G_omega_old(w,i1,i2,o));
            }
          }
        }
      }
      std::cout<<"convergence loop: max diff in dressed Green (Matsubara freq): "
      <<max_diff<<"\t(convergency criterion: "<<converged<<")"<<std::endl<<std::flush;
    }
    print_all_green_functions(basename, iteration_ctr, G0_omega_old, G_omega, G0_tau_old, G_tau, beta);
  }while ((max_diff > converged || iteration_ctr <= 1) && iteration_ctr < max_it);           
  // write G0 (to be read as an input Green function)
  G0_omega.write(input.g0_omega_output.c_str());
}


