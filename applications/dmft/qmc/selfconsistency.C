 /*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
 * THIS SOFTWARE NEEDS AN APPROPRIATE LICENSE BLOCK HERE
 *****************************************************************************/

/* $Id: selfconsistency.C 379 2009-10-07 13:58:36Z haase $ */

/// @file selfconsistency.C
/// @brief implements the selfconsistency loop functions

#include "selfconsistency.h"
#include "green_function.h"
#include "fouriertransform.h"
#include "types.h"
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
  int N = static_cast<int>(parms["N"]);
  int flavors = parms.value_or_default("FLAVORS", 2);
  double beta = static_cast<double>(parms["BETA"]);
  double h_old = static_cast<double>(parms["H"]);
  double h_init = parms.value_or_default("H_INIT", 0.);
  double h = h_old + h_init;
  (*const_cast<alps::Parameters*>(&parms))["H"] = h;
  double converged = static_cast<double>(parms["CONVERGED"]);
  bool paramagnet = parms.value_or_default("PARAMAGNET", false);
  //bool degenerate = parms.value_or_default("DEGENERATE", false);
  int max_it=static_cast<int>(parms.value_or_default("MAX_IT", 1000));

  itime_green_function_t G0_tau = hilbert.initial_G0(parms);
  itime_green_function_t G_tau = G0_tau;
  itime_green_function_t G0_tau_old(G0_tau);

  int iteration_ctr=0;
  double max_diff;	
  do {
    double mu = static_cast<double>(parms["MU"]);
    G0_tau_old = G0_tau;
    std::cout<<"running solver"<<std::endl;
    G_tau = solver.solve(G0_tau, parms);
    //std::cout<<"G after solver: "<<G_tau.to_multiple_vector()<<std::endl;
    std::cout<<"running Hilbert transform"<<std::endl;
    G_tau = hilbert.symmetrize(G_tau, paramagnet);
    G0_tau= hilbert(G_tau, mu, h, beta);
    std::cout<<"comparing old and new results"<<std::endl;
    max_diff=0;
    for(int f=0; f<flavors;++f){
      for(int i=0; i<N; i++) {
        if (fabs(G0_tau(i,f)-G0_tau_old(i,f)) > max_diff)
          max_diff = fabs(G0_tau(i,f)-G0_tau_old(i,f));
      }
    }	
    std::cout<<"maximum difference in G0 is: "<<max_diff<<std::endl;
    if (iteration_ctr == 0) 
      (*const_cast<alps::Parameters*>(&parms))["H"] = h_old;
    print_tau_green_functions(iteration_ctr++, G0_tau.to_multiple_vector(), G_tau.to_multiple_vector(), beta);
  } while (max_diff > converged  && iteration_ctr < max_it);
  std::cout<<"converged!"<<std::endl;
  // write G0 (to be read in as an input for a new simulation)
  G0_tau.write(parms.value_or_default("G0OMEGA_output", "G0omega_output").c_str());
}



void F_selfconsistency_loop(alps::Parameters& parms, ImpuritySolver& solver,  itime_green_function_t& G_tau)
{
  int N = static_cast<int>(parms["N"]);
  int flavors = parms.value_or_default("FLAVORS", 2);;
  double beta = static_cast<double>(parms["BETA"]);
  //double mu = static_cast<double>(parms["MU"]);
  //double h = parms.value_or_default("H", 0);;
  double converged = static_cast<double>(parms["CONVERGED"]);
  bool paramagnet = parms.value_or_default("PARAMAGNET", false);
  //bool degenerate = parms.value_or_default("DEGENERATE", false);
  matsubara_green_function_t G_omega(G_tau.ntime()-1, G_tau.nsite(), G_tau.nflavor());
  itime_green_function_t G_tau_old(G_tau.ntime(), G_tau.nsite(), G_tau.nflavor());
  int iteration_ctr=0;
  double max_diff;
  do {
    std::cout<<"running solver"<<std::endl;
    G_tau_old=G_tau;
    G_tau= solver.solve(G_tau, parms); //the Werner solver WANTS a G_tau as an input. It then makes an F function out of it.
                                       //symmetrize
    if(paramagnet){
      for(unsigned int f=0;f<G_tau.nflavor();f+=2){
        for(unsigned int i=0;i<G_tau.ntime();++i){
          G_tau(i,f  )=0.5*(G_tau(i,f)+G_tau(i,f+1));
          G_tau(i,f+1)=G_tau(i,f);
        }
      }
    }
    std::cout<<"comparing old and new results"<<std::endl;
    max_diff=0;
    for(int f=0; f<flavors;++f){
      for(int i=0; i<N; i++) {
        if (fabs(G_tau(i,f)-G_tau_old(i,f)) > max_diff)
          max_diff = fabs(G_tau(i,f)-G_tau_old(i,f));
      }
    }
    std::cout<<"maximum difference in G is: "<<max_diff<<std::endl;
    print_dressed_tau_green_functions(iteration_ctr++, G_tau, beta);
  } while (max_diff > converged);
  std::cout<<"converged!"<<std::endl;
  // write G (to be read in as an input for a new simulation)
  G_tau.write(parms.value_or_default("G0OMEGA_output", "G0omega_output").c_str());

}


/// @brief Run the self consistency loop for G and G0 mainly in Matsubara frequency space. Perform Fourier transformations where needed.
///
/// @param parms contains the ALPS parameters needed for the simulation.
/// @param solver is the impurity solver (e.g. Hirsch Fye) that creates G_omega out of G0_omega and G0_tau.
/// @param hilbert is the HilbertTransformer that solves the Dyson equation, i.e. generates G0_omega out of G_omega
/// @param G0_omega is the bare Green's function in Matsubara frequency. It has to be provided as an initial guess
/// @param G_omega is the dressed Green's function, it does not have to be initialized but reasonable values will be returned upon completion of the loop

void selfconsistency_loop_omega(const alps::Parameters& parms, MatsubaraImpuritySolver& solver, 
                                FrequencySpaceHilbertTransformer& hilbert) 
{
  unsigned int n_tau=boost::lexical_cast<unsigned int>(parms["N"]);
  unsigned int n_matsubara=boost::lexical_cast<unsigned int>(parms["NMATSUBARA"]);
  unsigned int n_orbital=parms.value_or_default("FLAVORS", 2);
  unsigned int n_site=parms.value_or_default("SITES", 1);
  
  double beta = static_cast<double>(parms["BETA"]);
  double h_old = static_cast<double>(parms["H"]);
  double h_init = parms.value_or_default("H_INIT", 0.);
  double h = h_old + h_init;
  (*const_cast<alps::Parameters*>(&parms))["H"] = h;
  double converged = static_cast<double>(parms["CONVERGED"]);
  bool paramagnet = parms.value_or_default("PARAMAGNET", false);
  //bool degenerate = parms.value_or_default("DEGENERATE", false);
  double relax_rate=static_cast<double>(parms.value_or_default("RELAX_RATE", 1.));
  int max_it=static_cast<int>(parms.value_or_default("MAX_IT", 1000));
  
  matsubara_green_function_t G0_omega = hilbert.initial_G0(parms);
  matsubara_green_function_t G_omega = G0_omega;

  int iteration_ctr = 0;
  //define multiple vectors
  matsubara_green_function_t G_old_omega(G0_omega);
  matsubara_green_function_t G_omega_shifted(G0_omega);
  matsubara_green_function_t G0_old_omega(G0_omega);
  itime_green_function_t G_tau(n_tau +1, n_site, n_orbital);
  itime_green_function_t G0_tau(n_tau+1, n_site, n_orbital);
  itime_green_function_t G_old_tau(n_tau+1, n_site, n_orbital);
  itime_green_function_t G0_old_tau(n_tau+1, n_site, n_orbital);
  //initialize space for multiple vectors
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  double max_diff;	
  do {
    iteration_ctr++;
    std::cout<<"starting iteration nr. "<<iteration_ctr<<std::endl;
    double mu = static_cast<double>(parms["MU"]);
    G_old_omega = G_omega;
    G0_old_omega = G0_omega;
    G_old_tau = G_tau;
    G0_old_tau = G0_tau;
    FourierTransformer::generate_transformer(parms, fourier_ptr);
    fourier_ptr->backward_ft(G0_tau, G0_omega);
    //fourier.setmu(mu);
    std::cout<<"running solver."<<std::endl<<std::flush;
    boost::tie(G_omega, G_tau) = solver.solve_omega(G0_omega,parms);
    G_tau = hilbert.symmetrize(G_tau, paramagnet);
    G_omega = hilbert.symmetrize(G_omega, paramagnet);
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
    print_all_green_functions(iteration_ctr, G0_omega, G_omega, G0_tau, G_tau, beta);
    if (iteration_ctr == 1) 
      (*const_cast<alps::Parameters*>(&parms))["H"] = h_old;
  }while ((max_diff > converged || iteration_ctr <= 1) && iteration_ctr < max_it);           
  // write G0 (to be read as an input Green function)
  G0_omega.write(parms.value_or_default("G0OMEGA_output", "G0omega_output").c_str());
}


