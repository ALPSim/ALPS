/*****************************************************************************
 *
 * Copyright (C) 2005-2015 Emanuel Gull <egull@umich.edu>
 *
 *****************************************************************************/

#include "green_function.h"

#include "fourier_transform.h"

#include "solver.h"
#include "externalsolver.h"
#include "selfconsistency.h"
#include "hilbert_transform.h"

void selfconsistency_loop_single_site(const alps::params& parms, ImpuritySolver& solver,
                                FrequencySpaceHilbertTransformer& hilbert) 
{
  unsigned int n_tau=parms["N"];
  unsigned int n_matsubara=parms["NMATSUBARA"];
  unsigned int n_orbital=parms["FLAVORS"]|2;
  unsigned int n_site=parms["SITES"]|1;
  
  double beta = parms["BETA"];
  double h = parms["H"];
  double converged = parms["CONVERGED"];
  bool paramagnet = parms["PARAMAGNET"]|false;
  double relax_rate=parms["RELAX_RATE"]|1.;
  int max_it=parms["MAX_IT"]|1000;
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
    double mu = parms["MU"];
    G_old_omega = G_omega;
    G0_old_omega = G0_omega;
    G_old_tau = G_tau;
    G0_old_tau = G0_tau;
    fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega);
    //fourier.setmu(mu);
    std::cout<<"running solver."<<std::endl<<std::flush;
    boost::tie(G_omega, G_tau) = solver.solve_matsubara(G0_omega,parms);
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
    print_all_green_functions(basename, iteration_ctr, G0_omega, G_omega, G0_tau, G_tau, beta);
  }while ((max_diff > converged || iteration_ctr <= 1) && iteration_ctr < max_it);           
  // write G0 (to be read as an input Green function)
  G0_omega.write(std::string(parms["G0OMEGA_output"]|"G0omega_output").c_str());
}



//construct the Hilbert transformer and run a single site loop
void perform_single_site_loop(const alps::params &parms){
  
  FrequencySpaceHilbertTransformer *transform_ptr;
  if(parms.defined("TWODBS"))
    transform_ptr=setup_2D_hilbert_transform(parms);
  else if (parms.defined("DOSFILE"))
    transform_ptr=setup_DOS_hilbert_transform(parms);
  else
    transform_ptr=setup_Bethe_hilbert_transform(parms);
  
  std::string launcher_path=parms["LAUNCHER"]|"";
  std::string solver_path=parms["SOLVER"];
  ImpuritySolver* solver_ptr=new ExternalSolver(launcher_path+" "+solver_path);
  selfconsistency_loop_single_site(parms, *solver_ptr, *transform_ptr);
  delete solver_ptr;
}

