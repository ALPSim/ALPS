/*
 * selfconsistency_DCAPlus.C
 *
 *  Created on: Feb 4, 2015
 *      Author: pakij
 */

/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005-2014 by Emanuel Gull <gullc@comp-phys.org>
 *
 * To be used with permission by the authors
 *
 *****************************************************************************/


#include "types.h"
#include "green_function.h"
#include "selfconsistency.h"
#include "fouriertransform.h"
#include "externalsolver.h"
#include <sys/types.h>
#include <boost/tuple/tuple.hpp>
#include <alps/parameter.h>

#include "cluster.h"
#include "DCA_transformer.hxx"
#include "SimpleDCA_transformer.hxx"
#include "AFMDCA_transformer.hxx"


void selfconsistency_loop_DCAPlus(const alps::params& parms, ImpuritySolver& solver,
                              DCATransformer& clustertrans)
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
  bool readwrite_real_space=parms["READWRITE_REAL_SPACE"]|false;
  bool phase=parms["PHASE"]|false;
  shape_t gf_shape=paramagnet?diagonal:blockdiagonal;

  //calculate initial G0
  matsubara_green_function_t G0_omega=clustertrans.initial_G0(parms, gf_shape);
  matsubara_green_function_t G_omega=G0_omega;

  //Print a file that lists the momenta involved in the calculation
  //They appear in the same order as the columns of the various output files
  std::string momenta_filename= "momenta";
  clustertrans.print_momenta(momenta_filename.c_str(), G_omega);

  int iteration_ctr = 0;
  matsubara_green_function_t G_old_omega(G0_omega);
  matsubara_green_function_t G0_old_omega(G0_omega);
  itime_green_function_t G_tau(n_tau +1, n_site, n_orbital, gf_shape);
  itime_green_function_t G0_tau(n_tau+1, n_site, n_orbital, gf_shape);
  itime_green_function_t G_old_tau(n_tau+1, n_site, n_orbital, gf_shape);
  itime_green_function_t G0_old_tau(n_tau+1, n_site, n_orbital, gf_shape);
  //initialize space for multiple vectors
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  FourierTransformer::generate_transformer(parms, fourier_ptr);
  double max_diff,norm2_diff;
  double mu;
  do {
    iteration_ctr++;
    std::cout<<"starting iteration nr. "<<iteration_ctr<<std::endl;
    (*const_cast<alps::params*>(&parms))["CURRENT_IT"] = iteration_ctr;
    mu = parms["MU"];
    G_old_omega = G_omega;
    G0_old_omega = G0_omega;
    G_old_tau = G_tau;
    G0_old_tau = G0_tau;
    fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega);
    //
    G0_omega = clustertrans.symmetrize(G0_omega, paramagnet);

    //std::cout<<"running solver in frequency space."<<std::endl<<std::flush;
    if(!readwrite_real_space){
      boost::tie(G_omega, G_tau) = solver.solve_matsubara(G0_omega, parms, &clustertrans);
    }else{
      matsubara_green_function_t G0_omega_rs=clustertrans.transform_into_real_space(G0_omega);
      matsubara_green_function_t G_omega_rs(n_matsubara, n_site, n_orbital, nondiagonal);
      itime_green_function_t G_tau_rs(n_tau+1, n_site, n_orbital, nondiagonal);
      boost::tie(G_omega_rs, G_tau_rs) = solver.solve_matsubara(G0_omega_rs, parms,&clustertrans);
      G_omega=clustertrans.transform_into_k_space(G_omega_rs);
      G_tau=clustertrans.transform_into_k_space(G_tau_rs);
    }
    G_omega = clustertrans.symmetrize(G_omega, paramagnet);
    if(paramagnet && !phase ) G_omega = clustertrans.symmetrize_k_space(G_omega);
    G0_omega = clustertrans(G_omega, G0_omega, mu, h, beta);
    G0_omega.write(std::string(parms["G0OMEGA_output"]|"G0omega_output").c_str());
    std::cout<<"comparing old and new result."<<std::endl<<std::flush;
    max_diff=0;
    norm2_diff=0;
    if(iteration_ctr>1){
      //comparison for the dressed Green's function in Matsubara freq.
      for(unsigned int w=0; w<n_matsubara; w++) {
        for(unsigned int i1=0; i1<n_site; i1++) {
          for(unsigned int i2=0; i2<n_site; i2++) {
            for(unsigned int o=0; o<n_orbital; o++) {
              norm2_diff += std::abs(G_omega(w,i1,i2,o)-G_old_omega(w,i1,i2,o))*std::abs(G_omega(w,i1,i2,o)-G_old_omega(w,i1,i2,o));
              if (std::abs(G_omega(w,i1, i2, o)-G_old_omega(w, i1, i2, o)) > max_diff)
                max_diff = std::abs(G_omega(w,i1,i2,o)-G_old_omega(w,i1,i2,o));
            }
          }
        }
      }
      std::cout<<"convergence loop: max diff in dressed Green (Matsubara freq): "
      <<max_diff<<"\tconverged: "<<converged<<std::endl;
      std::cout<<"norm2_diff = "<<std::sqrt(norm2_diff)<<std::endl;
    }
	std::cout << "writing Green's functions." << std::endl;

	if (parms["PRINTHEADERS"] | false) {
		G0_omega.set_header(clustertrans.generate_header_kspace(G0_omega));
		G_omega.set_header(clustertrans.generate_header_kspace(G_omega));
	}
	print_all_green_functions(basename, iteration_ctr, G0_omega, G_omega, G0_tau, G_tau, beta);
    if(!parms.defined("WRITE_SPARSE")){
      std::cout<<"writing Green's functions in real space"<<std::endl;
      matsubara_green_function_t real_G0_omega = clustertrans.transform_into_real_space(G0_omega);
      matsubara_green_function_t real_G_omega = clustertrans.transform_into_real_space(G_omega);
      itime_green_function_t real_G0_tau = clustertrans.transform_into_real_space(G0_tau);
      itime_green_function_t real_G_tau = clustertrans.transform_into_real_space(G_tau);
      print_all_green_functions(basename, iteration_ctr, real_G0_omega, real_G_omega, real_G0_tau, real_G_tau, beta, "_realspace");
    }
    std::cout<<"...continuing."<<std::endl;

/// TODO, rewriting to parms so that mu_orig=previous mu iteration, allows for damping in filled cases.
    if(parms.defined("MODIFY_PARAMETER_FILE")){
        std::string basename_str=parms["BASENAME"];
	std::ofstream os(basename_str.c_str());
       (*const_cast<alps::params*>(&parms))["MU"]=mu;

        os<<parms;
        os.close();
    }

  } while ((max_diff > converged || iteration_ctr <= 1) && iteration_ctr < max_it);
  // write G0
  std::string G0_output_string=parms["G0OMEGA_output"]|"G0omega_output";
  G0_omega.write(G0_output_string.c_str());
}

void perform_DCAPlus_loop(const alps::params &parms){
  std::cout<<"running a DCA calculation"<<std::endl;
  ClusterTransformer *transform_ptr=ClusterTransformer::generate_transformer(parms);
  std::string solver_path=parms["SOLVER"];
  ImpuritySolver* solver_ptr=new ExternalSolver(boost::filesystem::path(solver_path));
  if(!(parms["NAMBU"]|false)){
    selfconsistency_loop_DCA(parms, *solver_ptr, *(DCATransformer*)transform_ptr);
  }else{
	    std::cout<<"In perform_DCAPLus_loop:"<<std::endl;
	    boost::throw_exception(std::logic_error("Nambu is not implemented for DCAPlus"));
  }
  delete transform_ptr;
  delete solver_ptr;
}





