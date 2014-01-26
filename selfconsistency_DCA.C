/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005-2012 by Emanuel Gull <gullc@comp-phys.org>
 *                            Philipp Werner <werner@comp-phys.org>
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


#include "selfconsistency.h"
#include "green_function.h"
#include "fouriertransform.h"
#include "types.h"
#include <sys/types.h> //will this fail on Windows?
#include <boost/tuple/tuple.hpp> //will this fail on Windows?



void selfconsistency_loop_DCA(const alps::Parameters& parms, MatsubaraImpuritySolver& solver,
                              DCATransformer& clustertrans)
{
  unsigned int n_tau=boost::lexical_cast<unsigned int>(parms["N"]);
  unsigned int n_matsubara=boost::lexical_cast<unsigned int>(parms["NMATSUBARA"]);
  unsigned int n_orbital=parms.value_or_default("FLAVORS", 2);
  unsigned int n_site=parms.value_or_default("SITES", 1);
  
  double beta = static_cast<double>(parms["BETA"]);
  double h = static_cast<double>(parms["H"]);
  double h_init = parms.value_or_default("H_INIT", 0.);
  double converged = static_cast<double>(parms["CONVERGED"]);
  bool paramagnet = parms.value_or_default("PARAMAGNET", false);
  double relax_rate=static_cast<double>(parms.value_or_default("RELAX_RATE", 1.));
  int max_it=static_cast<int>(parms.value_or_default("MAX_IT", 1000));
  std::string basename=parms["BASENAME"];
  shape_t gf_shape=paramagnet?diagonal:blockdiagonal;
  
  //calculate initial G0
  matsubara_green_function_t G0_omega=clustertrans.initial_G0(parms, gf_shape);
  matsubara_green_function_t G_omega=G0_omega;
  
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
  double max_diff;
  double mu;
  do {
    iteration_ctr++;
    std::cout<<"starting iteration nr. "<<iteration_ctr<<std::endl;
    (*const_cast<alps::Parameters*>(&parms))["CURRENT_IT"] = iteration_ctr;
    mu = static_cast<double>(parms["MU"]);
    G_old_omega = G_omega;
    G0_old_omega = G0_omega;
    G_old_tau = G_tau;
    G0_old_tau = G0_tau;
    fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega);
    //
    //std::ofstream g0_before_symm("g0_before_symm.dat"); g0_before_symm<<std::make_pair(G0_omega, beta)<<std::endl;
    //G0_omega = clustertrans.symmetrize(G0_omega, paramagnet);
    //std::ofstream g0_after_symm("g0_after_symm.dat"); g0_after_symm<<std::make_pair(G0_omega, beta)<<std::endl;
    
    //std::cout<<"running solver in frequency space."<<std::endl<<std::flush;
    boost::tie(G_omega, G_tau) = solver.solve_omega(G0_omega, parms);
    //std::cout<<"applying spin symmetries"<<std::endl<<std::flush;
    //std::ofstream g_before_symm("g_before_symm.dat"); g_before_symm<<std::make_pair(G_omega, beta)<<std::endl;
    //G_omega = clustertrans.symmetrize(G_omega, paramagnet);
    //std::ofstream g_after_symm("g_after_symm.dat"); g_after_symm<<std::make_pair(G_omega, beta)<<std::endl;
    //exit(0);
    /*if(paramagnet){
      std::cout<<"applying point group symmetries"<<std::endl<<std::flush;
      G_tau = clustertrans.symmetrize(G_tau, paramagnet);
      G_omega = clustertrans.symmetrize_k_space(G_omega);
      G_tau = clustertrans.symmetrize_k_space(G_tau);
    }*/
    std::cout<<"running cluster transform"<<std::endl<<std::flush;
    G0_omega = clustertrans(G_omega, G0_omega, mu, h, beta);
    G0_omega.write(parms.value_or_default("G0OMEGA_output", "G0omega_output").c_str());
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
    std::cout<<"writing Green's functions."<<std::endl;
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
  } while ((max_diff > converged || iteration_ctr <= 1) && iteration_ctr < max_it);
  // write G0
  G0_omega.write(parms.value_or_default("G0OMEGA_output", "G0omega_output").c_str());
}

void selfconsistency_loop_DCA(const alps::Parameters& parms, NambuImpuritySolver& solver,
                              DCATransformer& clustertrans)
{
  unsigned int n_tau=boost::lexical_cast<unsigned int>(parms["N"]);
  unsigned int n_matsubara=boost::lexical_cast<unsigned int>(parms["NMATSUBARA"]);
  unsigned int n_site=parms.value_or_default("SITES", 1);
  
  double beta = static_cast<double>(parms["BETA"]);
  double h    = static_cast<double>(parms["H"]);
  double eta1 = static_cast<double>(parms["ETA1"]);
  double eta2 = static_cast<double>(parms["ETA2"]);
  double converged = static_cast<double>(parms["CONVERGED"]);
  bool paramagnet = parms.value_or_default("PARAMAGNET", true);
  if(!paramagnet) throw std::logic_error("please implement mixed AFM & Nambu!");
  int max_it=static_cast<int>(parms.value_or_default("MAX_IT", 1000));
  std::string basename=parms["BASENAME"];
  
  
  nambu_matsubara_green_function_t G0_omega=clustertrans.initial_G0_nambu(parms);
  nambu_matsubara_green_function_t G_omega=G0_omega;
  
  int iteration_ctr = 0;
  nambu_matsubara_green_function_t G_old_omega(G0_omega);
  nambu_matsubara_green_function_t G0_old_omega(G0_omega);
  nambu_matsubara_green_function_t selfenergy(G0_omega);
  nambu_itime_green_function_t G_tau(n_tau +1, n_site, 2, diagonal);
  nambu_itime_green_function_t G0_tau(n_tau+1, n_site, 2, diagonal);
  nambu_itime_green_function_t G_old_tau(n_tau+1, n_site, 2, diagonal);
  nambu_itime_green_function_t G0_old_tau(n_tau+1, n_site, 2, diagonal);
  //initialize space for multiple vectors
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  FourierTransformer::generate_transformer(parms, fourier_ptr);
  double max_diff;
  do {
    iteration_ctr++;
    std::cout<<"starting iteration nr. "<<iteration_ctr<<std::endl;
    double mu = static_cast<double>(parms["MU"]);
    G_old_omega = G_omega;
    G0_old_omega = G0_omega;
    G_old_tau = G_tau;
    G0_old_tau = G0_tau;
    fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega);
    std::cout<<"running NAMBU solver in omega."<<std::endl<<std::flush;
    //running additional symmetrizations:
    G0_omega = clustertrans.symmetrize_k_space(G0_omega);
    boost::tie(G_omega, G_tau) = solver.solve_nambu(G0_omega, parms);
    //std::ofstream G_tau_before_symmetrize("G_tau_before_symmetrize.dat"); G_tau_before_symmetrize<<std::make_pair(G_tau,beta);
    G_tau = clustertrans.symmetrize(G_tau, paramagnet);
    //std::ofstream G_tau_before_symmetrize_kspace("G_tau_before_symmetrize_kspace.dat"); G_tau_before_symmetrize_kspace<<std::make_pair(G_tau,beta);
    G_tau = clustertrans.symmetrize_k_space(G_tau); //<- this would symmetrize away the d-wave behavior!
    //std::ofstream G_tau_after_symmetrize_kspace("G_tau_after_symmetrize_kspace.dat"); G_tau_after_symmetrize_kspace<<std::make_pair(G_tau,beta);
    //std::ofstream G_omega_before_symmetrize("G_omega_before_symmetrize.dat"); G_omega_before_symmetrize<<std::make_pair(G_omega,beta);
    G_omega = clustertrans.symmetrize(G_omega, paramagnet);
    //std::ofstream G_omega_before_symmetrize_kspace("G_omega_before_symmetrize_kspace.dat"); G_omega_before_symmetrize_kspace<<std::make_pair(G_omega,beta);
    G_omega = clustertrans.symmetrize_k_space(G_omega);
    //std::ofstream G_omega_after_symmetrize_kspace("G_omega_after_symmetrize_kspace.dat"); G_omega_after_symmetrize_kspace<<std::make_pair(G_omega,beta);
    
    std::cout<<"running cluster transform"<<std::endl<<std::flush;
    G0_omega = clustertrans(G_omega, G0_omega, selfenergy, mu, h, eta1, eta2, beta, false);
    
    
    G0_omega.write(parms.value_or_default("G0OMEGA_output", "G0omega_output").c_str());
    std::cout<<"comparing old and new result."<<std::endl<<std::flush;
    max_diff=0;
    if(iteration_ctr>1){
      for(int f1=0;f1<2;++f1){
        for(int f2=0;f2<2;++f2){
          //comparison for the dressed Green's function in Matsubara freq.
          for(unsigned int w=0; w<n_matsubara; w++) {
            for(unsigned int i1=0; i1<n_site; i1++) {
              for(unsigned int i2=0; i2<n_site; i2++) {
                if (std::abs(G_omega(w,i1, i2, f1,f2)-G_old_omega(w, i1, i2, f1,f2)) > max_diff)
                  max_diff = std::abs(G_omega(w,i1,i2,f1,f2)-G_old_omega(w,i1,i2,f1,f2));
              }
            }
          }
        }
      }
      std::cout<<"convergence loop: max diff in dressed Green (Matsubara freq): "
      <<max_diff<<"\tconverged: "<<converged<<std::endl;
    }
    nambu_matsubara_green_function_t real_G0_omega = clustertrans.transform_into_real_space(G0_omega);
    nambu_matsubara_green_function_t real_G_omega = clustertrans.transform_into_real_space(G_omega);
    nambu_itime_green_function_t real_G0_tau = clustertrans.transform_into_real_space(G0_tau);
    nambu_itime_green_function_t real_G_tau = clustertrans.transform_into_real_space(G_tau);
    print_all_nambu_green_functions(basename, iteration_ctr, G0_omega, G_omega, G0_tau, G_tau, beta, diagonal);
    if(!parms.defined("WRITE_SPARSE"))
      print_all_nambu_green_functions(basename, iteration_ctr, real_G0_omega, real_G_omega, real_G0_tau, real_G_tau, beta, nondiagonal, "_realspace");
  } while ((max_diff > converged || iteration_ctr <= 1) && iteration_ctr < max_it);
  // write G0
  G0_omega.write(parms.value_or_default("G0OMEGA_output", "G0omega_output").c_str());
}

