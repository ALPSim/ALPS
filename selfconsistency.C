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
  double h = static_cast<double>(parms["H"]);
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
    print_tau_green_functions(iteration_ctr++, G0_tau, G_tau, beta);
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
  matsubara_green_function_t G_omega(G_tau.ntime()-1, G_tau.nsite(), G_tau.nflavor(), diagonal);
  itime_green_function_t G_tau_old(G_tau.ntime(), G_tau.nsite(), G_tau.nflavor(), diagonal);
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

/*void adjust_chemical_potential(const alps::Parameters &parms, double current_filling, double &mu){
 
 //static double filling_lower=-std::numeric_limits<double>::infinity(), filling_upper=std::numeric_limits<double>::infinity();
 static double mu_lower=parms.value_or_default("MU_LOWER_BOUND",-std::numeric_limits<double>::infinity());
 static double mu_upper=parms.value_or_default("MU_UPPER_BOUND",std::numeric_limits<double>::infinity());
 static double target_filling=static_cast<double>(parms["FILLING"]);
 static double filling_tolerance=parms.value_or_default("FILLING_TOLERANCE", 0.002);
 //static double filling_step=parms.value_or_default("FILLING_STEP", static_cast<double>(parms["t"]));
 static double mu_old=std::numeric_limits<double>::infinity();
 static double filling_old=std::numeric_limits<double>::infinity();
 if(parms.defined("OLD_MU")){
 mu_old=parms["OLD_MU"];
 filling_old=parms["OLD_FILLING"];
 }
 static bool keep_adjusting=parms.value_or_default("KEEP_ADJUSTING_MU",true);
 double mu = static_cast<double>(parms["MU"]);
 double mu_target;
 
 if(keep_adjusting==false){
 std::cout<<"chem. potential fixed at "<<mu<<" current filling is: "<<current_filling<<" no longer adjusting. "<<std::endl;
 std::cout<<"using new mu: "<<mu_old<<std::endl;
 return mu;
 }
 std::cout<<"current filling is: "<<current_filling<<std::endl;
 std::cout<<"current mu is: "<<mu<<std::endl;
 std::cout<<"old filling is: "<<filling_old<<" old mu: "<<mu_old<<std::endl;
 std::cout<<"lower: "<<mu_lower<<std::endl;
 std::cout<<"upper: "<<mu_upper<<std::endl;
 
 if (std::abs(current_filling-target_filling) > filling_tolerance){
 std::cout<<"need to adjust filling. difference is "<<std::abs(current_filling-target_filling)<<std::endl;
 if(filling_old==std::numeric_limits<double>::infinity()){ //first time we run this
 filling_old=current_filling;
 if(current_filling>target_filling){
 mu_target=(mu==mu_upper?mu_upper+0.001:mu_upper);
 }else{
 mu_target=(mu==mu_lower?mu_lower+0.001:mu_lower);
 }
 mu_old=mu;
 }else{ //we have a valid mu_old, mu, filling_old and current filling
 mu_target=(mu_old-mu)*(target_filling-current_filling)/(filling_old-current_filling)+mu;
 if(mu_target<mu_lower){
 mu_target=(mu==mu_lower?mu_lower+0.001:mu_lower);
 }else if(mu_target>mu_upper){
 mu_target=(mu==mu_upper?mu_upper+0.001:mu_upper);
 }
 filling_old=current_filling;
 mu_old=mu;
 }
 mu=mu_target;
 std::cout<<"using new mu: "<<mu<<std::endl;
 std::stringstream mustring; mustring<<mu;
 (const_cast<alps::Parameters*>(&parms))->operator[]("MU")=mustring.str();
 }else
 keep_adjusting=false;
 return mu;
 }*/

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


