/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005-2006 by Philipp Werner <werner@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
*                            Emanuel Gull <gullc@comp-phys.org>
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
#include <sys/types.h> //will this fail on Windows?
#include <boost/tuple/tuple.hpp> //will this fail on Windows?
/// @brief Run the self consistency loop for G and G0 mainly in imaginary time. Perform Fourier transformations if needed.
///
/// @param parms contains the ALPS parameters needed for the simulation.
/// @param solver is the impurity solver (e.g. Hirsch Fye) that creates G out of G0.
/// @param hilbert is the HilbertTransformer that solves the Dyson equation, i.e. generates G0 out of G
/// @param G0 is the bare Green's function in imaginary time. It has to be provided as an initial guess
/// @param G is the Green's function, it does not have to be initialized but reasonable values will be returned upon completion of the loop
void selfconsistency_loop(alps::Parameters& parms, ImpuritySolver& solver, HilbertTransformer& hilbert, itime_multiple_vector_t& G0, itime_multiple_vector_t& G)
{
  int N = static_cast<int>(parms["N"]);
  double beta = static_cast<double>(parms["BETA"]);
  double mu = static_cast<double>(parms["MU"]);
  double h = static_cast<double>(parms["H"]);
  double converged = static_cast<double>(parms["CONVERGED"]);
  bool paramagnet = static_cast<bool>(parms["PARAMAGNET"]);
  
  multiple_vector_type G0_old(G0);
  int iteration_ctr=0;
    
  double max_diff;	
  do {
    G0_old = G0;
    std::cout<<"running solver"<<std::endl;
    G = solver.solve(G0, parms);
    std::cout<<"running Hilbert transform"<<std::endl;
    G0 = hilbert(G, mu, h, beta, paramagnet);
    std::cout<<"comparing old and new results"<<std::endl;
	max_diff=0;
    for(int i=0; i<N; i++) {
	  if (fabs(G0.first[i]-G0_old.first[i]) > max_diff)
	    max_diff = fabs(G0.first[i]-G0_old.first[i]);
	  if (fabs(G0.second[i]-G0_old.second[i]) > max_diff)
	    max_diff = fabs(G0.second[i]-G0_old.second[i]);	
	}	
  std::cout<<"maximum difference in G0 is: "<<max_diff<<std::endl;
  print_tau_green_functions(iteration_ctr++, G0, G, beta);
  } while (max_diff > converged);
  std::cout<<"converged!"<<std::endl;
}

/// @brief Run the self consistency loop for G and G0 mainly in Matsubara frequency space. Perform Fourier transformations where needed.
///
/// @param parms contains the ALPS parameters needed for the simulation.
/// @param solver is the impurity solver (e.g. Hirsch Fye) that creates G_omega out of G0_omega and G0_tau.
/// @param hilbert is the HilbertTransformer that solves the Dyson equation, i.e. generates G0_omega out of G_omega
/// @param G0_omega is the bare Green's function in Matsubara frequency. It has to be provided as an initial guess
/// @param G_omega is the dressed Green's function, it does not have to be initialized but reasonable values will be returned upon completion of the loop

void selfconsistency_loop_omega_tau(const alps::Parameters& parms, FrequencyTimeImpuritySolver& solver, FrequencySpaceHilbertTransformer& hilbert,
									freq_multiple_vector_t& G0_omega, freq_multiple_vector_t& G_omega){
	int n_tau=boost::lexical_cast<int>(parms["N"]);
	int n_matsubara=boost::lexical_cast<int>(parms["N"]);
	double beta = static_cast<double>(parms["BETA"]);
	double mu = static_cast<double>(parms["MU"]);
	double h = static_cast<double>(parms["H"]);
	double converged = static_cast<double>(parms["CONVERGED"]);
	bool paramagnet = static_cast<bool>(parms["PARAMAGNET"]);
	double t=static_cast<double>(parms["t"]);
					
	double mu_up = mu+h/2.;
	double mu_down = mu-h/2.;
	double tsq=t*t;
	int iteration_ctr = 0;
	freq_multiple_vector_t G_old_omega(G0_omega);
    	itime_multiple_vector_t G0_tau;
	itime_multiple_vector_t G_tau;
	itime_multiple_vector_t G_old_tau;
	
	double max_diff;	
	do {
		iteration_ctr++;
		std::cout<<"starting iteration nr. "<<iteration_ctr<<std::endl;
		G_old_omega = G_omega;
		G_old_tau = G_tau;
		G0_tau.first.resize(n_tau+1);
		G0_tau.second.resize(n_tau+1);

		backward_FT(beta, -mu_up,   mu_up*mu_up+tsq,     G0_tau.first,  G0_omega.first); //make G0(tau) out of G0(iomega)
		backward_FT(beta, -mu_down, mu_down*mu_down+tsq, G0_tau.second, G0_omega.second);
		
		std::cout<<"running solver."<<std::endl<<std::flush;
		boost::tie(G_omega, G_tau) = solver.solve_omega(G0_omega, G0_tau,parms);
                //if Green_itime_m is defined we know G in tau but not in
                //omega. Compute G in omega via fourier transform
                if(parms.defined("DIRECT_GREEN_ITIME")|| parms.defined("GREEN_ITIME_M")){
                  if(paramagnet){
                    complex_vector G_omega_para(n_matsubara);
                    double_vector G_para(n_tau+1);
                    for (uint i=0; i<n_tau+1; i++){
                      G_para[i]=0.5*(G_tau.first[i]+G_tau.second[i]);
                    }
                    forward_FT(beta, -mu, mu*mu+tsq, G_para, G_omega_para); //make G(iomega) out of G(tau)
                    for(uint i=0;i<n_matsubara;++i){
                      G_omega.first[i]=G_omega_para[i];
                      G_omega.second[i]=G_omega_para[i];
                    }
                  }
                  else{
                    forward_FT(beta, -mu_up  , mu_up*mu_up+tsq    , G_tau.first, G_omega.first);
                    forward_FT(beta, -mu_down, mu_down*mu_down+tsq, G_tau.second, G_omega.second);
                                                                    
                  }
                }
		std::cout<<"running Hilbert transformer."<<std::endl<<std::flush;
		G0_omega = hilbert(G_omega, mu, h, beta, paramagnet);
                std::cout<<"kinetic energy is: "<<kinetic_energy(G_tau, beta, t);
		if(!parms.defined("GREEN_ITIME_M")){
			backward_FT(beta, -mu_up,   mu_up*mu_up+tsq,     G_tau.first,  G_omega.first);
			backward_FT(beta, -mu_down, mu_down*mu_down+tsq, G_tau.second, G_omega.second);
		}
                if(paramagnet){
                  for(uint i=0;i<n_tau+1;++i){
                    G_tau.first[i]=0.5*(G_tau.first[i]+G_tau.second[i]);
                    G0_tau.first[i]=0.5*(G0_tau.first[i]+G0_tau.second[i]);
                    G_tau.second[i]=G_tau.first[i];
                    G0_tau.second[i]=G0_tau.first[i];
                  }
                }
		std::cout<<"comparing old and new result."<<std::endl<<std::flush;
		max_diff=0;
		if(iteration_ctr>1){
			//comparison for the dressed Green's function in Matsubara freq.
			for(uint i=0; i<n_matsubara; i++) {
				if (fabs((G_omega.first[i]-G_old_omega.first[i]).imag()) > max_diff)
					max_diff = fabs((G_omega.first[i]-G_old_omega.first[i]).imag());
				if (fabs((G_omega.second[i]-G_old_omega.second[i]).imag()) > max_diff)
					max_diff = fabs((G_omega.second[i]-G_old_omega.second[i]).imag());	
			}
			std::cerr<<"convergence loop: max diff in dressed Green (Matsubara freq): "<<max_diff<<"\tconverged: "<<converged<<std::endl;
		}
		//write checkpoint - clean this up in the future!
		print_all_green_functions(iteration_ctr, G0_omega, G_omega, G0_tau, G_tau, beta, paramagnet);
	}while (max_diff > converged || iteration_ctr <= 1); 
	
}
