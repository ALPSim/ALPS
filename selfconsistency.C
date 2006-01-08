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
#include <sys/types.h> //will this fail on Windows?
#include <boost/tuple/tuple.hpp> //will this fail on Windows?
#include <unistd.h>
#include <fstream>
#include <sstream>
/// @brief Run the self consistency loop for G and G0 mainly in imaginary time. Perform Fourier transformations if needed.
///
/// @param parms contains the ALPS parameters needed for the simulation.
/// @param solver is the impurity solver (e.g. Hirsch Fye) that creates G out of G0.
/// @param hilbert is the HilbertTransformer that solves the Dyson equation, i.e. generates G0 out of G
/// @param G0 is the bare Green's function in imaginary time. It has to be provided as an initial guess
/// @param G is the Green's function, it does not have to be initialized but reasonable values will be returned upon completion of the loop
void selfconsistency_loop(alps::Parameters& parms, ImpuritySolver& solver, HilbertTransformer& hilbert, itime_multiple_vector_t& G0, itime_multiple_vector_t& G)
{
  
  double beta = static_cast<double>(parms["BETA"]);
  double mu = static_cast<double>(parms["MU"]);
  double h = static_cast<double>(parms["H"]);
  double converged = static_cast<double>(parms["CONVERGED"]);
  bool paramagnet = static_cast<bool>(parms["PARAMAGNET"]);
  
  multiple_vector_type G_old(G0);
  G = G0;	  
  double max_diff;
  	
  do {
    G_old = G;
    std::cout<<"running solver"<<std::endl;
    G = solver.solve(G0, parms);
    std::cout<<"running Hilbert transform"<<std::endl;
    G0 = hilbert(G, mu, h, beta, paramagnet);
    std::cout<<"comparing old and new results"<<std::endl;
	max_diff=0;
    for(int i=0; i<G.first.size(); i++) {
	  if (fabs(G.first[i]-G_old.first[i]) > max_diff)
	    max_diff = fabs(G.first[i]-G_old.first[i]);
	  if (fabs(G.second[i]-G_old.second[i]) > max_diff)
	    max_diff = fabs(G.second[i]-G_old.second[i]);	
	}	
    std::cout<<"maximum difference in G is: "<<max_diff<<std::endl;
  } while (max_diff > converged);
  
}

/// @brief Run the self consistency loop for G and F in imaginary time. Perform Fourier transform if needed.
///
/// @param parms contains the ALPS parameters needed for the simulation.
/// @param solver is the impurity solver that creates G out of F.
/// @param hilbert is the HilbertTransformer that solves the Dyson equation, i.e. generates F out of G
/// @param F is the F-function in imaginary time. The initial guess for F is computed from the initial guess for G0
/// @param G is the Green's function, it does not have to be initialized but reasonable values will be returned upon completion of the loop
void F_selfconsistency_loop(alps::Parameters& parms, ImpuritySolver& solver, HilbertTransformer& hilbert, itime_multiple_vector_t& G0, itime_multiple_vector_t& G)
{
  
  double beta = static_cast<double>(parms["BETA"]);
  double mu = static_cast<double>(parms["MU"]);
  double h = static_cast<double>(parms["H"]);
  double converged = static_cast<double>(parms["CONVERGED"]);
  bool paramagnet = static_cast<bool>(parms["PARAMAGNET"]);
  
  itime_multiple_vector_t F = hilbert(G0, mu, h, beta, paramagnet);
  itime_multiple_vector_t G_old(G0);
  G = G0;
  double max_diff;	

  do {
    G_old = G;
    std::cout<<"running solver"<<std::endl;
    G = solver.solve(F, parms);
    std::cout<<"running Hilbert transform"<<std::endl;
    F = hilbert(G, mu, h, beta, paramagnet); // mu, h, beta not needed for semi-circle, but for other densities of states
    std::cout<<"comparing old and new results"<<std::endl;
	max_diff=0;
    for(int i=0; i<G.first.size(); i++) {
	  if (fabs(G.first[i]-G_old.first[i]) > max_diff)
	    max_diff = fabs(G.first[i]-G_old.first[i]);
	  if (fabs(G.second[i]-G_old.second[i]) > max_diff)
	    max_diff = fabs(G.second[i]-G_old.second[i]);	
	}
    std::cout<<"maximum difference in G is: "<<max_diff<<std::endl;
  } while (max_diff > converged);
  
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
	int N=boost::lexical_cast<int>(parms["N"]);
	int NMATSUBARA = static_cast<int>(parms["NMATSUBARA"]);
	double beta = static_cast<double>(parms["BETA"]);
	double mu = static_cast<double>(parms["MU"]);
	double h = static_cast<double>(parms["H"]);
	double converged = static_cast<double>(parms["CONVERGED"]);
	int paramagnet = static_cast<int>(parms["PARAMAGNET"]);
	double t=static_cast<double>(parms["t"]);
					
	double mu_up = mu+h/2.;
	double mu_down = mu-h/2.;
	double tsq=t*t;
	int iteration_ctr = 0;
	pid_t process_id=getpid();
		
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
		G0_tau.first.resize(N);
		G0_tau.second.resize(N);

		backward_FT(beta, -mu_up,   mu_up*mu_up+tsq,     G0_tau.first,  G0_omega.first);
		backward_FT(beta, -mu_down, mu_down*mu_down+tsq, G0_tau.second, G0_omega.second);
		
		std::cout<<"running solver."<<std::endl<<std::flush;
		boost::tie(G_omega, G_tau) = solver.solve_omega(G0_omega, G0_tau,parms);
		std::cout<<"running Hilbert transformer."<<std::endl<<std::flush;
		G0_omega = hilbert(G_omega, mu, h, beta, paramagnet);
		std::cout<<"comparing old and new result."<<std::endl<<std::flush;
		max_diff=0;
		
		if(!parms.defined("DIRECT_GREEN_ITIME")){
			backward_FT(beta, -mu_up,   mu_up*mu_up+tsq,     G_tau.first,  G_omega.first);
			backward_FT(beta, -mu_down, mu_down*mu_down+tsq, G_tau.second, G_omega.second);
		}
		if(iteration_ctr>1){
			//comparison for the dressed Green's function in imaginary time.
		/*	for(int i=0; i<N; i++) {
				if (fabs((G_tau.first[i]-G_old_tau.first[i])) > max_diff)
					max_diff = fabs((G_tau.first[i]-G_old_tau.first[i]));
				if (fabs((G_tau.second[i]-G_old_tau.second[i])) > max_diff)
					max_diff = fabs((G_tau.second[i]-G_old_tau.second[i]));	
			}
			std::cerr<<"convergence loop: max diff in dressed Green (imag time): "<<max_diff<<"\tconverged: "<<converged<<std::endl;*/
			//comparison for the dressed Green's function in Matsubara freq.
			for(uint i=0; i<G_omega.first.size(); i++) {
				if (fabs((G_omega.first[i]-G_old_omega.first[i]).imag()) > max_diff)
					max_diff = fabs((G_omega.first[i]-G_old_omega.first[i]).imag());
				if (fabs((G_omega.second[i]-G_old_omega.second[i]).imag()) > max_diff)
					max_diff = fabs((G_omega.second[i]-G_old_omega.second[i]).imag());	
			}
			std::cerr<<"convergence loop: max diff in dressed Green (Matsubara freq): "<<max_diff<<"\tconverged: "<<converged<<std::endl;
		}
		//write checkpoint - clean this up in the future!
		std::ostringstream G0omega_name, G0tau_name, Gomega_name, Gtau_name, Gtau_name_2, Gomega_name_2;
		G0omega_name<<"G0_omega_"<<iteration_ctr;//<<"_"<<process_id;
		G0tau_name<<"G0_tau_"<<iteration_ctr;//<<"_"<<process_id;
		Gomega_name<<"G_omega_"<<iteration_ctr;//<<"_"<<process_id;
		Gtau_name<<"G_tau_"<<iteration_ctr;//<<"_"<<process_id;
		Gomega_name_2<<"G_omega";//<<"_"<<process_id;
		Gtau_name_2<<"G_tau";//<<"_"<<process_id;
		std::ofstream G0_omega_file(G0omega_name.str().c_str());
		std::ofstream G0_tau_file(G0tau_name.str().c_str());
		std::ofstream G_omega_file(Gomega_name.str().c_str());
		std::ofstream G_tau_file(Gtau_name.str().c_str());
		std::ofstream G_omega_file_2(Gomega_name_2.str().c_str());
		std::ofstream G_tau_file_2(Gtau_name_2.str().c_str());
		assert(G0_omega_file.is_open() && G0_tau_file.is_open() && G_omega_file.is_open() && G_tau_file.is_open());
		print_imag_green_matsubara(G0_omega_file, G0_omega,beta);
		print_green_itime(G0_tau_file,G0_tau,beta);
		print_imag_green_matsubara(G_omega_file,G_omega,beta);
		print_green_itime(G_tau_file,G_tau,beta);
		print_imag_green_matsubara(G_omega_file_2,G_omega,beta);
		print_green_itime(G_tau_file_2,G_tau,beta);
		G0_omega_file.close(); G_omega_file.close(); G0_tau_file.close(); G_tau_file.close(); G_tau_file_2.close(); G_omega_file_2.close();
	}while (max_diff > converged || iteration_ctr <= 1); 
	
}
