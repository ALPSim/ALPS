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

/// @file main.C
/// @brief main program of the DMFT program

#include "hirschfyesim.h"
#include "selfconsistency.h"
#include "externalsolver.h"
#include "hilberttransformer.h"
#include "../DMFT-solvers/Rubtsov/rubtsov.h"

#include <alps/parameters.h>
#include <alps/copyright.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdio.h>



/// @brief The DMFT main program
///
/// The program takes (at the moment) no command line options and reads the parameters 
/// of the simulation from the standard input. Parameters are read in the short ALPS text format where
/// each line contains one parameter definition in the format name=value .
///
/// At the moment the bare imaginary time green's function is assumed to be present in G_input_up and
/// G_input_down.
///
/// The following parameters have special meaning at this time: 
/// 
/// If SOLVER is Hirsch-Fye, a Hirsch-Fye solver is used, otherwise the value of SOLVER is assumed
/// to be the name of an external executable to be used as impurity solver.
/// @todo decent input routines needed!


int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
	try {
#endif
        // we need a factory to create Hirsch-Fye simulations
        alps::scheduler::BasicFactory<HirschFyeSim,HirschFyeRun> factory;
		
		std::cout << "ALPS DMFT solver for the single site impurity problem.\n\n";
		alps::print_copyright(std::cout);
		
		alps::Parameters parms;
		std::cin >> parms;
		//perform selfconsistency loop in...
		if(!parms.defined("OMEGA_LOOP")){
			//...imaginary time tau
			int N=boost::lexical_cast<int>(parms["N"]);
			multiple_vector_type G0_tau, G_tau;
			G0_tau.first.resize(N);
			G0_tau.second.resize(N);
			// read initial guess
			std::ifstream infile_up("G_input_up");		//Bare Green's function in Imaginary time for spin up
			std::ifstream infile_down("G_input_down");	//Bare Green's function in Imaginary time for spin down
			assert(infile_up.is_open() && infile_down.is_open());
			for (int i=0; i<N; i++) {
				infile_up >> G0_tau.first[i];
				infile_down >> G0_tau.second[i];
			}
			
			boost::shared_ptr<ImpuritySolver> solver_ptr;
			if (parms["SOLVER"]=="Hirsch-Fye")
				solver_ptr.reset(new alps::ImpuritySolver(factory,argc,argv));
			else 
				solver_ptr.reset(new ExternalSolver(boost::filesystem::path(parms["SOLVER"],boost::filesystem::native)));
		
			SemicircleHilbertTransformer transform(boost::lexical_cast<double>(parms["t"]));
		
			selfconsistency_loop(parms, *solver_ptr, transform, G0_tau, G_tau);
		
			// write final result
			FILE * outfile_up = fopen("G_output_up","w");
			FILE * outfile_down = fopen("G_output_down","w");
			assert((outfile_up != NULL) && (outfile_down !=NULL));
			for (int i=0; i<N; i++) {
				fprintf(outfile_up, "%5d   %9.6e\n", i, G_tau.first[i]);
				fprintf(outfile_down, "%5d   %9.6e\n", i, G_tau.second[i]);
			}  
		}else{
			//perform self consistency loop in Matsubara frequency omega
			int Nmatsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
			freq_multiple_vector_t G0_omega, G_omega;
			G0_omega.first.resize(Nmatsubara);
			G0_omega.second.resize(Nmatsubara);
			
			// read initial guess
			std::ifstream infile_up("G0omega_input_up");		//Bare Green's function in Matsubara freq for spin up
			std::ifstream infile_down("G0omega_input_down");	//Bare Green's function in Matsubara freq for spin down
			assert(infile_up.is_open() && infile_down.is_open());
			for (int i=0; i<Nmatsubara; i++) {
				assert(!infile_up.eof());
				double first, second;
				infile_up >> first;
				infile_down >> second;
				G0_omega.first[i]=std::complex<double>(0,first);
				G0_omega.second[i]=std::complex<double>(0,second);
			}

			boost::shared_ptr<FrequencyTimeImpuritySolver> solver_ptr;	//construct a solver that takes the GF in both imaginary time and Matsubara freq
			if (parms["SOLVER"]=="Rubtsov")
				solver_ptr.reset(new RubtsovSolver(parms));
			else{
				//throw exception or commit suicide or even better: create an external solver that takes its arguments in frequency space.
				//solver_ptr.reset(new ExternalSolver(boost::filesystem::path(parms["SOLVER"],boost::filesystem::native)));
			}
			FSSemicircleHilbertTransformer transform(boost::lexical_cast<double>(parms["t"]));
			
			selfconsistency_loop_omega_tau(parms, *solver_ptr, transform, G0_omega, G_omega);
				
			// write final result
			std::ofstream outfile_up("Gomega_output_up");
			std::ofstream outfile_down("Gomega_output_down");
			assert(outfile_up.is_open() && outfile_down.is_open());
			for (int i=0; i<Nmatsubara; i++) {
				outfile_up  <<i<<"\t"<<G_omega.first[i]<<std::endl;
				outfile_down<<i<<"\t"<<G_omega.second[i]<<std::endl;
			} 
			outfile_up.close(); outfile_down.close();
		}
#ifndef BOOST_NO_EXCEPTIONS
	}
	catch (std::exception& exc) {
		std::cerr<<exc.what()<<std::endl;
		return -1;
    }
	catch (...) {
		std::cerr << "Fatal Error: Unknown Exception!\n";
		return -2;
	}
#endif  
	return 0;  
}
