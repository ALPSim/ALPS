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
#include "../DMFT-solvers/Werner/impurity.h"

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
		
		std::cout << "ALPS DMFT solver for the single site impurity problem.\n\n";
		alps::print_copyright(std::cout);
		
		alps::Parameters parms;
		std::cin >> parms;
		//perform selfconsistency loop in...
		if(!parms.defined("OMEGA_LOOP")){
        		
			//...imaginary time tau
			int n_tau=boost::lexical_cast<int>(parms["N"]);
			int n_orbital=parms.value_or_default("FLAVORS", 2);
			itime_green_function_t G0_tau(n_tau+1, 1, n_orbital), G_tau(n_tau+1, 1, n_orbital);
			// read initial guess
                        G0_tau.read(parms.value_or_default("G0TAU_INPUT", "G0tau_input").c_str());
			std::cout<<"starting value for G is: "<<G0_tau<<std::endl;
			boost::shared_ptr<ImpuritySolver> solver_ptr;
			if (parms["SOLVER"]=="Hirsch-Fye"){
				std::cout<<"solving Hirsch Fye"<<std::endl;
				// we need a factory to create Hirsch-Fye simulations
				alps::scheduler::BasicFactory<HirschFyeSim,HirschFyeRun> factory;	
				solver_ptr.reset(new alps::ImpuritySolver(factory,argc,argv));
				SemicircleHilbertTransformer transform(boost::lexical_cast<double>(parms["t"]));
				selfconsistency_loop(parms, *solver_ptr, transform, G0_tau, G_tau);
			}
			else if (parms["SOLVER"]=="Werner") {
				std::cout<<"solving Werner"<<std::endl;
			    alps::scheduler::BasicFactory<WernerSim,WernerRun> factory;
				solver_ptr.reset(new alps::ImpuritySolver(factory,argc,argv));
				F_SemicircleHilbertTransformer transform(boost::lexical_cast<double>(parms["t"]));
				F_selfconsistency_loop(parms, *solver_ptr, transform, G0_tau, G_tau);
			}
			else {
				solver_ptr.reset(new ExternalSolver(boost::filesystem::path(parms["SOLVER"],boost::filesystem::native)));
		
                                SemicircleHilbertTransformer transform(boost::lexical_cast<double>(parms["t"]));
                                selfconsistency_loop(parms, *solver_ptr, transform, G0_tau, G_tau);
                        }
			// write final result
                        std::ofstream out_file("G_output");
			assert(out_file.is_open());
                        out_file<<G_tau<<std::endl;
		}else{
			alps::scheduler::BasicFactory<RubtsovSim,SSHubbardRubtsovRun> rubtsov_factory_ss; //for single site problems
			alps::scheduler::BasicFactory<RubtsovSim,MultiBandDensityHubbardRubtsovRun> rubtsov_factory_mbd; //for multi orbital problems, but only density density
			//perform self consistency loop in Matsubara frequency omega
			int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
			int n_orbital=parms.value_or_default("FLAVORS", 2);
			matsubara_green_function_t G0_omega(n_matsubara,1,n_orbital), G_omega(n_matsubara,1,n_orbital);
			
			// read initial guess
                        G0_omega.read(parms.value_or_default("G0OMEGA_INPUT", "G0omega_input").c_str());

			boost::shared_ptr<FrequencyTimeImpuritySolver> solver_ptr;	//construct a solver that takes the GF in both imaginary time and Matsubara freq
			if ((parms["SOLVER"]=="Rubtsov") && n_orbital==1){
                          std::cout<<"using single site Hubbard solver"<<std::endl;
				solver_ptr.reset(new alps::ImpuritySolver(rubtsov_factory_ss,argc,argv));
			}
                        else if ((parms["SOLVER"]=="Rubtsov") && (parms.value_or_default("SITES", "1") =="1")){
                          std::cout<<"using multiband Hubbard solver"<<std::endl;
				solver_ptr.reset(new alps::ImpuritySolver(rubtsov_factory_mbd,argc,argv));
			}
			else{
				//throw exception or commit suicide or even better: create an external solver that takes its arguments in frequency space.
				//solver_ptr.reset(new ExternalSolver(boost::filesystem::path(parms["SOLVER"],boost::filesystem::native)));
				boost::throw_exception(std::logic_error("there is no such solver!"));	
			}
			FrequencySpaceHilbertTransformer *transform_ptr;
                        if(!parms.defined("DOSFILE")){
                              transform_ptr= new FSSemicircleHilbertTransformer(boost::lexical_cast<double>(parms["t"]));
                              std::cout<<"using Bethe lattice Hilbert transform"<<std::endl;
                        }else{
                              std::cout<<"using DOS Hilbert transform"<<std::endl;
                              transform_ptr= new FSDOSHilbertTransformer(parms);
                        }
			selfconsistency_loop_omega_tau(parms, *solver_ptr, *transform_ptr, G0_omega, G_omega);
				
			// write final result
			std::ofstream outfile("Gomega_output");
			assert(outfile.is_open());
                        outfile<<G_omega<<std::endl;
			outfile.close();
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
