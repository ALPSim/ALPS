/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005-2007 by Philipp Werner <werner@comp-phys.org>,
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

/// @file main.C
/// @brief main program of the DMFT program

#include "hirschfyesim.h"
#include "selfconsistency.h"
#include "externalsolver.h"
#include "hilberttransformer.h"
#include "cluster.h"
#ifdef RUBTSOV
#include "../DMFT-solvers/Rubtsov/rubtsov.h"
#endif
#ifdef WERNER
#include "../DMFT-solvers/Werner/impurity.h"
#endif

#include <alps/parameter.h>
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
    if(argc <2){
      std::cout<<"use: program_name parameter_file. Under MPI to pass argument to nodes use"<<std::endl;
      std::cout<<"mpi-run mpi_options program_name -- parameter_file"<<std::endl;
      exit(1);
    }
    std::cout << "ALPS DMFT solver for the single site impurity problem.\n\n";
    alps::print_copyright(std::cout);
    alps::Parameters parms;
    {
      std::ifstream is(argv[1]);
      assert(is.is_open());
      is>>parms;
    }
    //std::cout<<parms<<std::endl;
    boost::shared_ptr<FourierTransformer> fourier_ptr;
    FourierTransformer::generate_transformer(parms, fourier_ptr);
    //perform selfconsistency loop in...
    if(!parms.defined("CLUSTER_LOOP")) {
      if(!parms.defined("OMEGA_LOOP")){
	//...imaginary time tau
	SemicircleHilbertTransformer transform(boost::lexical_cast<double>(parms["t"]));
	itime_green_function_t G0_tau = transform.initial_G0(parms);
	itime_green_function_t G_tau = G0_tau;
	boost::shared_ptr<ImpuritySolver> solver_ptr;
	if (parms["SOLVER"]=="Hirsch-Fye"){
	  std::cout<<"solving Hirsch Fye"<<std::endl;
	  // we need a factory to create Hirsch-Fye simulations
	  alps::scheduler::BasicFactory<HirschFyeSim,HirschFyeRun> factory;	
	  solver_ptr.reset(new alps::ImpuritySolver(factory,argc,argv));
	  selfconsistency_loop(parms, *solver_ptr, transform, *fourier_ptr, G0_tau, G_tau);
	}
#ifdef WERNER
	else if (parms["SOLVER"]=="Werner") {
	  std::cout<<"solving Werner"<<std::endl;
	  alps::scheduler::BasicFactory<WernerSimItime,WernerRun> factory;
	  solver_ptr.reset(new alps::ImpuritySolver(factory,argc,argv));
	  std::cout<<"interpreting input value as G_tau, not G0"<<std::endl;
	  //G_tau=G0_tau;
	  F_selfconsistency_loop(parms, *solver_ptr, G_tau);
	}
#endif
	else {
	  solver_ptr.reset(new ExternalSolver(boost::filesystem::path(parms["SOLVER"],boost::filesystem::native)));
	  selfconsistency_loop(parms, *solver_ptr, transform, *fourier_ptr, G0_tau, G_tau);
	}
	// write G0 (to be read in as an input for a new simulation)
	G0_tau.write(parms.value_or_default("G0OMEGA_output", "G0omega_output").c_str());
      } 
      else {
#ifdef RUBTSOV
	alps::scheduler::BasicFactory<RubtsovSim,HalfFillingHubbardRubtsovRun> rubtsov_factory_ss; //for single site problems
	alps::scheduler::BasicFactory<RubtsovSim,MultiBandDensityHubbardRubtsovRun> rubtsov_factory_mbd; 
	//for multi orbital problems, but only density density
#endif
#ifdef WERNER
	alps::scheduler::BasicFactory<WernerSimFrequency,WernerRun> werner_factory;
#endif
	//perform self consistency loop in Matsubara frequency omega
	FrequencySpaceHilbertTransformer *transform_ptr;
	if(!parms.defined("DOSFILE")){
	  transform_ptr= new FSSemicircleHilbertTransformer(boost::lexical_cast<double>(parms["t"]));
	  std::cout<<"using Bethe lattice Hilbert transform"<<std::endl;
	}
	else{
	  std::cout<<"using DOS Hilbert transform"<<std::endl;
	  transform_ptr= new FSDOSHilbertTransformer(parms);
	}
	matsubara_green_function_t G0_omega = transform_ptr->initial_G0(parms);
	matsubara_green_function_t G_omega = G0_omega;
	boost::shared_ptr<MatsubaraImpuritySolver> solver_ptr;	
#ifdef RUBTSOV
	if ((parms["SOLVER"]=="Rubtsov") && (parms.value_or_default("FLAVORS", "2")=="1")){
	  std::cout<<"using single site Hubbard solver"<<std::endl;
	  solver_ptr.reset(new alps::ImpuritySolver(rubtsov_factory_ss,argc,argv));
	}
	else if ((parms["SOLVER"]=="Rubtsov") && (parms.value_or_default("SITES", "1") =="1")){
	  std::cout<<"using multiband Hubbard solver"<<std::endl;
	  solver_ptr.reset(new alps::ImpuritySolver(rubtsov_factory_mbd,argc,argv));
	}
#endif
#ifdef WERNER
	else if (parms["SOLVER"]=="Werner") {
	  std::cout<<"solving Werner -- using omega loop"<<std::endl;
	  solver_ptr.reset(new alps::ImpuritySolver(werner_factory,argc,argv));
	}
	else{
	  //throw exception or commit suicide or even better: 
	  //create an external solver that takes its arguments in frequency space.
	  solver_ptr.reset(new ExternalSolver(boost::filesystem::path(parms["SOLVER"],boost::filesystem::native)));
	  //std::cout<<parms["SOLVER"]<<std::endl;
	  //boost::throw_exception(std::logic_error("there is no such solver!"));	
	}
#endif
	selfconsistency_loop_omega(parms, *solver_ptr, *transform_ptr, *fourier_ptr, G0_omega, G_omega);
      
	// write G0 (to be read as an input Green function)
	G0_omega.write(parms.value_or_default("G0OMEGA_output", "G0omega_output").c_str());
      }
    }
    else { //CLUSTER_LOOP
      if (parms["CLUSTER_LOOP"]=="DCA") {
#ifdef RUBTSOV
	alps::scheduler::BasicFactory<RubtsovSim,HalfFillingHubbardRubtsovRun> half_fill_rubtsov_factory; 
	alps::scheduler::BasicFactory<RubtsovSim,HubbardRubtsovRun> rubtsov_factory;
#endif
	//perform self consistency loop in Matsubara frequency omega
	// all gfs are formulated in k-space and are diagonal in k! 
	DCATransformer *transform_ptr;
	std::cout<<"using simple DCA transformer"<<std::endl;
	if (parms.defined("ANTIFERROMAGNET"))
	  transform_ptr = new AntiFerroDCATransformer(parms);
	else
	  transform_ptr = new SimpleDCATransformer(parms);
	matsubara_green_function_t G0_omega=transform_ptr->initial_G0(parms);
	matsubara_green_function_t G_omega=G0_omega;
	boost::shared_ptr<MatsubaraImpuritySolver> solver_ptr;	
#ifdef RUBTSOV
	if (parms["SOLVER"]=="Rubtsov" && (parms.value_or_default("FLAVORS", "2")=="1")){
	  std::cout<<"using Rubtsov Solver for half filling"<<std::endl;
	  solver_ptr.reset(new alps::ImpuritySolver(half_fill_rubtsov_factory,argc,argv));
	}
	else if (parms["SOLVER"]=="Rubtsov" && (parms.value_or_default("FLAVORS", "2")=="2")){
	  std::cout<<"using Rubtsov Solver"<<std::endl;
	  solver_ptr.reset(new alps::ImpuritySolver(rubtsov_factory,argc,argv));
	}
	else{
#endif
          std::cout<<"starting external solver "<<parms["SOLVER"]<<std::endl;
	  solver_ptr.reset(new ExternalSolver(boost::filesystem::path(parms["SOLVER"],boost::filesystem::native)));
	  //std::cout<<parms["SOLVER"]<<std::endl;
	  //boost::throw_exception(std::logic_error("there is no such cluster solver!"));	
#ifdef RUBTSOV
	}
#endif
	selfconsistency_loop_DCA(parms, *solver_ptr, *transform_ptr, *fourier_ptr, G0_omega, G_omega);
	// write G0
	G0_omega.write(parms.value_or_default("G0OMEGA_output", "G0omega_output").c_str());
      }
      else {
	std::cout<<parms["CLUSTER_LOOP"]<<std::endl;
	boost::throw_exception(std::logic_error("there is no such cluster method!"));
      }
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
