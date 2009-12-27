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


/// @file main.C
/// @brief main program of the DMFT program

#include "hirschfyesim.h"
#include "selfconsistency.h"
#include "externalsolver.h"
#include "hilberttransformer.h"

#include <alps/parameter.h>
#include <alps/copyright.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdio.h>

#ifdef BUILD_DMFT_QMC_HYBRIDIZATION
#include "hybridization/impurity.h"
#endif



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
      std::cerr<<"use: program_name parameter_file. Under MPI to pass argument to nodes use"<<std::endl;
      std::cerr<<"mpi-run mpi_options program_name -- parameter_file"<<std::endl;
      exit(1);
    }
    std::cout << "ALPS DMFT solver for the single site impurity problem.\n\n";
    alps::print_copyright(std::cout);
    alps::Parameters parms;
    {
      std::ifstream is(argv[1]);
      if(!is.is_open()){std::cerr<<"parameter file argv[1] "<<argv[1]<<"is not open! exiting!"<<std::endl; abort(); }
      is>>parms;
      is.close();
    }
    //perform selfconsistency loop in...
    if(!parms.defined("CLUSTER_LOOP")) {
      if(!parms.defined("OMEGA_LOOP")){
        //...imaginary time tau
        SemicircleHilbertTransformer transform(boost::lexical_cast<double>(parms["t"]));
        boost::shared_ptr<ImpuritySolver> solver_ptr;
        if (parms["SOLVER"]=="Hirsch-Fye"){
          std::cout<<"solving Hirsch Fye"<<std::endl;
          // we need a factory to create Hirsch-Fye simulations
          alps::scheduler::BasicFactory<HirschFyeSim,HirschFyeRun> factory;	
          solver_ptr.reset(new alps::ImpuritySolver(factory,argc,argv));
          selfconsistency_loop(parms, *solver_ptr, transform);
        }
#ifdef BUILD_DMFT_QMC_HYBRIDIZATION
        else if (parms["SOLVER"]=="Werner") {
          std::cout<<"Using Single Site Hybridization Expansion Solver"<<std::endl;
          alps::scheduler::BasicFactory<WernerSimItime,WernerRun> factory;
          solver_ptr.reset(new alps::ImpuritySolver(factory,argc,argv));
          std::cout<<"interpreting input value as G_tau, not G0"<<std::endl;
          itime_green_function_t G_tau = transform.initial_G0(parms);
          F_selfconsistency_loop(parms, *solver_ptr, G_tau);
        }
#endif
        else {
          boost::filesystem::path p(parms["SOLVER"],boost::filesystem::native);
          solver_ptr.reset(new ExternalSolver(boost::filesystem::complete(p)));
          selfconsistency_loop(parms, *solver_ptr, transform);
        }
      } 
      else {
#ifdef RUBTSOV
        alps::scheduler::BasicFactory<RubtsovSim,HalfFillingHubbardRubtsovRun> rubtsov_factory_sshf;
        alps::scheduler::BasicFactory<RubtsovSim,HubbardRubtsovRun> rubtsov_factory_ss; 
        alps::scheduler::BasicFactory<RubtsovSim,MultiBandDensityHubbardRubtsovRun> rubtsov_factory_mbd; 
#endif
#ifdef BUILD_DMFT_QMC_HYBRIDIZATION
        alps::scheduler::BasicFactory<WernerSimFrequency,WernerRun> werner_factory;
#endif
        //perform self consistency loop in Matsubara frequency omega
        FrequencySpaceHilbertTransformer *transform_ptr;
#ifdef IGOR
        if(parms.defined("IGORS_TRANSFORMER")) {
          std::cout<<"using Igor's Hilbert transformer\n";
          transform_ptr = new IgorsHilbertTransformer();
        }
        else 
#endif
          if(!parms.defined("DOSFILE") && !parms.defined("TWODBS")){
            transform_ptr= new FSSemicircleHilbertTransformer(boost::lexical_cast<double>(parms["t"]));
            std::cout<<"using Bethe lattice Hilbert transform"<<std::endl;
          }else if(parms.defined("TWODBS")){
            if(!parms.defined("ANTIFERROMAGNET") || ((bool)(parms.value_or_default("ANTIFERROMAGNET",false))==false)){
              if((bool)(parms.value_or_default("PARAMAGNET",true))==false){ 
                std::cerr<<"AFM and PM? redundant parameters! set paramagnet to false."<<std::endl; abort();
              }
              std::cerr<<"implement PM DOS integration!!"<<std::endl;
              abort();
            }else{
              if((bool)(parms["PARAMAGNET"])==true){ 
                std::cerr<<"AFM and PM? redundant parameters! set ANTIFERROMAGNET to false."<<std::endl; abort(); 
              }
              transform_ptr=new TwoDAFMHilbertTransformer(parms);
            }
          }
          else{
            std::cout<<"using DOS Hilbert transform"<<std::endl;
            if(!parms.defined("ANTIFERROMAGNET") || ((bool)(parms.value_or_default("ANTIFERROMAGNET",false))==false)){
              if((bool)(parms.value_or_default("PARAMAGNET",true))==false){ 
                std::cerr<<"AFM and PM? redundant parameters! set paramagnet to false."<<std::endl; abort();
              }
              transform_ptr= new FSDOSHilbertTransformer(parms);
            }else{
              if((bool)(parms["PARAMAGNET"])==true){ 
                std::cerr<<"AFM and PM? redundant parameters! set ANTIFERROMAGNET to false."<<std::endl; abort(); 
              }
              transform_ptr= new AFM_FSDOSHilbertTransformer(parms);
            }
          }
        boost::shared_ptr<MatsubaraImpuritySolver> solver_ptr;	
#ifdef RUBTSOV
        if ((parms["SOLVER"]=="Rubtsov") && (parms.value_or_default("FLAVORS", "2")=="1")){
          std::cout<<"using single site Hubbard solver for half filling"<<std::endl;
          solver_ptr.reset(new alps::ImpuritySolver(rubtsov_factory_sshf,argc,argv));
        }
        if ((parms["SOLVER"]=="Rubtsov") && (parms.value_or_default("FLAVORS", "2")=="2")){
          std::cout<<"using single site Hubbard solver"<<std::endl;
          solver_ptr.reset(new alps::ImpuritySolver(rubtsov_factory_ss,argc,argv));
        }
        else if ((parms["SOLVER"]=="Rubtsov") && (parms.value_or_default("SITES", "1") =="1")){
          std::cout<<"using multiband Hubbard solver"<<std::endl;
          solver_ptr.reset(new alps::ImpuritySolver(rubtsov_factory_mbd,argc,argv));
        }
#endif
#ifdef BUILD_DMFT_QMC_HYBRIDIZATION
        if (parms["SOLVER"]=="Werner") {
          std::cout<<"Using Single Site Hybridization Expansion Solver"<<std::endl;
          solver_ptr.reset(new alps::ImpuritySolver(werner_factory,argc,argv));
        }
        else{
          boost::filesystem::path p(parms["SOLVER"],boost::filesystem::native);
          solver_ptr.reset(new ExternalSolver(boost::filesystem::complete(p)));
        }
#endif
        selfconsistency_loop_omega(parms, *solver_ptr, *transform_ptr);
      }
    }
    else { //CLUSTER_LOOP
      throw std::logic_error("you should use the cluster framework for a cluster calculation.");
    }		
    {
      std::ofstream os(argv[1]);
      os<<parms;
      os.close();
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
