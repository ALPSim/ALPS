/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2026 by Emanuel Gull <gull@phys.columbia.edu>
 *               2005 - 2009 by Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *               2012 - 2013 by Jakub Imriska <jimriska@phys.ethz.ch>
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/


/// @file main.C
/// @brief main program of the DMFT program

#include "hirschfyesim.h"
#include "selfconsistency.h"
#include "externalsolver.h"
#include "hilberttransformer.h"
#include "interaction_expansion_choice.h"

#include <alps/parameter.h>
#include <alps/utility/copyright.hpp>
#include <alps/utility/citations.hpp>
#include <alps/utility/cli.hpp>

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdio.h>
#include <boost/filesystem/operations.hpp>

#include "interaction_expansion/interaction_expansion.hpp"



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
    if (alps::handle_cli_information(argc, argv, "dmft", [&] {
      std::cout << "Usage: " << argv[0] << " parameter_file\n";
    })) return 0;
    if (argc != 2) throw std::invalid_argument("Expected one parameter file");
    alps::cli_mpi_guard mpi(argc, argv);
    alps::Parameters parms;
    {
      std::ifstream is(argv[1]);
      if (!is.is_open()) throw std::runtime_error(std::string("Cannot open parameter file: ") + argv[1]);
      is>>parms;
      parms["BASENAME"]=std::string(argv[1]);
    }
    if (alps::cli_is_master()) {
      std::cout << "ALPS DMFT framework for the single site impurity problem.\n\n";
      alps::print_copyright(std::cout, "dmft");
    }
    // set working directory
    boost::filesystem::path p(static_cast<std::string>(parms["BASENAME"]));
    if (!p.parent_path().empty())
      boost::filesystem::current_path(p.parent_path());

    //perform selfconsistency loop in...
    if(!parms.defined("CLUSTER_LOOP")) {
      if(!parms.value_or_default("OMEGA_LOOP",false)){
        //...imaginary time tau
        SemicircleHilbertTransformer transform(parms);
        boost::shared_ptr<ImpuritySolver> solver_ptr;
        if (parms["SOLVER"]=="Hirsch-Fye"){
          std::cout<<"solving Hirsch Fye"<<std::endl;
          // we need a factory to create Hirsch-Fye simulations
          alps::scheduler::BasicFactory<HirschFyeSim,HirschFyeRun> factory;  
          solver_ptr.reset(new alps::ImpuritySolver(factory,argc,argv));
          // Built-in solvers print no notice of their own; external solvers do.
          if (alps::cli_is_master()) alps::print_citations(std::cout, "hirschfye");
          selfconsistency_loop(parms, *solver_ptr, transform);
        }
        else if (parms["SOLVER"]=="Hybridization") {
          throw std::invalid_argument("The internal hybridization solver has been replaced by a standalone hybridzation solver.\nPlease use the \'hybridization\' program");
        } else if (parms["SOLVER"]=="hybridization" || static_cast<bool>(parms.value_or_default("SC_WRITE_DELTA",false))) {
          std::string p(parms["SOLVER"]);
          std::cout<<"using external solver: "<<p<<std::endl;
          parms["SC_WRITE_DELTA"]=1; //we need the hybridization function for this solver
          solver_ptr.reset(new ExternalSolver(p));
          F_selfconsistency_loop(parms, *solver_ptr, transform);
        } else {
            /*boost::filesystem::path*/ std::string p(parms["SOLVER"]/**/);
          solver_ptr.reset(new ExternalSolver(/*boost::filesystem::absolute(*/p/*)*/));
          selfconsistency_loop(parms, *solver_ptr, transform);
        }
      }
      else {
        //perform self consistency loop in Matsubara frequency omega
        if (parms.value_or_default("ANTIFERROMAGNET",false) && parms.value_or_default("SYMMETRIZATION",false)) { 
          std::cerr<<"ERROR: incompatible parameters: ANTIFERROMAGNET==1 (true)  and  SYMMETRIZATION==1 (true)"<<std::endl
                   <<"  for simulation in PM phase set: SYMMETRIZATION=1, ANTIFERROMAGNET=0."<<std::endl
                   <<"  for simulation in FM phase set: SYMMETRIZATION=0, ANTIFERROMAGNET=0."<<std::endl
                   <<"  for simulation in AFM phase set: SYMMETRIZATION=0, ANTIFERROMAGNET=1."<<std::endl;
          throw std::logic_error("Incompatible parameters: ANTIFERROMAGNET==1,  SYMMETRIZATION==1.");
        }
          
        alps::scheduler::BasicFactory<InteractionExpansionSim,HubbardInteractionExpansionRun> interaction_expansion_factory_ss;
        alps::scheduler::BasicFactory<InteractionExpansionSim,MultiBandDensityHubbardInteractionExpansionRun> interaction_expansion_factory_mbd;

        boost::shared_ptr<FrequencySpaceHilbertTransformer> transform_ptr;
        if (parms.value_or_default("SEMICIRCLE_HILBERT",false)) {
          transform_ptr.reset(new SemicircleFSHilbertTransformer(parms));
        } else {
          transform_ptr.reset(new GeneralFSHilbertTransformer(parms));
        }
        
        // Exactly one impurity solver is selected. The CT-INT run type is
        // chosen by select_interaction_expansion (FLAVORS/SITES read as ints);
        // every other SOLVER value is handed to the external/hybridization
        // path. No solver is constructed-then-discarded.
        boost::shared_ptr<MatsubaraImpuritySolver> solver_ptr;
        if (parms["SOLVER"]=="Interaction Expansion") {
          const int flavors = static_cast<int>(parms.value_or_default("FLAVORS", 2));
          const int sites   = static_cast<int>(parms.value_or_default("SITES", 1));
          switch (select_interaction_expansion(flavors, sites)) {
            case interaction_expansion_choice::single_site_hubbard:
              std::cout<<"using single site Hubbard solver"<<std::endl;
              solver_ptr.reset(new alps::ImpuritySolver(interaction_expansion_factory_ss,argc,argv));
              break;
            case interaction_expansion_choice::multiband_density:
              std::cout<<"using multiband Hubbard solver"<<std::endl;
              solver_ptr.reset(new alps::ImpuritySolver(interaction_expansion_factory_mbd,argc,argv));
              break;
            case interaction_expansion_choice::unsupported:
              throw std::runtime_error("DMFT Interaction Expansion: unsupported (FLAVORS, SITES) "
                                       "combination; set FLAVORS=2 (single site) or SITES=1 (multiband).");
          }
          if (alps::cli_is_master()) alps::print_citations(std::cout, "interaction");
        }
        else if (parms["SOLVER"]=="Hybridization") {
          throw std::invalid_argument("The internal hybridization solver has been replaced by a standalone hybridzation solver.\nPlease use the \'hybridization\' program");
        }
        else {
          std::string p(parms["SOLVER"]);
          std::cout<<"using external solver: "<<p<<std::endl;
          if(parms["SOLVER"]=="hybridization") parms["SC_WRITE_DELTA"]=1; //we need the hybridization function for this solver
          solver_ptr.reset(new ExternalSolver(p));
        }
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
    {
      alps::hdf5::archive os(std::string(argv[1])+".h5", "a");
      os<<alps::make_pvp("/parameters",parms);
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
