/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *               2012 - 2013 by Jakub Imriska <jimriska@phys.ethz.ch>
 *
 *
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/


/// @file main.C
/// @brief main program of the DMFT program

#include "hirschfyesim.h"
#include "selfconsistency.h"
#include "externalsolver.h"
#include "hilberttransformer.h"

#include <alps/parameter.h>
#include <alps/utility/copyright.hpp>

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdio.h>
#ifdef BOOST_MSVC
#include <direct.h>
#endif

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
    if(argc <2){
      std::cerr<<"use: program_name parameter_file. Under MPI to pass argument to nodes use"<<std::endl;
      std::cerr<<"mpi-run mpi_options program_name -- parameter_file"<<std::endl;
      exit(1);
    }
    std::cout << "ALPS DMFT framework for the single site impurity problem.       "<<std::endl;
    std::cout << "  For further information see the ALPS DMFT paper:              "<<std::endl;
    std::cout << "  Computer Physics Communications 182, 1078 (2011)              "<<std::endl;
    std::cout << "                                                                "<<std::endl;
    std::cout << "  copyright (c) 2005-2010 by the ALPS collaboration.            "<<std::endl;
    std::cout << "  Consult the web page for license details.                     "<<std::endl;
    std::cout << std::endl;
    alps::print_copyright(std::cout);
    alps::Parameters parms;
    {
      std::ifstream is(argv[1]);
      if(!is.is_open()){std::cerr<<"parameter file argv[1] "<<argv[1]<<" is not open! exiting!"<<std::endl; abort(); }
      is>>parms;
      parms["BASENAME"]=std::string(argv[1]);
    }
    // set working directory
    boost::filesystem::path p(static_cast<std::string>(parms["BASENAME"]));
    chdir(p.parent_path().string().c_str());

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
          
        alps::scheduler::BasicFactory<InteractionExpansionSim,HalfFillingHubbardInteractionExpansionRun> interaction_expansion_factory_sshf;
        alps::scheduler::BasicFactory<InteractionExpansionSim,HubbardInteractionExpansionRun> interaction_expansion_factory_ss; 
        alps::scheduler::BasicFactory<InteractionExpansionSim,MultiBandDensityHubbardInteractionExpansionRun> interaction_expansion_factory_mbd; 
        
        boost::shared_ptr<FrequencySpaceHilbertTransformer> transform_ptr;
        if (parms.value_or_default("SEMICIRCLE_HILBERT",false)) {
          transform_ptr.reset(new SemicircleFSHilbertTransformer(parms));
        } else {
          transform_ptr.reset(new GeneralFSHilbertTransformer(parms));
        }
        
        boost::shared_ptr<MatsubaraImpuritySolver> solver_ptr;  
        if ((parms["SOLVER"]=="Interaction Expansion") && (parms.value_or_default("FLAVORS", "2")=="1")){
          std::cout<<"using single site Hubbard solver for half filling"<<std::endl;
          solver_ptr.reset(new alps::ImpuritySolver(interaction_expansion_factory_sshf,argc,argv));
        }
        if ((parms["SOLVER"]=="Interaction Expansion") && (parms.value_or_default("FLAVORS", "2")=="2")){
          std::cout<<"using single site Hubbard solver"<<std::endl;
          solver_ptr.reset(new alps::ImpuritySolver(interaction_expansion_factory_ss,argc,argv));
        }
        else if ((parms["SOLVER"]=="Interaction Expansion") && (parms.value_or_default("SITES", "1") =="1")){
          std::cout<<"using multiband Hubbard solver"<<std::endl;
          solver_ptr.reset(new alps::ImpuritySolver(interaction_expansion_factory_mbd,argc,argv));
        }
        else
          if (parms["SOLVER"]=="Hybridization") {
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
