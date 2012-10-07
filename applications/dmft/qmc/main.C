/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *               2012           Jakub Imriska <jimriska@phys.ethz.ch>
 *
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
    chdir(p.branch_path().string().c_str());

    //perform selfconsistency loop in...
    if(!parms.defined("CLUSTER_LOOP")) {
      if(!parms.defined("OMEGA_LOOP") || (bool)(parms["OMEGA_LOOP"])==false){
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
        }
        else {
            /*boost::filesystem::path*/ std::string p(parms["SOLVER"]/**/);
          
          solver_ptr.reset(new ExternalSolver(/*boost::filesystem::absolute(*/p/*)*/));
          selfconsistency_loop(parms, *solver_ptr, transform);
        }
      } 
      else {
        alps::scheduler::BasicFactory<InteractionExpansionSim,HalfFillingHubbardInteractionExpansionRun> interaction_expansion_factory_sshf;
        alps::scheduler::BasicFactory<InteractionExpansionSim,HubbardInteractionExpansionRun> interaction_expansion_factory_ss; 
        alps::scheduler::BasicFactory<InteractionExpansionSim,MultiBandDensityHubbardInteractionExpansionRun> interaction_expansion_factory_mbd; 
        //perform self consistency loop in Matsubara frequency omega
        FrequencySpaceHilbertTransformer *transform_ptr;
        if(!parms.defined("DOSFILE") && !parms.defined("TWODBS")){
          transform_ptr= new FSSemicircleHilbertTransformer(parms);
          std::cout<<"using Bethe lattice Hilbert transform"<<std::endl;
        } 
        else if(parms.defined("TWODBS")) {
          if(!parms.defined("ANTIFERROMAGNET") || ((bool)(parms.value_or_default("ANTIFERROMAGNET",false))==false)){
            if((bool)(parms["SYMMETRIZATION"])==false){ 
              std::cerr<<"ANTIFERROMAGNET==false  and  SYMMETRIZATION==false? redundant parameters! Set SYMMETRIZATION=1 (true) for a paramagnetic solution."<<std::endl; abort();
            }
            transform_ptr=new TwoDHilbertTransformer(parms);
          }else{
            if((bool)(parms["SYMMETRIZATION"])==true){ 
              std::cerr<<"ANTIFERROMAGNET==true  and  SYMMETRIZATION==true? redundant parameters! Set SYMMETRIZATION=0 (false) for an antiferromagnetic solution."<<std::endl; abort(); 
            }
            transform_ptr=new TwoDAFMHilbertTransformer(parms);
          }
        }
        else{
          std::cout<<"using DOS Hilbert transform"<<std::endl;
          if(!parms.defined("ANTIFERROMAGNET") || ((bool)(parms.value_or_default("ANTIFERROMAGNET",false))==false)){
            if((bool)(parms["SYMMETRIZATION"])==false){ 
              std::cerr<<"ANTIFERROMAGNET==false  and  SYMMETRIZATION==false? redundant parameters! Set SYMMETRIZATION=1 (true) for a paramagnetic solution."<<std::endl; abort();
            }
            transform_ptr= new FSDOSHilbertTransformer(parms);
          }else{
            if((bool)(parms["SYMMETRIZATION"])==true){ 
              std::cerr<<"ANTIFERROMAGNET==true  and  SYMMETRIZATION==true? redundant parameters! Set SYMMETRIZATION=0 (false) for an antiferromagnetic solution."<<std::endl; abort(); 
            }
            transform_ptr= new AFM_FSDOSHilbertTransformer(parms);
          }
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
          else
          {
            std::string p(parms["SOLVER"]);
            std::cout<<"using external solver: "<<p<<std::endl;
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
