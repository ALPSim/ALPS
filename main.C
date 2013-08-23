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

#include "selfconsistency.h"
#include "externalsolver.h"
#include "hilberttransformer.h"
#include "cluster.h"

#include <alps/parameter.h>
#include <alps/utility/copyright.hpp>

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdio.h>



/// @brief The DMFT main program
///
/// The program takes (at the moment) no command line options and reads the parameters
/// of the simulation from the standard input. Parameters are read in the short ALPS text format where
/// each line contains one parameter definition in the format name=value .


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
    std::cout << "ALPS DMFT self-consistency condition.\n\n";
    alps::print_copyright(std::cout);
    alps::Parameters parms;
    {
      std::ifstream is(argv[1]);
      if(!is.is_open()){std::cerr<<"parameter file argv[1]: "<<argv[1]<<" is not open! exiting!"<<std::endl; abort(); }
      is>>parms;
      is.close();
      parms["BASENAME"]=std::string(argv[1]);
    }
    //perform selfconsistency loop in...
    if(!parms.defined("CLUSTER_LOOP")) {
      //perform self consistency loop in Matsubara frequency omega
      FrequencySpaceHilbertTransformer *transform_ptr;
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
      selfconsistency_loop_omega(parms, *solver_ptr, *transform_ptr);
    }
    else { //CLUSTER_LOOP
      if (parms["CLUSTER_LOOP"]=="DCA") {
        if(parms.defined("NAMBU") && parms["NAMBU"]){
          std::cout<<"preparing a NAMBU DMFT calculation"<<std::endl;
          boost::shared_ptr<DCATransformer> transform_ptr;
          DCATransformer::generate_transformer(parms, transform_ptr);
          boost::shared_ptr<NambuImpuritySolver> solver_ptr;
          std::cout<<"setting up external NAMBU solver: "<<parms["SOLVER"]<<std::endl;
          solver_ptr.reset(new ExternalSolver(boost::filesystem::path(parms["SOLVER"])));
          selfconsistency_loop_DCA(parms, *solver_ptr, *transform_ptr);
        }else{
          boost::shared_ptr<DCATransformer> transform_ptr;
          DCATransformer::generate_transformer(parms, transform_ptr);
          boost::shared_ptr<MatsubaraImpuritySolver> solver_ptr;
          std::cout<<"starting external solver "<<parms["SOLVER"]<<std::endl;
          solver_ptr.reset(new ExternalSolver(boost::filesystem::path(parms["SOLVER"])));
          selfconsistency_loop_DCA(parms, *solver_ptr, *transform_ptr);
        }
      }else if (parms["CLUSTER_LOOP"]=="CDMFT") {
        if(parms.defined("NAMBU") && parms["NAMBU"]){
          std::cout<<"preparing a NAMBU DMFT calculation"<<std::endl;
          boost::shared_ptr<CDMFTTransformer> transform_ptr;
          CDMFTTransformer::generate_transformer(parms, transform_ptr);
          boost::shared_ptr<NambuImpuritySolver> solver_ptr;
          std::cout<<"setting up external NAMBU solver: "<<parms["SOLVER"]<<std::endl;
          solver_ptr.reset(new ExternalSolver(boost::filesystem::path(parms["SOLVER"])));
          selfconsistency_loop_CDMFT(parms, *solver_ptr, *transform_ptr);
        }else{
          boost::shared_ptr<CDMFTTransformer> transform_ptr;
          CDMFTTransformer::generate_transformer(parms, transform_ptr);
          boost::shared_ptr<MatsubaraImpuritySolver> solver_ptr;
          std::cout<<"starting external solver "<<parms["SOLVER"]<<std::endl;
          solver_ptr.reset(new ExternalSolver(boost::filesystem::path(parms["SOLVER"])));
          selfconsistency_loop_CDMFT(parms, *solver_ptr, *transform_ptr);
        }
      }
      else {
        std::cout<<parms["CLUSTER_LOOP"]<<std::endl;
        boost::throw_exception(std::logic_error("there is no such cluster method!"));
      }
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
