/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005-2014 by Emanuel Gull <gullc@comp-phys.org>
 *
 * To be used with permission by the authors
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

void read_parameters_from_file(alps::Parameters &parms, const std::string &filename){
  std::ifstream is(filename.c_str());
  if(!is.is_open()){std::cerr<<"parameter file: "<<filename<<" is not open! exiting!"<<std::endl; exit(1); }
  is>>parms;
  is.close();
  parms["BASENAME"]=std::string(filename);
}



int main(int argc, char** argv)
{
  try {
    if(argc <2){
      std::cerr<<"use: program_name parameter_file. Under MPI to pass argument to nodes use"<<std::endl;
      std::cerr<<"mpi-run mpi_options program_name -- parameter_file"<<std::endl;
      exit(1);
    }
    std::cout << "ALPS DMFT self-consistency condition.\n\n";
    alps::print_copyright(std::cout);
    alps::Parameters parms;
    read_parameters_from_file(parms, argv[1]);
    
    //perform selfconsistency loop for single site or cluster...
    if(!parms.defined("CLUSTER_LOOP"))
      perform_single_site_loop(parms);
    else
      perform_cluster_loop(parms);
    if(parms.defined("MODIFY_PARAMETER_FILE")){
      std::ofstream os(argv[1]);
	
      (*const_cast<alps::Parameters*>(&parms))["CURRENT_IT"]=0;// sloppy fix.  Sets current_it to zero for continuing program later
      os<<parms;
      os.close();
    }
  }
  catch (std::exception& exc) {
    std::cerr<<exc.what()<<std::endl;
    return -1;
  }
  catch (...) {
    std::cerr << "Fatal Error: Unknown Exception!\n";
    return -2;
  }
  return 0;
}
