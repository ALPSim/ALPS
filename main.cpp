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
#include <iostream>
#include <fstream>

#include <alps/params.hpp>

#include "selfconsistency.h"

/// @brief The DMFT main program

alps::params read_parameters_from_file(const std::string &filename){
  {
    std::ifstream is(filename.c_str());
    if(!is.is_open()){std::cerr<<"parameter file: "<<filename<<" is not open! exiting!"<<std::endl; exit(1); }
  }
  alps::params parms(filename);
  parms["BASENAME"]=std::string(filename);
  return parms;
}



int main(int argc, char** argv)
{
  try {
    if(argc <2){
      std::cerr<<"use: program_name parameter_file. Under MPI to pass argument to nodes use"<<std::endl;
      std::cerr<<"mpi-run mpi_options program_name -- parameter_file"<<std::endl;
      exit(1);
    }
    std::cout << "#################################################"<<std::endl;
    std::cout << "### ALPS Dynamical Mean Field Framework.      ###"<<std::endl;
    std::cout << "### Copyright Emanuel Gull, 2005 - 2015       ###"<<std::endl;
    std::cout << "### To be used with permission by the authors ###"<<std::endl;
    std::cout << "#################################################"<<std::endl;
    alps::params parms=read_parameters_from_file(argv[1]);
    //perform selfconsistency loop for single site or cluster...
    if(!parms.defined("CLUSTER_LOOP")){
      perform_single_site_loop(parms);
    }else{
//      perform_cluster_loop(parms);
    }
    if(parms.defined("MODIFY_PARAMETER_FILE")){
      std::ofstream os(argv[1]);
	
      parms["CURRENT_IT"]=0;// sloppy fix.  Sets current_it to zero for continuing program later
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
