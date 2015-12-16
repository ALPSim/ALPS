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
void define_parameters(alps::params &p){
  p.define<std::string>("SC_LOOP", "switch between single_site, DCA, and CDMFT");
}
void define_dca_parameters(alps::params &p);

alps::params read_parameters_from_file(int argc, const char ** argv){
  alps::params parms(argc, argv);
  define_selfconsistency_parameters(parms);
  define_solver_parameters(parms);
  define_dca_parameters(parms);
  define_parameters(parms);

  if (parms.help_requested(std::cout)) {
    exit(0);
  }
  return parms;
}



int main(int argc, char** argv)
{
  try {
    if(argc <2){
      std::cerr<<"use: program_name parameter_file. Under MPI to pass argument to nodes use"<<std::endl;
      std::cerr<<"mpirun mpi_options program_name -- parameter_file"<<std::endl;
      exit(1);
    }
    std::cout << "#################################################"<<std::endl;
    std::cout << "### ALPS Dynamical Mean Field Framework.      ###"<<std::endl;
    std::cout << "### Copyright Emanuel Gull, 2005 - 2015       ###"<<std::endl;
    std::cout << "### To be used with permission by the authors ###"<<std::endl;
    std::cout << "#################################################"<<std::endl;

    alps::params parms=read_parameters_from_file(argc, (const char **)argv);
    std::string param_file_name(argv[1]);

    if(parms["SC_LOOP"].as<std::string>()=="SINGLE_SITE"){
      perform_single_site_loop(param_file_name, parms);
    }else if(parms["SC_LOOP"].as<std::string>()=="DCA"){
      perform_cluster_loop(param_file_name, parms);
    }else
      throw std::invalid_argument("self-consistency loop: "+parms["SC_LOOP"].as<std::string>()+"unknown");
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
