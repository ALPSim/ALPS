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

/// @file main.cpp
/// @brief main program of the DCA program
#include <iostream>
#include <fstream>

#include <alps/params.hpp>

#include "selfconsistency.h"

void define_dca_parameters(alps::params &p);

alps::params read_parameters_from_file(int argc, const char ** argv){
  alps::params parms(argc, argv);
  define_selfconsistency_parameters(parms);
  define_solver_parameters(parms);
  define_dca_parameters(parms);

  if (parms.help_requested(std::cout)) {
    exit(0);
  }
  return parms;
}



int main(int argc, char** argv)
{
  try {
    if(argc <2){
      std::cerr<<"use: program_name parameter_file."<<std::endl;
      exit(1);
    }
    std::cout << "#################################################"<<std::endl;
    std::cout << "### ALPS Dynamical Mean Field Framework.      ###"<<std::endl;
    std::cout << "### Copyright Emanuel Gull, 2005 - 2016       ###"<<std::endl;
    std::cout << "### To be used with permission by the authors ###"<<std::endl;
    std::cout << "#################################################"<<std::endl;

    alps::params parms=read_parameters_from_file(argc, (const char **)argv);
    std::string param_file_name(argv[1]);

    perform_cluster_loop(param_file_name, parms);
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
