/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005-2026 by Emanuel Gull
 *
 * ALPS Project: https://alps.comp-phys.org/
 * SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

/* $Id$ */

/// @file main.cpp
/// @brief main program of the DCA program
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <alps/cli.hpp>
#include <alps/params.hpp>
#include <alps/utilities/citations.hpp>

#include "selfconsistency.h"

void define_dca_parameters(alps::params &p);

void define_all_parameters(alps::params &parms){
  define_selfconsistency_parameters(parms);
  define_solver_parameters(parms);
  define_dca_parameters(parms);
}

alps::params read_parameters_from_file(int argc, const char ** argv){
  alps::params parms(argc, argv);
  define_all_parameters(parms);
  return parms;
}



int main(int argc, char** argv)
{
  try {
    // --help, --license and --citations need no parameter file.
    if (alps::handle_cli_information(argc, argv, "dca", [] {
          alps::params parms;
          parms.description("dca: dynamical cluster approximation self-consistency loop.\n"
                            "usage: dca <parameter file> [--NAME=value ...]");
          define_all_parameters(parms);
          parms.print_help(std::cout);
        }))
      return 0;
    if(argc <2){
      std::cerr<<"use: program_name parameter_file."<<std::endl;
      exit(1);
    }
    std::cout << "#################################################"<<std::endl;
    std::cout << "### ALPS Dynamical Cluster Approximation      ###"<<std::endl;
    std::cout << "### Copyright Emanuel Gull, 2005 - 2026       ###"<<std::endl;
    std::cout << "### Licensed under the MIT license            ###"<<std::endl;
    std::cout << "#################################################"<<std::endl;
    alps::print_citations(std::cout, "dca");

    alps::params parms=read_parameters_from_file(argc, (const char **)argv);
    std::string param_file_name(argv[1]);

    // the --NAME=value overrides after the parameter file, for the solver
    const std::vector<std::string> arguments(argv + 2, argv + argc);
    perform_cluster_loop(param_file_name, parms, arguments);
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
