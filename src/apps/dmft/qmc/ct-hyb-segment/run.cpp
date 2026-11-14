// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
//
// Library entry point of the segment-picture CT-HYB solver, shared by the
// ct-hyb-segment program, the Python binding and the tests.
#include <alps/solvers.hpp>
#include "hyb.hpp"
#include "hybevaluate.hpp"
#include <alps/utilities/signal.hpp>
#include <ctime>
#include <iostream>
#include <stdexcept>
#ifdef ALPS_HAVE_MPI
#include <alps/mc/mpiadapter.hpp>
#endif

namespace {
#ifdef ALPS_HAVE_MPI
using sim_type = alps::mcmpiadapter<hybridization>;
#else
using sim_type = hybridization;
#endif

//stops the simulation if time > end_time or if signals received.
bool stop_requested(std::time_t end_time) {
  static alps::signal signal;
  return !signal.empty() || std::time(nullptr) > end_time;
}

//post processing: collect Green functions and write them into the output
//file; computes the vertex at the very end
void master_final_tasks(const alps::results_type<hybridization>::type &results,
                        const alps::params &parameters,
                        const std::string &output_name){
  alps::hdf5::archive solver_output(output_name, "a");

  evaluate_basics(results,parameters,solver_output);
  evaluate_gtau(results,parameters,solver_output);
  evaluate_freq(results,parameters,solver_output);
  evaluate_legendre(results,parameters,solver_output);
  evaluate_nnt(results,parameters,solver_output);
  evaluate_nnw(results,parameters,solver_output);
  evaluate_sector_statistics(results,parameters,solver_output);
  evaluate_2p(results, parameters, solver_output);
}
}

int global_mpi_rank;

void alps::solvers::cthyb_segment(alps::params const& input, std::string const& output_file) {
  alps::params parameters(input);
  // Alternative inputs (U or U_MATRIX, MU or MU_VECTOR, ...) are defined
  // without defaults and checked with exists(), so has_missing() does not apply.
  sim_type::define_parameters(parameters);

  const std::time_t deadline = std::time(nullptr) + static_cast<std::time_t>(parameters["MAX_TIME"].as<std::size_t>());
#ifndef ALPS_HAVE_MPI
  global_mpi_rank=0;
  sim_type s(parameters,global_mpi_rank);
#else
  alps::mpi::communicator c;
  c.barrier();
  global_mpi_rank=c.rank();
  sim_type s(parameters, c);
#endif
  if (global_mpi_rank==0)
    std::cout << "Parameters : " << std::endl << parameters << std::endl;

  //run the simulation
  s.run([deadline] { return stop_requested(deadline); });

  //on the master: collect MC results and store them in file, then postprocess
  if (global_mpi_rank==0){
    alps::results_type<hybridization>::type results = collect_results(s);
    std::string output_path = parameters["cthyb.BASEPATH"].as<std::string>()+"/simulation/results";
    alps::save_results(results, parameters, output_file, output_path);
    master_final_tasks(results, parameters, output_file);
  }
#ifdef ALPS_HAVE_MPI
  else{ //on any slave: send back results to master.
    collect_results(s);
  }
  c.barrier();
#endif
}
