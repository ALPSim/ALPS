/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                       Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
 *
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#include <ctime>
#include "hyb.hpp"
#include "hybevaluate.hpp"
#include <boost/bind/bind.hpp>

#ifdef ALPS_HAVE_MPI
#include <alps/mc/mpiadapter.hpp>
typedef alps::mcmpiadapter<hybridization> sim_type;
#else
typedef hybridization sim_type;
#endif

//stops the simulation if time > end_time or if signals received.
bool stop_callback(const time_t end_time) {
  static alps::signal signal;
  return !signal.empty() || time(0) > end_time;
}
void master_final_tasks(const alps::results_type<hybridization>::type &results, const alps::parameters_type<hybridization>::type &parameters, const std::string &output_name);
int global_mpi_rank;

int main(int argc, char* argv[]){

#ifdef ALPS_HAVE_MPI
    //boot up MPI environment
    MPI_Init(&argc, &argv);
#endif

  //read in command line options
  alps::parameters_type<hybridization>::type parameters(argc, (const char**)argv, "/parameters");
  sim_type::define_parameters(parameters);
  if (parameters.help_requested(std::cout)) {
    exit(0);
  }

    try {
      unsigned long max_time=parameters["MAX_TIME"];

#ifndef ALPS_HAVE_MPI
      global_mpi_rank=0;
      sim_type s(parameters,global_mpi_rank);
#else
      alps::mpi::communicator c;
      c.barrier();
      global_mpi_rank=c.rank();
      if (!global_mpi_rank) std::cout << "Parameters : " << std::endl << parameters << std::endl;
      sim_type s(parameters, c);
#endif
      //run the simulation
      s.run(boost::bind(&stop_callback, time(0)+max_time));

      //on the master: collect MC results and store them in file, then postprocess
      if (global_mpi_rank==0){
        alps::results_type<hybridization>::type results = collect_results(s);
        std::string output_path = parameters["cthyb.BASEPATH"].as<std::string>()+std::string("/simulation/results");
        alps::save_results(results, parameters, "sim.h5", output_path); //"/simulation/results");
        master_final_tasks(results, parameters, "sim.h5");
#ifdef ALPS_HAVE_MPI
      } else{ //on any slave: send back results to master.
        collect_results(s);
      }
      c.barrier();
      MPI_Finalize();
#else
    }
#endif
    }
    catch(std::exception& exc){
      std::cerr<<exc.what()<<std::endl;
      return -1;
    }
    catch(...){
      std::cerr << "Fatal Error: Unknown Exception!\n";
      return -2;
    }
}


void master_final_tasks(const alps::results_type<hybridization>::type &results,
                        const alps::parameters_type<hybridization>::type &parameters,
                        const std::string &output_name){
  //do some post processing: collect Green functions and write
  //them into hdf5 files; calls compute vertex at the very end

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
