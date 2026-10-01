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

#include "hyb.hpp"
#include "hybevaluate.hpp"
#include <alps/utility/copyright.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#ifdef ALPS_HAVE_MPI
#include <alps/mcmpiadapter.hpp>
typedef alps::mcmpiadapter<hybridization> sim_type;
#else
typedef hybridization sim_type;
#endif

//stops the simulation if time > end_time or if signals received.
bool stop_callback(boost::posix_time::ptime const & end_time) {
  static alps::ngs::signal signal;
  return !signal.empty() || boost::posix_time::second_clock::local_time() > end_time;
}
void master_final_tasks(const alps::results_type<hybridization>::type &results, const alps::parameters_type<hybridization>::type &parms, const std::string &output_name);
int global_mpi_rank;

#ifdef BUILD_PYTHON_MODULE
#include "dict_to_params.hpp"
#include "scoped_signal_handlers.hpp"
namespace nb = nanobind;

void solve(nb::dict const & parms_){
  pyalps::scoped_signal_handlers signal_handlers;
  alps::parameters_type<hybridization>::type parms = pyalps::params_from_dict(parms_);
  std::string output_file = boost::lexical_cast<std::string>(parms["BASENAME"]|"results")+std::string(".out.h5");
#else
int main(int argc, char** argv){
  //read in command line options
  alps::mcoptions options(argc, argv);
  if (options.valid) {
    std::string output_file = options.output_file;

#ifdef ALPS_HAVE_MPI
    //boot up MPI environment
    boost::mpi::environment env(argc, argv);
#endif

    //create ALPS parameters from hdf5 parameter file
    alps::parameters_type<hybridization>::type parms(alps::hdf5::archive(options.input_file, alps::hdf5::archive::READ));
    try {
      if(options.time_limit!=0)
        throw std::invalid_argument("time limit is passed in the parameter file!");
      if(!parms.defined("MAX_TIME")) throw std::runtime_error("parameter MAX_TIME is not defined. How long do you want to run the code for? (in seconds)");
#endif

#ifndef ALPS_HAVE_MPI
      global_mpi_rank=0;
      sim_type s(parms,global_mpi_rank);
#else
      boost::mpi::communicator c;
      c.barrier();
      global_mpi_rank=c.rank();
      sim_type s(parms, c);
#endif
      if (global_mpi_rank==0) {
        alps::print_copyright(std::cout);
        std::cout << "****************************************************************"<<std::endl;
        std::cout << "* Recommended citation in scientific publications:             *"<<std::endl;
        std::cout << "* We used the ALPS [1] implementation [2] of the `segment'     *"<<std::endl;
        std::cout << "* CT-QMC solver [3,4,5].                                       *"<<std::endl;
        std::cout << "* [1] JSTAT (2011) P05001; [2] CPC 182, 1078 (2011);           *"<<std::endl;
        std::cout << "* [3] PRL 97, 076405 (2006); [4] RMP 83, 349 (2011);           *"<<std::endl;
        std::cout << "* [5] PRB 84, 075145 (2011).                                   *"<<std::endl;
        std::cout << "****************************************************************"<<std::endl;
      }
      //run the simulation
      s.run(boost::bind(&stop_callback, boost::posix_time::second_clock::local_time() + boost::posix_time::seconds((int)parms["MAX_TIME"])));

      //on the master: collect MC results and store them in file, then postprocess
      if (global_mpi_rank==0){
        alps::results_type<hybridization>::type results = collect_results(s);
        std::string output_path = boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/simulation/results";
        save_results(results, parms, output_file, output_path); //"/simulation/results");
        master_final_tasks(results, parms, output_file);
#ifdef ALPS_HAVE_MPI
      } else{ //on any slave: send back results to master.
        collect_results(s);
      }
      c.barrier();
#else
      }
#endif
#ifdef BUILD_PYTHON_MODULE
    return;
#else
    }
    catch(std::exception& exc){
      std::cerr<<exc.what()<<std::endl;
      return -1;
    }
    catch(...){
      std::cerr << "Fatal Error: Unknown Exception!\n";
      return -2;
    }
  }//options.valid
  return 0;
#endif

}


void master_final_tasks(const alps::results_type<hybridization>::type &results,
                        const alps::parameters_type<hybridization>::type &parms,
                        const std::string &output_name){
  //do some post processing: collect Green functions and write
  //them into hdf5 files; calls compute vertex at the very end

  alps::hdf5::archive solver_output(output_name, "a");

  evaluate_basics(results,parms,solver_output);
  evaluate_time(results,parms,solver_output);
  evaluate_freq(results,parms,solver_output);
  evaluate_legendre(results,parms,solver_output);
  evaluate_nnt(results,parms,solver_output);
  evaluate_nnw(results,parms,solver_output);
  evaluate_sector_statistics(results,parms,solver_output);
  evaluate_2p(results, parms, solver_output);
}

#ifdef BUILD_PYTHON_MODULE
NB_MODULE(cthyb, m) {
    m.def("solve", solve);
}
#endif


