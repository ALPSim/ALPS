/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2010 by Emanuel Gull <gull@phys.columbia.edu>,
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#include "interaction_expansion.hpp"
#include "fouriertransform.h"
#include <alps/utility/copyright.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#ifdef ALPS_HAVE_MPI
#include <alps/mcmpiadapter.hpp>
typedef alps::mcmpiadapter<HubbardInteractionExpansion> sim_type;
#else
typedef HubbardInteractionExpansion sim_type;
#endif

bool stop_callback(boost::posix_time::ptime const & end_time) {
//stops the simulation if time > end_time or if signals received.
  static alps::ngs::signal signal;
  return !signal.empty() || boost::posix_time::second_clock::local_time() > end_time;
}
void compute_greens_functions(const alps::results_type<HubbardInteractionExpansion>::type &results, const alps::parameters_type<HubbardInteractionExpansion>::type& parms, const std::string &output_file);
int global_mpi_rank;

#ifdef BUILD_PYTHON_MODULE
#include "dict_to_params.hpp"
#include "scoped_signal_handlers.hpp"
namespace nb = nanobind;

void solve(nb::dict const & parms_){
    pyalps::scoped_signal_handlers signal_handlers;
    alps::parameters_type<HubbardInteractionExpansion>::type parms = pyalps::params_from_dict(parms_);
    std::string output_file = boost::lexical_cast<std::string>(parms["BASENAME"]|"results")+std::string(".out.h5");
#else
int main(int argc, char** argv)
{
  alps::mcoptions options(argc, argv);
  if (options.valid) {
    std::string output_file = options.output_file;

#ifdef ALPS_HAVE_MPI
    //boot up MPI environment
    boost::mpi::environment env(argc, argv);
#endif

    //create ALPS parameters from hdf5 parameter file
    alps::parameters_type<HubbardInteractionExpansion>::type parms(alps::hdf5::archive(options.input_file, alps::hdf5::archive::READ));
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
        std::cout << "* We used the ALPS [1] implementation [2] of the CT-INT        *"<<std::endl;
        std::cout << "* interaction expansion CT-QMC [3,4] solver.                   *"<<std::endl;
        std::cout << "* [1] JSTAT (2011) P05001; [2] CPC 182, 1078 (2011);           *"<<std::endl;
        std::cout << "* [3] PRB 72, 035122 (2005); [4] RMP 83, 349 (2011).           *"<<std::endl;
        std::cout << "****************************************************************"<<std::endl;
      }
      //run the simulation
      s.run(boost::bind(&stop_callback, boost::posix_time::second_clock::local_time() + boost::posix_time::seconds((int)parms["MAX_TIME"])));

      //on the master: collect MC results and store them in file, then postprocess
      if (global_mpi_rank==0){
        alps::results_type<HubbardInteractionExpansion>::type results = collect_results(s);
        save_results(results, parms, output_file, "/simulation/results");
        //compute the output Green's function and Fourier transform it, store in the right path
        compute_greens_functions(results, parms, output_file);
#ifdef ALPS_HAVE_MPI
      } else{
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

#ifdef BUILD_PYTHON_MODULE
    NB_MODULE(ctint, m) {
        m.def("solve", solve);
    }
#endif
    
