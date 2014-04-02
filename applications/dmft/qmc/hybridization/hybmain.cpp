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

#include "hyb.hpp"
#include "hybevaluate.hpp"
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
//compile it as a python module (requires boost::python library)
using namespace boost::python;

void solve(boost::python::dict parms_){
  alps::parameters_type<hybridization>::type parms(parms_);
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
BOOST_PYTHON_MODULE(cthyb)
{
    def("solve",solve);//define python-callable run method
};
#endif




