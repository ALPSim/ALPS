/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2010 by Emanuel Gull <gull@phys.columbia.edu>,
 *
 *
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#include "interaction_expansion.hpp"
#include "fouriertransform.h"
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
//compile it as a python module (requires boost::python library)
using namespace boost::python;

void solve(boost::python::dict parms_){
    alps::parameters_type<HubbardInteractionExpansion>::type parms(parms_);
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
    BOOST_PYTHON_MODULE(ctint)
    {
        def("solve",solve);//define python-callable run method
    };
#endif
    
