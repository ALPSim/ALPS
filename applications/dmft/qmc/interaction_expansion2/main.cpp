/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2010 by Emanuel Gull <gull@phys.columbia.edu>,
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

#include "interaction_expansion.hpp"
#include "fouriertransform.h"

#include <alps/ngs/scheduler/mpi_adapter.hpp>

bool stop_callback(boost::posix_time::ptime const & end_time) {
  static alps::ngs::signal signal;
  return !signal.empty() || boost::posix_time::second_clock::local_time() > end_time;
}
void compute_greens_functions(const alps::results_type<HubbardInteractionExpansion>::type &results, const alps::parameters_type<HubbardInteractionExpansion>::type& parms, const std::string &output_file);

int main(int argc, char** argv)
{
  alps::mcoptions options(argc, argv);
  if (options.valid) {
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator c;
    alps::parameters_type<HubbardInteractionExpansion>::type parms(alps::hdf5::archive(options.input_file, "r"));
    alps::mpi_adapter<HubbardInteractionExpansion> s(parms, c);
    if(options.time_limit!=0)
      throw std::invalid_argument("time limit is passed in the parameter file!");
    if(!parms.defined("MAX_TIME")) throw std::runtime_error("parameter MAX_TIME is not defined. How long do you want to run the code for? (in seconds)");
    s.run(boost::bind(&stop_callback, boost::posix_time::second_clock::local_time() + boost::posix_time::seconds((int)parms["MAX_TIME"])));
    if (c.rank()==0){
      alps::results_type<HubbardInteractionExpansion>::type results = collect_results(s);
      save_results(results, parms, options.output_file, "/simulation/results");
      
      //compute the output Green's function and Fourier transform it, store in the right path
      compute_greens_functions(results, parms, options.output_file);
    } else{
      collect_results(s);
    }
  }
}
