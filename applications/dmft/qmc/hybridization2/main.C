/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2010 by Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Emanuel Gull <gull@phys.columbia.edu>,
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

#include "impurity.h"
#include "fouriertransform.h"
#include "alps/ngs/make_deprecated_parameters.hpp"
#include<numeric>
typedef alps::mcmpisim<hybridization> sim_type;

bool stop_callback(boost::posix_time::ptime const & end_time) {
  static alps::ngs::signal signal;
  return !signal.empty() || boost::posix_time::second_clock::local_time() > end_time;
}

int main(int argc, char** argv)
{
  alps::mcoptions options(argc, argv);
  if (options.valid) {
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator c;
    alps::parameters_type<hybridization>::type parms(alps::hdf5::archive(options.input_file, "r"));
    alps::mcmpisim<hybridization> s(parms, c);
    if(options.time_limit!=0)
      throw std::invalid_argument("time limit is passed in the parameter file!");
    if(!parms.defined("MAX_TIME")) throw std::runtime_error("parameter MAX_TIME is not defined. How long do you want to run the code for? (in seconds)");
    s.run(boost::bind(&stop_callback, boost::posix_time::second_clock::local_time() + boost::posix_time::seconds((int)parms["MAX_TIME"])));
    if (c.rank()==0){
      alps::results_type<hybridization>::type results = collect_results(s);
      save_results(results, parms, options.output_file, "/simulation/results");
      
      if(!(parms.defined("EXTERNAL_INPUT_DATA") && ((bool)(parms["EXTERNAL_INPUT_DATA"])==true))){//if external input data is used, we skip this step
      //compute the output Green's function and Fourier transform it, store in the right path
      int N = static_cast < int >(parms["N"]);
      int FLAVORS = static_cast < int >(parms["FLAVORS"]);
        itime_green_function_t green_itime(N + 1, 1, FLAVORS);
        matsubara_green_function_t green_matsubara(N, 1, FLAVORS);
        boost::shared_ptr < FourierTransformer > fourier_ptr;
        std::vector<double> G       = results["Greens"].mean<std::vector<double> >();
        std::vector<double> G_error = results["Greens"].error<std::vector<double> >();
        std::vector<double> n       = results["n"].mean<std::vector<double> >();
        std::vector<double> n_error = results["n"].mean<std::vector<double> >();
        {
          //write matrix sizes and overlap to file
          std::ofstream matrix_size("matrix_size", std::ios::app);
          std::vector<double> order_mean =results["matrix_size"].mean <std::vector<double> >();
          std::vector<double> order_error=results["matrix_size"].error<std::vector<double> >();
          for(unsigned int i=0;i<order_mean.size();++i){
            matrix_size << order_mean[i] << " "<<order_error[i]<<" ";
          }
          double total_size=std::accumulate(order_mean.begin(), order_mean.end(), 0.);
          double total_error=std::accumulate(order_error.begin(), order_error.end(),0.);
          matrix_size<<total_size<<" "<<std::sqrt(order_error.size())*total_error<<std::endl; //independence assumption. Not sure if this is a good assumption.
          /*std::vector<double> overlap_mean =  results["overlap"].mean <std::vector<double> >() ;
          std::vector<double> overlap_error = results["overlap"].error<std::vector<double> >() ;
          std::ofstream overlap_file("overlap", std::ios::app);
          for (unsigned i = 0; i < overlap_mean.size(); ++i) {
            overlap_file << overlap_mean[i]<<" "<<overlap_error[i]<< "\t";
          }
          overlap_file << std::endl;*/
        }
        std::cout<<"doing green itime conversion"<<std::endl;

        for (int f = 0; f < FLAVORS; f++) {
          for (int i = 0; i < N + 1; i++) {
            green_itime(i, f) = -G[f * (N + 1) + i];
            green_itime.error(i, f) = G_error[f * (N + 1) + i];
          }
          green_itime(N, f) = -n[f];
          green_itime(0, f) = -1 + n[f];
          green_itime.error(N, f) = n_error[f];
          green_itime.error(0, f) = n_error[f];
        }
        {
          std::cout<<"converting ALPS parameters"<<std::endl;
          alps::Parameters p(alps::make_deprecated_parameters(parms));
          FourierTransformer::generate_transformer_U(p, fourier_ptr, n); //still takes old alps parameter class.
          fourier_ptr->forward_ft(green_itime, green_matsubara);      
        }
        {
          alps::hdf5::archive ar(options.output_file, "a");
          green_matsubara.write_hdf5(ar, "/G_omega");
          green_itime.write_hdf5(ar, "/G_tau");
        }
      }
    } else{//c.rank!=0
      collect_results(s);
    }
  }
}
