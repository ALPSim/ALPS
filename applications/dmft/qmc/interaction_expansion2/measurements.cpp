/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
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
#include <complex>
#include <alps/alea.h>
#include "alps/ngs/make_deprecated_parameters.hpp"

void evaluate_selfenergy_measurement_matsubara(const alps::results_type<HubbardInteractionExpansion>::type &results, 
                                                                        matsubara_green_function_t &green_matsubara_measured,
                                                                        const matsubara_green_function_t &bare_green_matsubara, 
                                                                        std::vector<double>& densities,
                                                                        const double &beta, std::size_t n_site, 
                                                                        std::size_t n_flavors, std::size_t n_matsubara);
void evaluate_selfenergy_measurement_itime_rs(const alps::results_type<HubbardInteractionExpansion>::type &results, 
                                                                       itime_green_function_t &green_result,
                                                                       const itime_green_function_t &green0, 
                                                                       const double &beta, const int n_site, 
                                                                       const int n_flavors, const int n_tau, const int n_self);



void compute_greens_functions(const alps::results_type<HubbardInteractionExpansion>::type &results, const alps::parameters_type<HubbardInteractionExpansion>::type& parms, const std::string &output_file)
{
  std::cout<<"getting result!"<<std::endl;
  unsigned int n_matsubara = parms["NMATSUBARA"]|parms["N_MATSUBARA"];
  unsigned int n_matsubara_measurements=parms["NMATSUBARA_MEASUREMENTS"] | (int)n_matsubara;
  unsigned int n_tau=parms["N"]|parms["N_TAU"];
  unsigned int n_self=parms["NSELF"] | (int)(10*n_tau);
  spin_t n_flavors(parms["FLAVORS"] | (parms["N_ORBITALS"]| 2));
  unsigned int n_site(parms["SITES"] | 1);
  double beta(parms["BETA"]);
  itime_green_function_t green_itime_measured(n_tau+1, n_site, n_flavors);
  matsubara_green_function_t green_matsubara_measured(n_matsubara, n_site, n_flavors);
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  boost::shared_ptr<FourierTransformer> fourier_ptr_g0;
  FourierTransformer::generate_transformer(alps::make_deprecated_parameters(parms), fourier_ptr_g0);
  //find whether our data is in imaginary time or frequency:
  bool measure_in_matsubara=true;
  if(parms["HISTOGRAM_MEASUREMENT"] | false) 
    measure_in_matsubara=false;
  std::vector<double> mean_order=results["PertOrder"].mean<std::vector<double> >();
  
  std::cout<<"average matrix size was: "<<std::endl;
  std::ofstream matrix_size("matrix_size", std::ios::app);
  for(unsigned int i=0;i<n_flavors;++i){
    std::cout<<mean_order[i]<<"\t";
    matrix_size<<mean_order[i]<<"\t";
  }
  std::cout<<std::endl;
  matrix_size<<std::endl;
  std::cout<<"average sign was: "<<results["Sign"].mean<double>()<<" error: "<<results["Sign"].error<double>()<<std::endl;
  //single particle Green function measurements
  matsubara_green_function_t bare_green_matsubara(n_matsubara, n_site, n_flavors);
  std::vector<double> densities(n_flavors);
  {
    alps::hdf5::archive ar(parms["INFILE"],"r");
    if(parms.defined("DMFT_FRAMEWORK") && static_cast<bool>(parms["DMFT_FRAMEWORK"])){
      //read in as green_function
      bare_green_matsubara.read_hdf5(ar,"/G0");
      
    } else { //plain hdf5
      std::vector<std::complex<double> > tmp(n_matsubara);
      for(std::size_t j=0; j<n_flavors; j++){
        std::stringstream path; path<<"/G0_"<<j;
        ar>>alps::make_pvp(path.str(),tmp);
        for(std::size_t i=0; i<n_matsubara; i++)
          bare_green_matsubara(i,0,0,j)=tmp[i];
      }
      tmp.clear();
    }
  }
  if(measure_in_matsubara) {
    evaluate_selfenergy_measurement_matsubara(results, green_matsubara_measured, 
                                              bare_green_matsubara, densities, 
                                              beta, n_site, n_flavors, n_matsubara_measurements);
  } 
  else {
    itime_green_function_t bare_green_itime(n_tau+1, n_site, n_flavors);
    fourier_ptr_g0->backward_ft(bare_green_itime, bare_green_matsubara);
    evaluate_selfenergy_measurement_itime_rs(results, green_itime_measured, bare_green_itime, 
                                             beta, n_site, n_flavors, n_tau, n_self);
  }
  //Fourier transformations
  if (!measure_in_matsubara) {
    for (unsigned int z=0; z<n_flavors; ++z) {
      densities[z] = 0;
      for (unsigned int i=0; i<n_site; ++i)
        densities[z] -= green_itime_measured(n_tau,i,i,z);
      densities[z] /= n_site;
    }
  }
  FourierTransformer::generate_transformer_U(alps::make_deprecated_parameters(parms), fourier_ptr, densities);
  if (measure_in_matsubara) {
    fourier_ptr->append_tail(green_matsubara_measured, bare_green_matsubara, n_matsubara_measurements);
    fourier_ptr->backward_ft(green_itime_measured, green_matsubara_measured);
  }
  else 
    fourier_ptr->forward_ft(green_itime_measured, green_matsubara_measured);
  {
    alps::hdf5::archive ar(output_file, "a");
    green_matsubara_measured.write_hdf5(ar, "/G_omega");
    green_itime_measured.write_hdf5(ar, "/G_tau");
  }
} 
