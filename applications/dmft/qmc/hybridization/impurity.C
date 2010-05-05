 /*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Philipp Werner <werner@itp.phys.ethz.ch>
 *                              Emanuel Gull <gull@phys.columbia.edu>
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

#include "types.h"
#include <alps/lattice.h>
#include "impurity.h"
#include "moves.h"
#include "fouriertransform.h"
#include "xml.h"
#include <alps/alea.h>
//#include "matrix_builder.h"
using namespace std;
using namespace alps;


HybridizationRun::HybridizationRun(const alps::ProcessList & where, const alps::Parameters & p, int node)
:alps::scheduler::MCRun(where, p, node), sweeps(0),
thermalization_sweeps(static_cast < int >(parms["THERMALIZATION"])),
total_sweeps(static_cast < int >(parms["SWEEPS"])), start_time(time(NULL)),
max_time((int) (parms.value_or_default("MAX_TIME", 86400))),
mu(static_cast < double >(parms["MU"])),
BETA(static_cast < double >(parms["BETA"])),
u(parms),
full_line(static_cast < int >(parms["FLAVORS"]), 0),
sign(static_cast < int >(parms["FLAVORS"])), 
G_meas(static_cast < int >(parms["FLAVORS"]) * (static_cast < int >(parms["N"]) + 1))
{
  int N = static_cast < int >(parms["N"]);
  int FLAVORS = static_cast < int >(parms["FLAVORS"]);
  int Np1 = N + 1;
  // define interactions between flavors

  // for semi-circular density of states we can directly compute F_{\sigma}(\tau) = t^2*G_{-\sigma}(-\tau)
  std::vector < std::vector < double > >G;

  G.resize(FLAVORS);
  F.resize(FLAVORS);
  for (int i = 0; i < FLAVORS; i++) {
    G[i].resize(Np1);
    F[i].resize(Np1);
  }
  itime_green_function_t f_itime(Np1, 1, FLAVORS);
  double mu_shift = u.mu_shift();
  mu = mu + mu_shift;
  //std::cout << "mu_shift to half filling is: " << mu_shift << std::endl;
  //general case: SC loop in omega, we have to convert G0(iomega) into F(tau).
  if (parms.defined("OMEGA_LOOP") && (bool)(parms["OMEGA_LOOP"])==true) {
    matsubara_green_function_t bare_green_matsubara(N, 1, FLAVORS);
    matsubara_green_function_t f_twiddle_matsubara(N, 1, FLAVORS);
    itime_green_function_t f_twiddle_itime(Np1, 1, FLAVORS);
    std::cout << "U is: " << u << std::endl;
    //find the second moment of the band structure
    double epssqav ;
    if (parms.defined("DOSFILE")) {
      if (!parms.defined("EPSSQAV")) {
        throw std::logic_error("error: you specify a DOS file, please also specify the second moment of the band structure EPSSQAV!");
      } else {
        epssqav = parms["EPSSQAV"];
      }
    }else{
      t=parms["t"]; //this is essentially an energy unit.
      epssqav = t * t;
    }
    std::istringstream in_omega(parms["G0(omega)"]);
    read_freq(in_omega, bare_green_matsubara);
    //std::cout << "absolute mu is: " << mu << std::endl;
    FFunctionFourierTransformer Fourier(BETA , 0, epssqav , FLAVORS, 1);
    for (int f = 0; f < FLAVORS; ++f) {
      for (int i = 0; i < N; ++i) {
        f_twiddle_matsubara(i, f) = -1. / bare_green_matsubara(i, f) + (std::complex < double >(mu - mu_shift, (2. * i + 1) * M_PI / BETA));
      }
    }
    Fourier.backward_ft(f_twiddle_itime, f_twiddle_matsubara);
    for (int f = 0; f < FLAVORS; ++f) {
      for (int i = 0; i < Np1; ++i) {
        f_itime(i, f) = -f_twiddle_itime(N - i, f);
      }
    }

    //std::cout << "starting WERNER simulation with F function:\n" << f_itime<< std::endl;
  } else {                      //no omega loop. Use Bethe Lattice SC as described in the paper.
    std::istringstream in_tau(parms["G0"]);     //this is not G0 but G(tau). See solver.solve in F_selfconsistency_loop()
    itime_green_function_t green_itime(Np1, 1, FLAVORS);
    read_itime(in_tau, green_itime);
    for (int f = 0; f < FLAVORS; ++f) {
      int orbital=f/2;
      std::stringstream tname; tname<<"t"<<orbital;
      if(parms.defined(tname.str())) {
        t=parms[tname.str()]; 
        std::cout<<"orbital: "<<f/2<<" flavor: "<<f<<" using t: "<<t<<std::endl;
      }
      else 
        t=parms["t"];
      for (int i = 0; i < Np1; ++i) {
        f_itime(i, f) = -t * t * green_itime(N - i, f);  //this is the self consistency loop, for Bethe lattice!
      }
    }
    //std::cout<<green_itime<<std::endl;
    std::cout << "U is: " << u << std::endl;
  }
  //std::cout<<f_itime<<std::endl;
  for (int i = 0; i < Np1; i++) {
    for (int j = 0; j < FLAVORS; j++) {
      F[j][i] = f_itime(i, j);
    }
  }

  segments.resize(FLAVORS);
  M.resize(FLAVORS, blas::matrix(0));
  // initialize list of segments
  for (int i = 0; i < FLAVORS; i++) {
    segments[i].clear();
    M[i].clear();
  }

  // create measurement objects
  measurements << RealVectorObservable("n");
  measurements << RealVectorObservable("order");
  measurements << vec_obs_t("Greens");
  measurements << vec_obs_t("Greens_4point");
  measurements << vec_obs_t("Greens_4point_ud");
  measurements << RealObservable("sign");
  measurements << RealObservable("MatrixSize");

  measurements << RealVectorObservable("overlap");

}

bool HybridizationRun::change_parameter(const std::string & name, const alps::StringValue & value)
{
  if (name == "SWEEPS")
    total_sweeps = static_cast < alps::uint32_t > (value);
  else if (name == "THERMALIZATION" && !is_thermalized())
    thermalization_sweeps = static_cast < alps::uint32_t > (value);
  else
    return false;               // cannot do it
  return true;                  // could do it
}

bool HybridizationRun::is_thermalized() const
{
  return (sweeps >= thermalization_sweeps);
}

double HybridizationRun::work_done() const
{
  if (time(NULL) - start_time > max_time)
    return 1.;
  return (is_thermalized()? (sweeps - thermalization_sweeps) / double (total_sweeps) : 0.);

}

void HybridizationRun::dostep()
{

  // increment sweep count
  sweeps++;

  static int N(static_cast < int >(parms["N"]));
  static int FLAVORS(static_cast < int >(parms["FLAVORS"]));
  static int N_order(static_cast < int >(parms["N_ORDER"]));
  static int N_meas(static_cast < int >(parms["N_MEAS"]));
  static int N_shift(static_cast < int >(parms["N_SHIFT"]));
  
  int overlap = true;           //static_cast<int>(parms["OVERLAP"]);
  std::valarray < double >overlap_meas(FLAVORS * (FLAVORS - 1) / 2);

  times full_segment(0, BETA);

  segment_container_t::const_iterator it1, it2;


  std::valarray < double >order_meas(N_order * FLAVORS);
  std::valarray < double >n_meas(FLAVORS);


  double sign_meas = 0, s = 1;
  G_meas = 0;
  n_meas = 0;

  for (int i = 0; i < N_meas; i++) {
    for (int j = 0; j < FLAVORS; j++) {
      if (segments[j].size() == 0) {
        // insert or remove full line
        insert_remove_full_line(random_01, mu, u, BETA, full_line[j], segments, full_line, j);
      }

      insert_remove_antisegment(random_01, BETA * random_01(), BETA, mu, u,
                                F[j], full_line[j], segments[j], M[j], sign[j], segments, full_line, j);

      if (!full_line[j]) {
        // local update
        insert_remove_segment(random_01, BETA * random_01(), BETA, mu, u, F[j], segments[j], M[j], sign[j], segments, full_line, j);

        // shift segments
        for (int k = 0; k < N_shift; k++)
          shift_segment(random_01, segments[j], BETA, mu, u, F[j], M[j], sign[j], segments, full_line, j);

      }

      if ((int) (segments[j].size()) < N_order)
        order_meas[j * FLAVORS + segments[j].size()] += 1;

      if (segments[j].size() > 0) {
        double N_div_beta = N / BETA;
        double *start_times=new double[segments[j].size()];
        double *end_times=new double[segments[j].size()];
        int i = 0;
        for (it1 = segments[j].begin(); it1 != segments[j].end(); ++it1) {
          start_times[i] = it1->t_start();
          end_times[i] = it1->t_end();
          ++i;
        }
        for (int i = 0; i < (int) M[j].size1(); i++) {
          for (int k = 0; k < (int) M[j].size1(); k++) {
            if (M[j] (k, i) != 0.) {
              double argument = end_times[i] - start_times[k];
              double bubble_sign = 1.;
              if (argument < 0) {
                bubble_sign = -1.;
                argument += BETA;
              }
              int index = (int) (argument * N_div_beta + 0.5);
              double beta_G = M[j] (k, i) * bubble_sign; //note the missing 1/beta that is added at the end
              G_meas[j * (N + 1) + index] += beta_G;
            }                   //if M != 0
          }                     //for k
        }                       //for i
        delete[] start_times;
        delete[] end_times;
      }                         // if there are segments.

      s *= sign[j];
      n_meas[j] += compute_overlap(full_segment, segments[j], full_line[j], BETA) / BETA;

    }

    sign_meas += s;

  }
  if (parms.defined("MEASURE_FOURPOINT"))
    measure_fourpoint();

  if (overlap) {
    int orbital_index = 0;
    for (int j = 1; j < FLAVORS; ++j) {        //overlap between orbital j
      for (int k = 0; k < j; ++k) {    //and orbital k
        if (segments[j].size() > 0) {
          for (it1 = segments[j].begin(); it1 != segments[j].end(); it1++) {
            overlap_meas[orbital_index] += compute_overlap(*it1, segments[k], full_line[k], BETA) / BETA;
          }
        } else if (full_line[j]) {
          overlap_meas[orbital_index] += compute_overlap(full_segment, segments[k], full_line[k], BETA) / BETA;
        }
        orbital_index++;
      }
    }
  }

  order_meas /= N_meas;
  measurements.get < RealVectorObservable > ("order") << order_meas;

  G_meas *= (1. * N) / (N_meas*BETA*BETA); //one beta from the binning x-axis size, one beta from the 1/beta \sum M...
  measurements.get < vec_obs_t > ("Greens") << G_meas;  //*sign;
  sign_meas /= N_meas;
  measurements.get < RealObservable > ("sign") << sign_meas;

  n_meas /= N_meas;
  measurements.get < RealVectorObservable > ("n") << n_meas;

  if (overlap)
    measurements.get < RealVectorObservable > ("overlap") << overlap_meas;

  measurements.get < RealObservable > ("MatrixSize") << M[0].size1();
}


std::pair < matsubara_green_function_t, itime_green_function_t > HybridizationSimFrequency::get_result() 
{

  int N = static_cast < int >(parms["N"]);
  int FLAVORS = static_cast < int >(parms["FLAVORS"]);

  itime_green_function_t multiple_G(N + 1, 1, FLAVORS);
  matsubara_green_function_t green_matsubara(N, 1, FLAVORS);
  boost::shared_ptr < FourierTransformer > fourier_ptr;
  alps::RealVectorObsevaluator G = get_measurements()["Greens"];
  alps::RealVectorObsevaluator n = get_measurements()["n"];
  //std::cout << "total measurements done: " << G.count() << std::endl;
  if(parms.defined("CHECKPOINT")){
    alps::Parameters *p=const_cast<alps::Parameters*>(&parms);
    p->erase(std::string("G0(omega)"));
    alps::Parameters::iterator it=p->begin();
    while(it!=p->end()){
      if(strncmp(it->key().c_str(),"G0(omega)",9)==0){
        it->value()="...ignored, nested...";
      }
      it++;
    }
    //std::cout<<"end of sim, checkpointing"<<std::endl;
    std::string fns=parms["CHECKPOINT"];
    fns+=".xml";
    boost::filesystem::path fn(fns);
    checkpoint(fn);
  }
  {
    //write matrix sizes to file
    std::ofstream matrix_size("matrix_size", std::ios::app);
    matrix_size << RealObsevaluator(get_measurements()["MatrixSize"]).mean() << std::endl;
    RealVectorObsevaluator overlap = get_measurements()["overlap"];
    std::ofstream overlap_file("overlap", std::ios::app);
    for (unsigned i = 0; i < overlap.mean().size(); ++i) {
      overlap_file << overlap.mean()[i] << "\t";
    }
    overlap_file << std::endl;
  }
  if(parms.defined("CHECKPOINT")){
    alps::Parameters *p=const_cast<alps::Parameters*>(&parms);
    p->erase(std::string("G0(omega)"));
    alps::Parameters::iterator it=p->begin();
    while(it!=p->end()){
      if(strncmp(it->key().c_str(),"G0(omega)",9)==0){
        it->value()="...ignored, nested...";
      }
      it++;
    }
    //std::cout<<"end of sim, checkpointing"<<std::endl;
    std::string fns=parms["CHECKPOINT"];
    fns+=".xml";
    boost::filesystem::path fn(fns);
    checkpoint(fn);
  }

  for (int f = 0; f < FLAVORS; f++) {

    for (int i = 0; i < N + 1; i++) {
      multiple_G(i, f) = -G.mean()[f * (N + 1) + i];
      multiple_G.error(i, f) = G.error()[f * (N + 1) + i];
    }
    multiple_G(N, f) = -n.mean()[f];
    multiple_G(0, f) = -1 + n.mean()[f];
    multiple_G.error(N, f) = n.error()[f];
    multiple_G.error(0, f) = n.error()[f];
  }
  std::vector<double> densities;
  densities.resize(FLAVORS);   
  for (int f = 0; f < FLAVORS; ++f) {
    densities[f] = n.mean()[f];
  }
  FourierTransformer::generate_transformer_U(parms, fourier_ptr, densities);
  fourier_ptr->forward_ft(multiple_G, green_matsubara);
  if (parms.defined("MEASURE_FOURPOINT")) {
    write_fourpoint();
  }
  return std::make_pair(green_matsubara, multiple_G);

}

itime_green_function_t HybridizationSimItime::get_result() const
{

  int N = static_cast < int >(parms["N"]);
  int FLAVORS = static_cast < int >(parms["FLAVORS"]);

  itime_green_function_t multiple_G(N + 1, 1, FLAVORS);
  alps::RealVectorObsevaluator G = get_measurements()["Greens"];
  alps::RealVectorObsevaluator n = get_measurements()["n"];

  {
    //write matrix sizes to file
    std::ofstream matrix_size("matrix_size", std::ios::app);
    matrix_size << RealObsevaluator(get_measurements()["MatrixSize"]).mean() << std::endl;
    RealVectorObsevaluator overlap = get_measurements()["overlap"];
    std::ofstream overlap_file("overlap", std::ios::app);
    overlap_file << setprecision(20);
    for (unsigned i = 0; i < overlap.mean().size(); ++i) {
      overlap_file << overlap.mean()[i] << "\t";
    }
    overlap_file << std::endl;
  }
  for (int f = 0; f < FLAVORS; f++) {

    for (int i = 0; i < N + 1; i++) {
      multiple_G(i, f) = -G.mean()[f * (N + 1) + i];
      multiple_G.error(i, f) = G.error()[f * (N + 1) + i];
    }
    multiple_G(N, f) = -n.mean()[f];
    multiple_G(0, f) = -1 + n.mean()[f];
    multiple_G.error(N, f) = n.error()[f];
    multiple_G.error(0, f) = n.error()[f];
  }
  return multiple_G;


}
void print_segments(const std::vector<segment_container_t> &segments){
  for (std::vector<segment_container_t>::const_iterator it=segments.begin();it!=segments.end();++it){
    if(it->size() ==0 ){std::cout<<"empty segment"<<std::endl;}
    else{
      for (segment_container_t::const_iterator it2=it->begin();it2!=it->end();++it2){
        std::cout<<it2->t_start()<<"->"<<it2->t_end()<<"; ";
      }
      std::cout<<std::endl;
    }
  }
}

