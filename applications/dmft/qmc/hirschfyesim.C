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

/* $Id: hirschfyesim.C 340 2009-01-28 16:20:03Z fuchs $ */

/// @file hirschfyesim.h
/// @brief the actual Hirsch-Fye simulation

#include "hirschfyesim.h"
#include "xml.h"
#include "fouriertransform.h"
#include <alps/alea.h>
#include <cmath>
#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>

#include <alps/config.h> // needed to set up correct bindings
#include <boost/filesystem/operations.hpp>
#include <boost/numeric/bindings/lapack/driver/gesv.hpp>
using namespace std;
using namespace alps;

/// compute largest element of a matrix
double norm_max(dense_matrix & m) {
  
  double max(0);
  for (dense_matrix::const_iterator1 it = m.begin1(); it != m.end1(); ++it)
    if (std::fabs((*it)) > max)
      max = std::fabs((*it));
  
  return max;
  
}

// TODO remove once dense_matrix is removed
alps::numeric::matrix<double> alps_matrix_from_dense_matrix(dense_matrix const& m) {
    alps::numeric::matrix<double> r(m.size1(),m.size2());
    // Matrices are column major
    for(std::size_t j=0; j < m.size2(); ++j)
        for(std::size_t i=0; i < m.size1(); ++i)
            r(i,j) = m(i,j);
    return r;
}

/// compute Green's function for given spin configuration from Green0
void update_from_zero(dense_matrix & Green, dense_matrix & Green0, std::vector<int> & spins, double l) {
  
  alps::numeric::matrix<double> M(alps_matrix_from_dense_matrix(Green)); // for sign problem cases: figure out determinant...
  int N(spins.size());
  std::vector<double> e_ls(N);
  for(int i=0;i<N;++i){
    e_ls[i]=exp(l*spins[i]);
  }
  // calculate a = 1 + (1-g)*(exp(v')-1) 
  dense_matrix a(N, N);
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      a(i,j) = -Green0(i,j)*(e_ls[j]-1);
    }
    a(i,i) += e_ls[i]; 
  }
  // solve a*Green=Green0
  // gesv solves A*X=B, input for B is Green0, output (=solution X) is Green   
  Green = Green0; 
  boost::numeric::ublas::vector<fortran_int_t> ipivot(N);
  boost::numeric::bindings::lapack::gesv(a,ipivot, Green);
}


HirschFyeRun::HirschFyeRun(const alps::ProcessList& where,const alps::Parameters& p,int node)
: alps::scheduler::MCRun(where,p,node),  
sweeps(0),
thermalization_sweeps(static_cast<int>(parms["THERMALIZATION"])),
total_sweeps(static_cast<int>(parms["SWEEPS"])),
beta(static_cast<double>(parms["BETA"])),
N(static_cast<int>(parms["N"])),	
n_site(static_cast<int>(parms.value_or_default("SITES", 1))),	
u(static_cast<double>(parms["U"])),
N_check(1000), // value is automatically adjusted during the simulation
check_counter(0),
fp_interval(static_cast<int>(parms.value_or_default("FOURPOINT_INTERVAL", 1000))),
measure_fourpoint(static_cast<bool>(parms.value_or_default("MEASURE_FOURPOINT_FUNCTION", false))),
tolerance(static_cast<double>(parms["TOLERANCE"])),	
Green0_up(N*n_site, N*n_site),
Green_up(N*n_site, N*n_site),
Green0_down(N*n_site, N*n_site),
Green_down(N*n_site, N*n_site),			
spins(static_cast<int>(N*n_site)),
max_time(static_cast<int>(parms.value_or_default("MAX_TIME", 86400))),
start_time(time(NULL)),
sign(1),
bare_green_tau(N+1, n_site, 2),
green_tau(N+1, n_site, 2)
{
  std::cout<<"starting Hirsch Fye Simulation."<<std::endl;
  std::cout<<"The Hirsch Fye legacy code is published exclusively for pedagogical purposes."<<std::endl;
  
  unsigned int n_matsubara=parms["NMATSUBARA"];
  matsubara_green_function_t bare_green_matsubara(n_matsubara, n_site, 2);
  itime_green_function_t bgf_sc_convention_k(N+1,n_site, 2); //bare green function in k_space
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  FourierTransformer::generate_transformer(parms, fourier_ptr);
  
  std::istringstream in_omega(parms["G0(omega)"]);
  read_freq(in_omega, bare_green_matsubara);
  
  if (n_site > 1) {
    throw std::logic_error("please use the cluster solver version of this solver for clusters!");
  } 
  fourier_ptr->backward_ft(bare_green_tau, bare_green_matsubara);
  
  for(int i=0;i<2;++i){  //change convention to positive GF
    for(int s1=0;s1<n_site;++s1){
      for(int s2=0;s2<n_site;++s2){
        for(int t=0;t<=N;++t){
          bare_green_tau(t,s1, s2, i)*=-1;
        }
      }
    }
  }
  
  // initialize lambda
  double tmp = exp(beta/N*u/2);
  lambda = log(tmp + sqrt(tmp*tmp-1));
  
  // initialize bath Green's function matrices
  green_matrix_from_vector(bare_green_tau, Green0_up, Green0_down);
  // choose random spin configuration
  for (unsigned int i=0; i<spins.size(); i++)
    spins[i] = (random_01()<0.5 ? 1 : -1);
  
  
  // initialize Green's function matrices
  update_from_zero(Green_up, Green0_up, spins, lambda);
  update_from_zero(Green_down, Green0_down, spins, -lambda);
  
  
  // create measurement objects
  measurements << RealObservable("Sign");
  measurements <<SignedObservable<alps::RealVectorObservable>("G_meas_up");
  measurements <<SignedObservable<alps::RealVectorObservable>("G_meas_down");  
#ifdef FOURPOINT
  measurements << SimpleRealVectorObservable("G_fourpoint_uu");  
  measurements << SimpleRealVectorObservable("G_fourpoint_dd");  
  measurements << SimpleRealVectorObservable("G_fourpoint_ud");  
  measurements << SimpleRealVectorObservable("G_fourpoint_du");  
#endif //FOURPOINT
  
  measurements.reset(true);
  
}


// in principle, need only save spins (and do update_from_zero after loading)
void HirschFyeRun::load(alps::IDump& dump) {
  throw std::logic_error("implement load first!");
}


void HirschFyeRun::save(alps::ODump& /*dump*/) const
{
}


bool HirschFyeRun::is_thermalized() const
{
  return (sweeps >= thermalization_sweeps);
}

double HirschFyeRun::work_done() const
{
  if(time(NULL)-start_time> max_time){
    std::cout<<"we ran out of time!"<<std::endl<<"sweeps: "<<sweeps-thermalization_sweeps<<std::endl;
    return 1;
  }
  
  return (is_thermalized() ? (sweeps-thermalization_sweeps)/double(total_sweeps) : 0.);
}


void HirschFyeRun::dostep()
{
  
  // increment sweep count
  sweeps++;
  check_counter++;
  // do local update
  for (int i=0; i<N; i++){
    update_single_spin(random_01, Green_up, Green_down, spins, lambda, sign);
  }
  green_vector_from_matrix(green_tau, Green_up, Green_down);
  
  // measure Greens function
  std::valarray<double> G_up(N*n_site*n_site);
  std::valarray<double> G_down(N*n_site*n_site);
  for(int s1=0;s1<n_site;++s1){
    for(int s2=0;s2<n_site;++s2){
      for(int t=0;t<N;++t){
        G_up  [t+N*s1+N*n_site*s2]=green_tau(t,s1,s2,0);
        G_down[t+N*s1+N*n_site*s2]=green_tau(t,s1,s2,1);
      }
    }
  }
  if(is_thermalized()){
  	measurements.get<SignedObservable<alps::RealVectorObservable> >("G_meas_up") << (G_up*(double)sign);
  	measurements.get<SignedObservable<alps::RealVectorObservable> >("G_meas_down") << (G_down*(double)sign);
    measurements.get<RealObservable>("Sign")<<sign;
    
  }
#ifdef FOURPOINT
  throw std::logic_error("do you know what you're doing? This is the fourpoint measurement...");
  if(is_thermalized() && measure_fourpoint && (sweeps %fp_interval==0)){
    int s=Green_up.size1();
    int s3=s*s*s;
    int s2=s*s;
    static alps::numeric::vector<double> G_ijkl_uu(s*s*s, 0);
    static alps::numeric::vector<double> G_ijkl_dd(s*s*s, 0);
    static alps::numeric::vector<double> G_ijkl_du(s*s*s, 0);
    static alps::numeric::vector<double> G_ijkl_ud(s*s*s, 0);
    static int fpsteps=0;
    fpsteps++;
    for(int j=0;j<s;++j){
      double guji=Green_up(j, 0);
      double gdji=Green_down(j, 0);
      for(int k=0;k<s;++k){
        double gujk=Green_up(j, k);
        double gdjk=Green_down(j, k);
        for(int l=0;l<s;++l){
          double gulk=Green_up(l, k);
          double gdlk=Green_down(l, k);
          double guli=Green_up(l, 0);
          double gdli=Green_down(l, 0);
          G_ijkl_uu(j*s2+k*s+l)+=guji*gulk-gujk*guli; //check review formula 154
          G_ijkl_dd(j*s2+k*s+l)+=gdji*gdlk-gdjk*gdli;
          G_ijkl_ud(j*s2+k*s+l)+=guji*gdlk;;
          G_ijkl_du(j*s2+k*s+l)+=gdji*gulk;
        }
      }
    }
    int meas_int=100;
    if(fpsteps% meas_int == 0){
      std::valarray<double> G_ijkl_uu_va(s*s*s);
      std::valarray<double> G_ijkl_dd_va(s*s*s);
      std::valarray<double> G_ijkl_ud_va(s*s*s);
      std::valarray<double> G_ijkl_du_va(s*s*s);
      for(int i=0;i<s*s*s;++i){
        G_ijkl_uu_va[i]=G_ijkl_uu(i)/meas_int;
        G_ijkl_ud_va[i]=G_ijkl_ud(i)/meas_int;
        G_ijkl_du_va[i]=G_ijkl_du(i)/meas_int;
        G_ijkl_dd_va[i]=G_ijkl_dd(i)/meas_int;
      }
      
      measurements.get<SimpleRealVectorObservable>("G_fourpoint_uu")<<G_ijkl_uu_va;
      measurements.get<SimpleRealVectorObservable>("G_fourpoint_dd")<<G_ijkl_dd_va;
      measurements.get<SimpleRealVectorObservable>("G_fourpoint_ud")<<G_ijkl_ud_va;
      measurements.get<SimpleRealVectorObservable>("G_fourpoint_du")<<G_ijkl_du_va;
      
      G_ijkl_uu=0;
      G_ijkl_du=0;
      G_ijkl_ud=0;
      G_ijkl_dd=0;
    }
  }
#endif
	
  // check if precision has deteriorated
  if (check_counter%N_check == 0) {
    
    dense_matrix Green_diff_up(Green_up), Green_diff_down(Green_down);
    
    update_from_zero(Green_up, Green0_up, spins, lambda);
    update_from_zero(Green_down, Green0_down, spins, -lambda);
    
    Green_diff_up -= Green_up;
    Green_diff_down -= Green_down;
    
    if (norm_max(Green_diff_up) > tolerance || norm_max(Green_diff_down) > tolerance)
      N_check /= 2;
    else 
      N_check *= 2;
	  
    check_counter = 0;
    
  }
  
}


std::pair<matsubara_green_function_t, itime_green_function_t>HirschFyeSim::get_result() {
  int N=parms["N"];
  //double beta=parms["BETA"];
  int n_matsubara=parms["NMATSUBARA"];
  int n_site=parms["SITES"];

  itime_green_function_t green_result(N+1, n_site, 2); //one site, two spins, n+1 time points including first and last
  matsubara_green_function_t green_result_matsubara(n_matsubara, n_site, 2); //n matsubara frequencies
  
  alps::RealVectorObsevaluator G_meas_up=get_measurements()["G_meas_up"];
  alps::RealVectorObsevaluator G_meas_down=get_measurements()["G_meas_down"];  
  
  std::valarray<double> g_up_mean_vector(G_meas_up.mean());
  std::valarray<double> g_down_mean_vector(G_meas_down.mean());
  
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
    std::cout<<"end of sim, checkpointing"<<std::endl;
    std::string fns=parms["CHECKPOINT"];
    boost::filesystem::path fn(fns);
    checkpoint(fn);
  }
  for(int i=0;i<n_site;++i){
    for(int j=0;j<n_site;++j){
      for(int t=0;t<N;++t){
        green_result(t,i,j,0)=g_up_mean_vector[t+N*i+N*n_site*j];
        green_result(t,i,j,1)=g_down_mean_vector[t+N*i+N*n_site*j];
      }
      green_result(N,i,j,0)=(i==j?1-green_result(0,i,j,0):-green_result(0,i,j,0));
      green_result(N,i,j,1)=(i==j?1-green_result(0,i,j,1):-green_result(0,i,j,1));
    }
  }
  std::cout<<"green result is: "<<green_result<<std::endl;
  for(int i=0;i<n_site;++i){ //change to SC condition convention
    for(int j=0;j<n_site;++j){
      for(int t=0;t<=N;++t){
        green_result(t,i,j,0)*=-1;
        green_result(t,i,j,1)*=-1;
      }
    }
  }

  std::vector<double> densities;
  densities.resize(2);
  for(int i=0;i<n_site;++i){
    densities[0] -= green_result(N,i,i,0);
    densities[1] -= green_result(N,i,i,1);
  }
  densities[0] /= n_site;
  densities[1] /= n_site;

  boost::shared_ptr<FourierTransformer> fourier_ptr;
  FourierTransformer::generate_transformer_U(parms, fourier_ptr, densities);
  fourier_ptr->forward_ft(green_result, green_result_matsubara);
  alps::hdf5::archive ar(static_cast<std::string>(parms["OUTFILE"]), "a");
  green_result.write_hdf5(ar, "/G_tau");
  green_result_matsubara.write_hdf5(ar, "/G_omega");
  return std::make_pair(green_result_matsubara, green_result);
}

///take a Green's function in vector form G(\tau-\tau') and create a matrix G(\tau, \tau') out of it
void HirschFyeRun::green_matrix_from_vector(const itime_green_function_t & bare_green_function,\
                                            dense_matrix & green_matrix_up, dense_matrix & green_matrix_down)const{
  
  green_matrix_up  .clear();
  green_matrix_down.clear();
  for (int s1=0; s1<n_site; s1++) {
    for (int s2=0; s2<n_site; s2++) {	  
      for (int i=0; i<N; i++) {
        for (int j=i; j<N; j++) {	  
          green_matrix_up  (s1*N+j,s2*N+i) = bare_green_function(j-i, s1, s2, 0); //green_vector[j-i];
          green_matrix_down(s1*N+j,s2*N+i) = bare_green_function(j-i, s1, s2, 1); //green_vector[j-i];
        }
        for (int j=i+1; j<N; j++) {//consider antisymmetry (Fermions)
          green_matrix_up  (s1*N+i,s2*N+j) = -bare_green_function(N-(j-i), s1, s2, 0); //-green_vector[N-(j-i)];
          green_matrix_down(s1*N+i,s2*N+j) = -bare_green_function(N-(j-i), s1, s2, 1); //-green_vector[N-(j-i)];
        }
      }
    }
  }
}
///take a Green's function in matrix form and transform it back into vector form.
void HirschFyeRun::green_vector_from_matrix(itime_green_function_t &green_tau, const dense_matrix & green_matrix_up,\
                                            const dense_matrix & green_matrix_down)const{
  green_tau.clear();
  
  for (int s1=0;s1<n_site;++s1){
    for(int s2=0;s2<n_site;++s2){
      for (int i=0; i<N; i++) {
        for (int j=i; j<N; j++) {
          green_tau(j-i,s1,s2,0) += green_matrix_up  (j+s1*N,i+s2*N);
          green_tau(j-i,s1,s2,1) += green_matrix_down(j+s1*N,i+s2*N);
          if (j>i) {//consider antisymmetry (Fermions)
            green_tau(N-(j-i),s1,s2,0) -= green_matrix_up  (i+s1*N,j+s2*N);
            green_tau(N-(j-i),s1,s2,1) -= green_matrix_down(i+s1*N,j+s2*N);
          }
        }
      }
      for(int i=0;i<N;++i){
        green_tau(i, s1, s2, 0)/=N;
        green_tau(i, s1, s2, 1)/=N;
      }
    }
  }
}


