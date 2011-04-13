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
#include "moves.h"
#include "fouriertransform.h"
#include <alps/alea.h>
#include <alps/hdf5.hpp>

using namespace std;
using namespace alps;


hybridization::hybridization(const parameters_type &parms, int crank)
: alps::mcbase(parms, crank),
sweeps(0),                                                                                       //Sweeps currently done
thermalization_sweeps(static_cast<int>(parms["THERMALIZATION"])),                                //Sweeps to be done for thermalization
total_sweeps(static_cast<int>(parms["SWEEPS"])),                                                 //Sweeps to be done in total
FLAVORS(static_cast<int>(parms["FLAVORS"])),                                                     //number of orbitals
mu_e(static_cast<int>(parms["FLAVORS"])),                                                        //level energy of each segment
u(static_cast<int>(parms["FLAVORS"]),static_cast<int>(parms["FLAVORS"])),                        //interaction matrix
t(static_cast<double>(parms["t"])),                                                              //overall energy unit.
full_line(static_cast<int>(parms["FLAVORS"]),0),                                                 //vector, is full_line[j] if orbital j is fully occuppied for all times.
sign(static_cast<int>(parms["FLAVORS"]), 1.),                                                    //fermionic sign. plus or minus one.
BETA(static_cast<double>(parms["BETA"])),                                                        //inverse temperature
N(static_cast<int>(parms["N"])),                                                                 //number of points on which Green's function is known
N_order(static_cast<int>(parms["N_ORDER"])),                                                     //number of orders that are measured for the order histogram
G_meas(static_cast<int>(parms["FLAVORS"])*(static_cast<int>(parms["N"])+1)),                     //Measured Green's function
G_meas_imp(static_cast<int>(parms["FLAVORS"])*(static_cast<int>(parms["N"])+1)),                 ///---PHILIPP CHECK///
N_corr(static_cast<int>(parms.value_or_default("N_CORR", 0))),                                   //number of tau-points on which density density correlator is measured
N_meas(static_cast<int>(parms["N_MEAS"])),                                                       //number of updates per measurement
N_shift(static_cast<int>(parms.value_or_default("N_SHIFT", 0))),                                 //number of times a segment start / end point is shifted per insertion/removal update
N_swap(static_cast<int>(parms.value_or_default("N_SWAP", 0))),                                   //number of times an orbital swap move is attempted.
F(static_cast<int>(parms["FLAVORS"]), std::vector<double>(static_cast<int>(parms["N"])+1,0.))
{
  
  /*reading parameters from the self consistency*/
  if(parms.defined("EXTERNAL_INPUT_DATA") && ((bool)(parms["EXTERNAL_INPUT_DATA"])==true)){
    read_external_input_data(parms);
  }else if(parms.defined("OMEGA_LOOP") && (bool)(parms["OMEGA_LOOP"])==true){
    read_alps_framework_input_data_omega(parms);
  }else{
    read_alps_framework_input_data_tau(parms);
  }
  
  if(crank==0){
    std::cout<<"*****************The interaction (\"U\")-matrix is: ************************"<<std::endl;
    for(int i=0;i<FLAVORS;++i){ for(int j=0;j<FLAVORS;++j){ std::cout<<u(i,j)<<" "; } std::cout<<std::endl;}
    std::cout<<"****************************************************************************"<<std::endl;
    
    std::cout<<"*****************The effective chemical potential is: **********************"<<std::endl;
    for(int i=0;i<FLAVORS;++i){ std::cout<<"effective mu/t for band: "<<i<<" "<<mu_e[i]<<std::endl; }
    std::cout<<"****************************************************************************"<<std::endl;
    
    std::cout<<"*****************The hybridization function is: ****************************"<<std::endl;
    for(int i=0;i<10;++i){ std::cout<<i<<" "; for(int j=0;j<FLAVORS;++j){ std::cout<<F[j][i]<<" ";} std::cout<<std::endl; }
    std::cout<<"... *** etc *** ...\n";
    std::cout<<"****************************************************************************"<<std::endl;
  }
  segments.resize(FLAVORS);
  M.resize(FLAVORS);
  // initialize list of segments
  for (int i=0; i<FLAVORS; i++) {
    segments[i].clear();
    M[i].clear();
  }
  
  // create measurement objects
  results.create_RealVectorObservable("n"); 
  results.create_RealVectorObservable("order");
  results.create_RealVectorObservable("Greens");
  results.create_RealVectorObservable("Greens_imp");
  results.create_RealVectorObservable("nn_corr");
  results.create_RealVectorObservable("nn_corr_equalt");
  //results.create_RealVectorObservable("overlap");
  results.create_RealObservable("sign"); 
  results.create_RealVectorObservable("matrix_size"); 
  
  results.reset(true);
  
  //resizing measurement vectors
  n_meas.resize(FLAVORS,0.);
  order_meas.resize(N_order*FLAVORS,0.);
  nn_corr_meas.resize(FLAVORS*(FLAVORS+1)/2*(N_corr+1),0.);
  nn_corr_equalt_meas.resize(FLAVORS*FLAVORS,0.);
  n_vectors.resize(FLAVORS);
  matrix_size.resize(FLAVORS, 0.);
  
  std::cout<<"starting simulation"<<std::endl;	
}


bool hybridization::is_thermalized() const
{
  return (sweeps >= thermalization_sweeps);
}

double hybridization::work_done() const
{
  return (is_thermalized() ? (sweeps-thermalization_sweeps)/double(total_sweeps) : 0.);
  
}

void hybridization::do_update()
{
  
  // increment sweep count
  sweeps++;

  segment_container_t::iterator it1, it2;
  double s=1;
  for (int i=0; i<N_meas; i++) {
   //do an update for each FLAVORS
    for (int j=0; j<FLAVORS; j++) {
      if (segments[j].size() == 0) {
        // insert or remove full line
        insert_remove_full_line(random, mu_e[j], BETA, u,full_line[j], segments, full_line, j);
      }
      insert_remove_antisegment(random, BETA*random(), BETA, mu_e[j], u,F[j], full_line[j], segments[j], M[j], sign[j], segments, full_line, j);
      if (!full_line[j]) {
        // local update
        insert_remove_segment(random, BETA*random(), BETA, mu_e[j], u,F[j], segments[j], M[j], sign[j], segments, full_line, j);
        // shift segments
        for (int k=0; k<N_shift; k++) 
          shift_segment(random, segments[j], BETA, mu_e[j], u, F[j], M[j], sign[j], segments, full_line, j);
      }
      
      // measure perturbation order
      if (segments[j].size()<(std::size_t)N_order)
        order_meas[j*N_order+segments[j].size()]++;
      
      
      // measure Green functions
      if (segments[j].size()>0) {
        for (int i=0; i<M[j].size1(); i++) {
          (i==0 ? it1 = segments[j].begin() : it1++);
          for (int k=0; k<M[j].size1(); k++) {
            (k==0 ? it2 = segments[j].begin() : it2++);
            if (M[j](k,i)!=0) {
              double argument = it1->t_end()-it2->t_start();
              double bubble_sign=1;
              if (argument > 0) {
                bubble_sign = 1;
              }
              else {
                bubble_sign = -1;
                argument += BETA;
              }
              int index = (int)(argument/BETA*N+0.5);
              G_meas[j*(N+1)+index] += M[j](k,i)*bubble_sign;
            }
          }
        }
      }
      
      s *= sign[j];
      times full_segment(0,BETA);
      n_meas[j] += compute_overlap(full_segment, segments[j], full_line[j], BETA)/BETA;
      
    }
    
    //these are measurement things that are done every single step (and buffered before stored in ALPS.)
    sign_meas += s;
    for(int i=0;i<FLAVORS;++i){
      matrix_size[i]+=M[i].size();
    }

    // measure density correlation functions 
    segment_container_t::iterator it;	
    for (int flavor=0; flavor<FLAVORS; ++flavor) {
      n_vectors[flavor].resize(N_corr+1, 1);
      if (segments[flavor].size()==0) {
        if (full_line[flavor]==0) {
          for (std::size_t i=0; i<n_vectors[flavor].size(); ++i)
            n_vectors[flavor][i]=0;
        }
      }
      else {
        it=segments[flavor].end(); it--;
        if (it->t_end()<it->t_start()) 
          n_vectors[flavor][0]=1;
        else
          n_vectors[flavor][0]=0;
        
        // mark segment start and end points
        int index;
        for (it=segments[flavor].begin(); it!=segments[flavor].end(); it++) {
          index = (int)(it->t_start()/BETA*N_corr+1);
          n_vectors[flavor][index] *= -1;
          index = (int)(it->t_end()/BETA*N_corr+1);
          n_vectors[flavor][index] *= -1;
        }
        // fill vector with occupation number
        for (std::size_t i=1; i<n_vectors[flavor].size(); i++) {
          if (n_vectors[flavor][i]==-1)
            n_vectors[flavor][i]=1-n_vectors[flavor][i-1];
          else
            n_vectors[flavor][i]=n_vectors[flavor][i-1];
        }
      }  
    }
    
    // compute n(\tau)n(0)
    memset(&(nn_corr_meas[0]), 0, FLAVORS*(FLAVORS+1)/2*(N_corr+1)*sizeof(double));
    memset(&(nn_corr_equalt_meas[0]), 0, FLAVORS*FLAVORS*sizeof(double));
    for (int flavor1=0; flavor1<FLAVORS; ++flavor1) {
      for (int flavor2=0; flavor2<FLAVORS; ++flavor2) {
        for(int i=0;i<N_corr;++i)
          nn_corr_equalt_meas[flavor1*FLAVORS+flavor2] += n_vectors[flavor1][i]*n_vectors[flavor2][i];
      }
    }
    int position=0;
    for (int flavor1=0; flavor1<FLAVORS; ++flavor1) {
      for (int flavor2=0; flavor2<=flavor1; ++flavor2) {
        for (int i=0; i<N_corr+1; ++i) {
          int index = 0; // measure only equal time correlations here
          int j=i+index;
          if (j>N_corr) j -= N_corr;   
          nn_corr_meas[position+index] += n_vectors[flavor1][i]*n_vectors[flavor2][j];
        }
        position += (N_corr+1);
      }
    }
  }
}
void hybridization::do_measurements(){
  if(is_thermalized()){
    nn_corr_meas /= (N_corr+1);
    results["nn_corr"] << nn_corr_meas;
    nn_corr_equalt_meas/=(double)N_corr;
    results["nn_corr_equalt"] << (nn_corr_equalt_meas);
    
    order_meas /= N_meas;
    results["order"] << order_meas;
    
    G_meas *= (1.*N)/N_meas/(BETA*BETA);
    results["Greens"] << (G_meas); //*sign;
    
    G_meas_imp *= (1.*N)/N_meas/(BETA*BETA);
    results["Greens_imp"] << (G_meas_imp); //*sign;
    
    sign_meas /= N_meas;
    results["sign"] << sign_meas;
    //if(sign_meas != 1.) throw std::runtime_error("negative sign encountered. The current code is not able to deal with this (in the segment representation such a situation should not arise for diagonal hybridizations). Do you know what you're doing?!?");
    
    n_meas /= N_meas;
    results["n"]<< n_meas;
    
    matrix_size /= N_meas;
    results["matrix_size"]<< matrix_size;
  }
  memset(&(G_meas[0]),0, G_meas.size()*sizeof(double));
  memset(&(G_meas_imp[0]),0, G_meas.size()*sizeof(double));
  memset(&(n_meas[0]),0, n_meas.size()*sizeof(double));
  memset(&(matrix_size[0]),0, matrix_size.size()*sizeof(double));
  sign_meas=0;
  
}


std::ostream &operator<<(std::ostream & os, segment_container_t segments){
  for(segment_container_t::iterator it=segments.begin(); it !=segments.end();++it){
    os<<it->t_start()<<" "<<it->t_end()<<" ";
  }
  std::cout<<std::endl;
  return os;
}

void hybridization::read_external_input_data(const parameters_type &parms){
  // define interactions between flavors. Read in interaction file.
  u.resize(FLAVORS, FLAVORS);
  ifstream infile_u(boost::lexical_cast<std::string>(parms["U_MATRIX"]).c_str());
  if(!infile_u.good()){
    throw(std::invalid_argument(std::string("U-matrix parameter ") + std::string(parms["U_MATRIX"])+" invalid. File could not be opened."));
  }
  for (int i=0; i<FLAVORS; i++) {
    if(!infile_u.good()){
      throw(std::invalid_argument(std::string("U-matrix file ") + std::string(parms["U_MATRIX"])+" not good. Probably too few lines."));
    }
    for (int j=0; j<FLAVORS; j++) {
      if(!infile_u.good()){
        throw(std::invalid_argument(std::string("U-matrix file ") + std::string(parms["U_MATRIX"])+" not good. Probably too few columns."));
      }
      infile_u >> u(i,j);
    }
  }
  ifstream mu_dc_solver("mu_dc_solver");
  int dummy;
  if(!mu_dc_solver.good()) throw(std::invalid_argument("file mu_dc_solver could not be found. Aborting"));
  for(int i=0;i<FLAVORS;++i){
    if(!mu_dc_solver.good()) throw(std::invalid_argument("mu dc solver file has not enough lines?"));
    mu_dc_solver>>dummy>>mu_e[i]>>std::ws;
    mu_e[i]*=t; //we need to adjust this to the right onints
  }
  
  
  // read F from file
  if(parms.defined("F")){
    ifstream infile(boost::lexical_cast<std::string>(parms["F"]).c_str());
    if(!infile.good()){
      throw(std::invalid_argument(std::string("could not open file ") +std::string(parms["F"])+"for F function"));
    }
    for (int i=0; i<N+1; i++) {
      if(!infile.good()){
        throw(std::invalid_argument(std::string("bad file ") +std::string(parms["F"])+"probably wrong number of lines"));
      }
      double dummy;
      infile >> dummy; 
      for (int j=0; j<FLAVORS; j++){
        if(!infile.good()){
          throw(std::invalid_argument(std::string("bad file ")+std::string(parms["F"])+"probably wrong number of columns"));
        }
        infile >> F[j][i];
      }
    }
  }
}
void hybridization::read_alps_framework_input_data_omega(const parameters_type &parms){
  
  if(!parms.defined("U")) throw std::runtime_error("parameter U (for the interaction strength) is missing!");
  double U=parms["U"];
  double mu=parms["MU"];
  for(int i=0;i<FLAVORS;++i){
    for(int j=0;j<FLAVORS;++j){
      u(i,j)=(i==j?0.:U);
    }
    mu_e[i]=mu+U/2; //discuss about more general Us and Js and general U shift.
  }
  
  matsubara_green_function_t bare_green_matsubara(N, 1, FLAVORS);
  matsubara_green_function_t f_twiddle_matsubara(N, 1, FLAVORS);
  itime_green_function_t f_twiddle_itime(N+1, 1, FLAVORS);
  {
    alps::hdf5::archive ar(parms["INFILE"], alps::hdf5::archive::READ);
    bare_green_matsubara.read_hdf5(ar,"/G0");
  }
  //find the second moment of the band structure
  double epssqav ;
  if (parms.defined("DOSFILE")) {
    if (!parms.defined("EPSSQAV")) {
      throw std::logic_error("error: you specify a DOS file, please also specify the second moment of the band structure EPSSQAV!");
    } else {
      epssqav = parms["EPSSQAV"];
    }
  }else{
    t=parms["t"]; //this is the default: semicircular density of states
    epssqav = t * t; //...and its moment.
  }
  
  FFunctionFourierTransformer Fourier(BETA , 0, epssqav , FLAVORS, 1);
  for (int f = 0; f < FLAVORS; ++f) {
    for (int i = 0; i < N; ++i) {
      f_twiddle_matsubara(i, f) = -1. / bare_green_matsubara(i, f) + (std::complex < double >(mu, (2. * i + 1) * M_PI / BETA));
    }
  }
  Fourier.backward_ft(f_twiddle_itime, f_twiddle_matsubara);
  for (int f = 0; f < FLAVORS; ++f) {
    for (int i = 0; i < N+1; ++i) {
      F[f][i] = -f_twiddle_itime(N - i, f);
    }
  }
}
/*************************************************************************************************/
/**This is a very special case that makes the self consistency easy: in the case of a BETHE    ***/
/**lattice (semicircular density of states) we can compute the hybridization function F        ***/
/**directly from G via F = -t*t*G (note that the usual convention has Delta(omega)=-F(-omega)  ***/
/**hence Delta =  t*t*G (t is the hopping matrix element of the semicircular DOS, the quarter  ***/
/**bandwidth.
 /*************************************************************************************************/
void hybridization::read_alps_framework_input_data_tau(const parameters_type &parms){
  itime_green_function_t green_itime(N+1, 1, FLAVORS);
  {
    alps::hdf5::archive ar(parms["INFILE"], alps::hdf5::archive::READ);
    green_itime.read_hdf5(ar,"/G0");
  }
  for (int f = 0; f < FLAVORS; ++f) {
    int orbital=f/2;
    std::stringstream tname; tname<<"t"<<orbital;
    if(parms.defined(tname.str())) {
      t=parms[tname.str()]; 
      std::cout<<"orbital: "<<f/2<<" flavor: "<<f<<" using hopping parameter (quarter bandwidth) t: "<<t<<std::endl;
    }
    else{
      t=parms["t"];
      if(!parms.defined("t")) throw std::runtime_error("parameter t (for the quarter bandwidth) is missing!");
    }
    for (int i = 0; i<= N; ++i) {
      F[f][i] = -t * t * green_itime(N - i, f);  //this is the self consistency condition, for Bethe lattice!
    }
  }
  
  if(!parms.defined("U")) throw std::runtime_error("parameter U (for the interaction strength) is missing!");
  double U=parms["U"];
  double mu=parms["MU"];
  for(int i=0;i<FLAVORS;++i){
    for(int j=0;j<FLAVORS;++j){
      u(i,j)=(i==j?0.:U);
    }
    mu_e[i]=mu+U/2; //discuss about more general Us and Js and general U shift.
  }
  std::cout << "U is: " << u << std::endl;
}
