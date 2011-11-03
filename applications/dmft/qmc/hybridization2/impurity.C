/************************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2011 by Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Emanuel Gull <gull@phys.columbia.edu>,
 *                              Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
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
 ************************************************************************************/

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
N_nn(static_cast<int>(parms.value_or_default("N_nn", 0))),                                       //number of tau-points on which density density correlator is measured
N_w(static_cast<int>(parms.value_or_default("N_w", 0))),                                         //number of Matsubara frequencies for gw
N_w2(static_cast<int>(parms.value_or_default("N_w2", 0))),                                       //number of Matsubara frequencies for two-particle Green's function
N_W(static_cast<int>(parms.value_or_default("N_W", 0))),                                         //number of bosonic frequncies for the two-particle Green's function
N_l(static_cast<int>(parms.value_or_default("N_l", 0))),                                         //number of Legendre coefficients
N_order(static_cast<int>(parms["N_ORDER"])),                                                     //number of orders that are measured for the order histogram
G_meas(static_cast<int>(parms["FLAVORS"])*(static_cast<int>(parms["N"])+1)),                     //Measured Green's function
N_meas(static_cast<int>(parms["N_MEAS"])),                                                       //number of updates per measurement
N_shift(static_cast<int>(parms.value_or_default("N_SHIFT", 0))),                                 //number of times a segment start / end point is shifted per insertion/removal update
N_swap(static_cast<int>(parms.value_or_default("N_SWAP", 0))),                                   //number of times an orbital swap move is attempted.
F(static_cast<int>(parms["FLAVORS"]), std::vector<double>(static_cast<int>(parms["N"])+1,0.)),   //hybridization function
MEASURE_gw(static_cast<int>(parms.value_or_default("MEASURE_gw", 0))),                           //measure Green's function on Matsubara frequencies
MEASURE_fw(static_cast<int>(parms.value_or_default("MEASURE_fw", 0))),                           //measure Sigma*G on Matsubara frequencies
MEASURE_gl(static_cast<int>(parms.value_or_default("MEASURE_gl", 0))),                           //measure Green's function in Legendre basis
MEASURE_fl(static_cast<int>(parms.value_or_default("MEASURE_fl", 0))),                           //measure Sigma*G in Legendre basis
MEASURE_g2w(static_cast<int>(parms.value_or_default("MEASURE_g2w", 0))),                         //measure two-particle Green's function
MEASURE_hw(static_cast<int>(parms.value_or_default("MEASURE_hw", 0))),                           //measure correlator for vertex function
MEASURE_nnt(static_cast<int>(parms.value_or_default("MEASURE_nnt", 0))),                         //measure density-density correlation function
MEASURE_nn(static_cast<int>(parms.value_or_default("MEASURE_nn", 0)))                            //measure density-density correlation function at equal times
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
    std::cout<<"*******************The interaction (\"U\")-matrix is: ************************"<<std::endl;
    for(int i=0;i<FLAVORS;++i){ for(int j=0;j<FLAVORS;++j){ std::cout<<u(i,j)<<" "; } std::cout<<std::endl;}
    std::cout<<"****************************************************************************"<<std::endl;
    
    std::cout<<"*****************The effective chemical potential is: **********************"<<std::endl;
    for(int i=0;i<FLAVORS;++i){ std::cout<<"effective mu/t for band: "<<i<<" "<<mu_e[i]<<std::endl; }
    std::cout<<"****************************************************************************"<<std::endl;
    
    std::cout<<"*****************The hybridization function is: ****************************"<<std::endl;
    for(int i=0;i<10;++i){ std::cout<<i<<" "; for(int j=0;j<FLAVORS;++j){ std::cout<<F[j][i]<<" ";} std::cout<<std::endl; }
    std::cout<<"... *** etc *** ...\n";
    std::cout<<F[0].size()-1<<" "; for(int j=0;j<FLAVORS;++j){ std::cout<<F[j].back()<<" ";} std::cout<<std::endl;
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
  create_measurements();
  measurements.reset(true);
  
  resize_measurement_vectors(crank);
  std::cout<<"process " << crank << " starting simulation"<<std::endl;
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
    }
  //these are measurements that are done every single step (and buffered before stored in ALPS.)
  set_measurement_vectors();
  }
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
    throw(std::invalid_argument(std::string("U-matrix parameter \'") + parms["U_MATRIX"].str() + "\' not specified or not pointing to a valid file. File could not be opened."));
  }
  for (int i=0; i<FLAVORS; i++) {
    if(!infile_u.good()){
      throw(std::invalid_argument(std::string("U-matrix file ") + parms["U_MATRIX"].str() + " not good. Probably too few lines."));
    }
    for (int j=0; j<FLAVORS; j++) {
      if(!infile_u.good()){
        throw(std::invalid_argument(std::string("U-matrix file ") + parms["U_MATRIX"].str() + " not good. Probably too few columns."));
      }
      infile_u >> u(i,j);
    }
  }
  std::string mu_dc_solver_string=parms.value_or_default("EPSILOND_VECTOR","mu_dc_solver");
  ifstream mu_dc_solver(mu_dc_solver_string.c_str());
  int dummy;
  if(!mu_dc_solver.good()) throw(std::invalid_argument("file for level shifts "+mu_dc_solver_string+" could not be found. Check paramter EPSILOND_VECTOR. Aborting"));
  for(int i=0;i<FLAVORS;++i){
    if(!mu_dc_solver.good()) throw(std::invalid_argument("mu dc solver file "+mu_dc_solver_string+" does not have enough lines? check format of mu file"));
    mu_dc_solver>>dummy>>mu_e[i]>>std::ws;
  }

  // read F from file
  if(parms.defined("F")){
    std::string fname=boost::lexical_cast<std::string>(parms["F"]);
    size_t found=fname.find(".h5",fname.size()-3);
    if(found!=string::npos){//attempt to read from h5 archive
      alps::hdf5::archive ar(fname, alps::hdf5::archive::READ);
      int N_, FLAVORS_;
      ar>>alps::make_pvp("/N",N_);
      ar>>alps::make_pvp("/FLAVORS",FLAVORS_);
      if(N_!=N) throw::std::invalid_argument(std::string("bad file ") + parms["F"].str() + "wrong number of time slices");
      if(FLAVORS_!=FLAVORS) throw::std::invalid_argument(std::string("bad file ") + parms["F"].str() + "wrong number of flavors");
      std::vector<double> tmp(FLAVORS*(N+1));
      ar>>alps::make_pvp("/data",tmp);
      for(int j=0; j<FLAVORS; j++)
        for(int i=0; i<N+1; i++)
          F[j][i]=tmp[j*(N+1)+i];
      tmp.clear();
    }
    else{
      ifstream infile(fname.c_str());
      if(!infile.good()){
        throw(std::invalid_argument(std::string("could not open file ") + parms["F"].str() + "for F function"));
      }
      for (int i=0; i<N+1; i++) {
        if(!infile.good()){
          throw(std::invalid_argument(std::string("bad file ") + parms["F"].str() + "probably wrong number of lines"));
        }
        double dummy;
        infile >> dummy; 
        for (int j=0; j<FLAVORS; j++){
          if(!infile.good()){
            throw(std::invalid_argument(std::string("bad file ") + parms["F"].str() + "probably wrong number of columns"));
          }
          infile >> F[j][i];
        }
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
