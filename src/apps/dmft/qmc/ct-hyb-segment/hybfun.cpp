/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                   
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#include"hybfun.hpp"
#include <sstream>

// Internal to ALPS::cthyb_segment: keeps its classes apart from ALPS::cthyb's.
namespace cthyb_segment {
//construct a hybridization function. ntime: number of time slices.
//noffidag_orbitals: number of offdiagonal orbitals. ndiag_orbitals: number
//of diagonal orbitals. 
hybfun::hybfun(const alps::params &p):
green_function<double>(p["N"].as<int>()+1, 1, p["FLAVORS"])
{
  beta_=p["BETA"];

  //read in Green's function from a file
  read_hybridization_function(p);
  hybridization_function_sanity_check();

  extern int global_mpi_rank;
  if(global_mpi_rank==0){
    std::cout<<*this<<std::endl;
  }
}

void hybfun::hybridization_function_sanity_check(void){
for(std::size_t i=0; i<ntime();++i)
  for(std::size_t j=0; j<nflavor();++j)
    if(operator()(i,j)>0.) throw std::invalid_argument("Problem with hybridization function: Delta(\\tau) > 0. Delta should always be negative!");
}

//this routine reads in the hybridization function, either from a text file or from an hdf5 file (for easy passing of binary data).
//In case of text files the file format is index - hyb_1 - hyb2 - hyb3 - ... in columns that go from time=0 to time=beta. Note that
//the hybridization function is in imaginary time and always positive between zero and \beta.
void hybfun::read_hybridization_function(const alps::params &p){
  if(!p.exists("cthyb.DELTA")){
    if(p.exists("cthyb.DMFT_FRAMEWORK") && p["cthyb.DMFT_FRAMEWORK"] && p.exists("solver.INFILE_H5GF")) {
      read_hybridization_from_h5gf(p);
      return;
    }
    else {
      throw(std::invalid_argument(std::string("Parameter DELTA missing, filename for hybridization function not specified.")));
    }
  }
  std::string fname=p["cthyb.DELTA"];
  if(p.exists("cthyb.DELTA_IN_HDF5") && p["cthyb.DELTA_IN_HDF5"]){//attempt to read from h5 archive
    alps::hdf5::archive ar(fname, alps::hdf5::archive::READ);
    if(p.exists("cthyb.DMFT_FRAMEWORK") && p["cthyb.DMFT_FRAMEWORK"]){//read in as green_function
      read_hdf5(ar,"/Delta");
    }
    else if(ar.is_data("/Delta")){//plain hdf5: all flavors in one array
      std::vector<double> tmp(nflavor()*ntime());
      ar>>alps::make_pvp("/Delta",tmp);
      for(std::size_t j=0; j<nflavor(); j++)
        for(std::size_t i=0; i<ntime(); i++)
          operator()(i,j)=tmp[j*ntime()+i];
      tmp.clear();
    }
    else{//plain hdf5: one array /Delta_j per flavor, as ALPS' hybridization solver reads it
      std::vector<double> tmp(ntime());
      for(std::size_t j=0; j<nflavor(); j++){
        std::stringstream path; path<<"/Delta_"<<j;
        ar>>alps::make_pvp(path.str(),tmp);
        for(std::size_t i=0; i<ntime(); i++)
          operator()(i,j)=tmp[i];
      }
    }
  }
  else{//read from text file
    std::ifstream infile(fname.c_str());
    if(!infile.good()){
      throw(std::invalid_argument(std::string("could not open hybridization file (text format) ") + p["cthyb.DELTA"].as<std::string>() + " for Delta function"));
    }
    for (std::size_t i=0; i<ntime(); i++) {
      if(!infile.good()){
        throw(std::invalid_argument(std::string("could not read hybridization file (text format) ") + p["cthyb.DELTA"].as<std::string>() + ". probably wrong number of lines"));
      }
      double dummy;
      infile >> dummy;
      //std::cout<<i<<" "<<ntime()<<" "<<"dummy: "<<dummy<<std::endl;
      for (std::size_t j=0; j<nflavor(); j++){
        if(!infile.good()){
          throw(std::invalid_argument(std::string("could not read hybridization file (text format) ") + p["cthyb.DELTA"].as<std::string>() + ". probably wrong number of columns"));
        }
        double delta;
        infile >> delta;
        operator()(i,j)=delta;
      }
    }
  }
}

void hybfun::read_hybridization_from_h5gf(const alps::params &p) {
  double mu = p["MU"];
  int n_tau = p["N"];
  int n_orbitals = p["FLAVORS"];
  int n_matsubara = p["NMATSUBARA"];
  alps::hdf5::archive ar(p["solver.INFILE_H5GF"], alps::hdf5::archive::READ);
  alps::gf::omega_sigma_gf_with_tail G0_omega(alps::gf::omega_sigma_gf(alps::gf::matsubara_positive_mesh(beta_, n_matsubara), alps::gf::index_mesh(n_orbitals)));
  alps::gf::omega_sigma_gf_with_tail Delta_w(alps::gf::omega_sigma_gf(alps::gf::matsubara_positive_mesh(beta_, n_matsubara), alps::gf::index_mesh(n_orbitals)));
  alps::gf::itime_sigma_gf_with_tail Delta(alps::gf::itime_sigma_gf(alps::gf::itime_mesh(beta_, n_tau+1), alps::gf::index_mesh(n_orbitals)));
  G0_omega.load(ar, "/G0");
  std::cout<<"successfully loaded G0."<<std::endl;

  for (alps::gf::matsubara_index i(0); i<G0_omega.mesh1().extent(); ++i) {
    for (alps::gf::index s(0); s < G0_omega.mesh2().extent(); ++s) {
      std::complex<double> iw(0., (2 * i() + 1) * M_PI / beta_);
      Delta_w(i, s) = iw + mu - 1.0/G0_omega(i, s);
    }
  }
  
  //compute the tail of the hybridization function from the tail of the Green's function.
  //use: Delta = Simplify[  i\[Omega] + \[Mu] - \[Epsilon] -  1/(1/i\[Omega] + c2/i\[Omega]^2 + c3/i\[Omega]^3 + c4/i\[Omega]^4)]
  //to get: Series[Delta,{i\[Omega],Infinity,3}]=(c2-\[Epsilon]+\[Mu])+(-c2^2+c3)/i\[Omega]+(c2^3-2 c2 c3+c4)/i\[Omega]^2+(-c2^4+3 c2^2 c3-c3^2-2 c2 c4)/i\[Omega]^3+O[1/i\[Omega]]^4
  //thus if we have only c2 and c3 we can get c1 of Delta, but we do need to have c3.
  typedef alps::gf::one_index_gf<double, alps::gf::index_mesh> density_matrix_type;
  density_matrix_type tail=density_matrix_type(alps::gf::index_mesh(2));
  tail.initialize();
  for (alps::gf::index s(0); s < G0_omega.mesh2().extent(); ++s) {
    double gf_c2=G0_omega.tail(2)(s);
    double gf_c3=G0_omega.tail(3)(s);
    tail(s) = -gf_c2*gf_c2+gf_c3;
    std::cout<<" spin: "<<s()<<" tail1: "<<G0_omega.tail(1)(s)<<" tail2: "<<G0_omega.tail(2)(s)<<" tail3: "<<G0_omega.tail(3)(s)<<" Delta tail: "<<tail(s)<<std::endl;
  }
  Delta_w.set_tail(1,tail);
  alps::gf::fourier_frequency_to_time(Delta_w, Delta);
  for (alps::gf::itime_index i(0); i<Delta.mesh1().extent(); ++i) {
    for (alps::gf::index s(0); s<Delta.mesh2().extent(); ++s) {
      operator()(i(), s()) = Delta(i, s);
    }
  }

}

std::ostream &operator<<(std::ostream &os, const hybfun &hyb){
  os<<"the hybridization function is: "<<std::endl;
  for(int i=0;i<10;++i){ std::cout<<i<<" "; for(std::size_t j=0;j<hyb.nflavor();++j){ std::cout<<hyb(i,j)<<" ";} std::cout<<std::endl; }
  os<<"... *** etc *** ...\n";
  os<<hyb.ntime()-1<<" "; for(std::size_t j=0;j<hyb.nflavor();++j){ std::cout<<hyb(hyb.ntime()-1,j)<<" ";} std::cout<<std::endl;
  return os;
}

//linear interpolation of the hybridization function. 
double hybfun::interpolate(double time, int orbital) const{

  ///TODO: do a spline here, write a better/faster interpolation routine.
  double sign=1;
  if (time<0) {
    time += beta_;
    sign=-1;
  }

  //this is the overall flip of Delta (in comparison to F)
  sign*=-1;
  //the code takes Delta as input, but internally works with F(tau)=-Delta(beta-tau). See Philipp's old paper
  time=beta_-time;

  double n = time/beta_*(ntime()-1);
  int n_lower = (int)n; // interpolate linearly between n_lower and n_lower+1

  return sign*(operator()(n_lower,orbital) + (n-n_lower)*(operator()(n_lower+1,orbital)-operator()(n_lower,orbital)));
}

}  // namespace cthyb_segment
