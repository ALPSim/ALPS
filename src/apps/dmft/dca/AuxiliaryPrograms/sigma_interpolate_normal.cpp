// Copyright (C) 2005-2026 by Emanuel Gull
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT

#include<iostream>
#include<fstream>
#include<stdexcept>
#include<vector>
#include<complex>
#include<cmath>
#include<cstdlib>


#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "cluster.h"
#include<Eigen/Dense>
#include"blasheader.h"

//take a vector and interpolate it using GSL in the middle, add constants to the back, and extrapolate parabolically for the first few frequencies if they are outside of the known interval
std::vector<double> interpolate_vector(const std::vector<double> &oldxvalues, const std::vector<double> &oldyvalues, const std::vector<double> &newxvalues);

void define_parameters_commandline(alps::params &p){
  p.define<std::string>("pfile_old", "old parameter file");
  p.define<std::string>("pfile_new", "new parameter file");
  p.define<std::string>("sigma", "self energy file (in GFTools hdf5 format)");
  p.define<std::string>("path", "/Sigma","hdf5 path to self energy");
  p.define<std::string>("g0", "hdf5 path to the G0 obtained with the interpolated self-energy");
}
void define_dca_parameters(alps::params &p);
void define_selfconsistency_parameters(alps::params &p);
void define_parameters_pfile(alps::params &p){
  define_selfconsistency_parameters(p);
  define_dca_parameters(p);
}

int main(int argc, char **argv){
  alps::params p(argc, (const char**)argv);
  define_parameters_commandline(p);
  if (p.help_requested(std::cout)) {
    exit(0);
  }

  //parse parameters: open old and new parameter file
  alps::params parms_old(p["pfile_old"].as<std::string>());
  alps::params parms_new(p["pfile_new"].as<std::string>());
  define_parameters_pfile(parms_old);
  define_parameters_pfile(parms_new);
  
  //general error check
  double beta_old=parms_old["BETA"];
  double mu_old=parms_old["MU"];
  double h_old=parms_old["H"];
  double nsites_old=parms_old["dca.SITES"];
  int nfreq_old = parms_old["NMATSUBARA"];
  
  double beta_new=parms_new["BETA"];
  double mu_new=parms_new["MU"];
  double h_new=parms_new["H"];
  int nfreq_new=parms_new["NMATSUBARA"];
  double nsites_new=parms_new["dca.SITES"];
  
  std::cout<<"main parameters: "<<std::endl;
  std::cout<<"new beta: "<<beta_new<<"(old: "<<beta_old<<")"<<std::endl;
  std::cout<<"new mu: "<<mu_new<<"(old: "<<mu_old<<")"<<std::endl;
  std::cout<<"new h: "<<h_new<<"(old: "<<h_old<<")"<<std::endl;

  //create old and new cluster transformers
  NormalStateClusterTransformer CT_old(parms_old);
  NormalStateClusterTransformer CT_new(parms_new);
  
  //definition of Green's functions
  cluster_matsubara_kspace_gf selfenergy_orig(alps::gf::omega_k_sigma_gf(alps::gf::matsubara_positive_mesh(beta_old,nfreq_old),CT_old.momentum_mesh(), alps::gf::index_mesh(2)));
  cluster_matsubara_kspace_gf selfenergy_kspace_interp(alps::gf::omega_k_sigma_gf(alps::gf::matsubara_positive_mesh(beta_old,nfreq_old),CT_new.momentum_mesh(), alps::gf::index_mesh(2)));
  cluster_matsubara_kspace_gf selfenergy_new(alps::gf::omega_k_sigma_gf(alps::gf::matsubara_positive_mesh(beta_new,nfreq_new),CT_new.momentum_mesh(), alps::gf::index_mesh(2)));
  cluster_matsubara_kspace_gf green_new(alps::gf::omega_k_sigma_gf(alps::gf::matsubara_positive_mesh(beta_new,nfreq_new),CT_new.momentum_mesh(), alps::gf::index_mesh(2)));
  cluster_matsubara_kspace_gf green0_new(alps::gf::omega_k_sigma_gf(alps::gf::matsubara_positive_mesh(beta_new,nfreq_new),CT_new.momentum_mesh(), alps::gf::index_mesh(2)));
  
  //open and read in self-energy
  {
    alps::hdf5::archive ar(p["sigma"].as<std::string>());
    selfenergy_orig.load(ar, p["path"].as<std::string>());
    std::ofstream s_file("sigma_orig.dat"); s_file<<selfenergy_orig;
  }
  
  //interpolate cluster geometry in k-space
  std::string cluster_name_new=parms_new["dca.LATTICE"];
  std::string cluster_name_old=parms_old["dca.LATTICE"];
  if(cluster_name_new != cluster_name_old){
    throw std::invalid_argument("you are trying to interpolate different clusters. Please reimplement using gftools");
  }
  selfenergy_kspace_interp=selfenergy_orig;
  
  //interpolate for new temperature
  double nfreq_max=std::max(nfreq_old, nfreq_new);
  std::vector<double> matsubara_freqs_old(nfreq_max);
  std::vector<double> matsubara_freqs_new(nfreq_max);
  for(int i=0;i<nfreq_max;++i){
    matsubara_freqs_old[i]=(2.*i+1)*M_PI/beta_old;
    matsubara_freqs_new[i]=(2.*i+1)*M_PI/beta_new;
  }

  for(alps::gf::momentum_index j(0);j<nsites_new;++j){
    for(alps::gf::index flavor(0);flavor<2;++flavor){
      std::vector<double> sigma_real(nfreq_new);
      std::vector<double> sigma_imag(nfreq_new);
      for(alps::gf::matsubara_index i(0);i<nfreq_new;++i){
        sigma_real[i()]=selfenergy_kspace_interp(i,j,flavor).real();
        sigma_imag[i()]=selfenergy_kspace_interp(i,j,flavor).imag();
        std::cout<<"setting from: "<<std::complex<double>(sigma_real[i()], sigma_imag[i()])<<std::endl;
      }
      std::cout<<"interpolating vector of size: "<<sigma_real.size()<<std::endl;
      std::vector<double> sigma_real_interp=interpolate_vector(matsubara_freqs_old, sigma_real, matsubara_freqs_new);
      std::vector<double> sigma_imag_interp=interpolate_vector(matsubara_freqs_old, sigma_imag, matsubara_freqs_new);
      for (alps::gf::matsubara_index i(0); i<nfreq_new; ++i){
        selfenergy_new(i,j,flavor) = std::complex<double>(sigma_real_interp[i()], sigma_imag_interp[i()]);
        std::cout<<"setting to: "<<std::complex<double>(sigma_real_interp[i()], sigma_imag_interp[i()])<<std::endl;
      }
    }
  }
  //compute the self consistency: G0 from the interpolated sigma
  {
    std::ofstream s_file("sigma_in.dat"); s_file<<selfenergy_new;
  }
  std::cout<<"...running self-consistency..."<<std::endl;
  CT_new.matsubara_k_sum(green_new, selfenergy_new);
  CT_new.dyson_get_green0(green0_new, selfenergy_new, green_new);
 
  { 
    std::string filename_new=p["g0"];
    alps::hdf5::archive ar(filename_new, alps::hdf5::archive::WRITE);
    green0_new    .save(ar, "/G0");
    green_new     .save(ar, "/G");
    selfenergy_new.save(ar, "/Sigma");
  }
}

std::vector<double> solve_parabola(const std::vector<double> &xvalues, const std::vector<double> &yvalues){
  int size=3, one=1;
  Eigen::Matrix< double , 3, 3> M;
  M(0,0)=xvalues[0]*xvalues[0]; M(0,1)=xvalues[0]; M(0,2)=1;
  M(1,0)=xvalues[1]*xvalues[1]; M(1,1)=xvalues[1]; M(1,2)=1;
  M(2,0)=xvalues[2]*xvalues[2]; M(2,1)=xvalues[2]; M(2,2)=1;
  
  std::vector<double> Z(size);
  Z[0]=yvalues[0]; Z[1]=yvalues[1]; Z[2]=yvalues[2];
  
  std::vector<int> ipiv(size);
  int info=0;

  lapack::dgesv_(&size, &one, &(M(0,0)),&size, &(ipiv[0]), &(Z[0]), &size, &info);
  if(info) throw std::runtime_error("dgesv info is nonzero.");
  
  return Z;
}
double evaluate_parabola(const std::vector<double> coefficients, double xvalue){
  return coefficients[0]*xvalue*xvalue+coefficients[1]*xvalue+coefficients[2];
}
std::vector<double> interpolate_vector(const std::vector<double> &oldxvalues, const std::vector<double> &oldyvalues, const std::vector<double> &newxvalues){
  std::vector<double> newyvalues(newxvalues.size());
  
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, oldxvalues.size());
  gsl_spline_init (spline, &(oldxvalues[0]), &(oldyvalues[0]), oldxvalues.size());
  
  std::vector<double> parabola_coefficients=solve_parabola(oldxvalues, oldyvalues);
  for(int i=0;i<newxvalues.size();++i){
    if(newxvalues[i]<oldxvalues[0]){ //need to do parabolic extrapolation
      newyvalues[i]=evaluate_parabola(solve_parabola(oldxvalues, oldyvalues), newxvalues[i]);
    }else if(newxvalues[i]> oldxvalues.back()){
      newyvalues[i]=oldyvalues.back(); //we take the last value of the interpolated function, typically makes sense for AFM or SC with a real part
    }else{
      newyvalues[i]=gsl_spline_eval (spline, newxvalues[i], acc);
    }
  }
  
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  
  return newyvalues;
}
