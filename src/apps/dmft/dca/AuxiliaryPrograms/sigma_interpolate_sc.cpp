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

#include <Eigen/Dense>

void define_parameters_commandline(alps::params &p){
  p.define<std::string>("pfile_old", "old parameter file");
  p.define<std::string>("pfile_new", "new parameter file");
  p.define<std::string>("sigma", "self energy file (in GFTools hdf5 format)");
  p.define<std::string>("path", "/Sigma","hdf5 path to self energy");
  p.define<std::string>("g0", "hdf5 path to the G0 obtained with the interpolated self-energy");
  p.define<double>("mult_anomalous", 1., "factor applied to the interpolated anomalous self-energy");
}
void define_dca_parameters(alps::params &p);
void define_selfconsistency_parameters(alps::params &p);
void define_parameters_pfile(alps::params &p){
  define_selfconsistency_parameters(p);
  define_dca_parameters(p);
}

std::vector<double> solve_parabola(const std::vector<double> &xvalues, const std::vector<double> &yvalues){
  Eigen::Matrix3d M;
  M << xvalues[0]*xvalues[0], xvalues[0], 1,
       xvalues[1]*xvalues[1], xvalues[1], 1,
       xvalues[2]*xvalues[2], xvalues[2], 1;
  const Eigen::Vector3d Z(yvalues[0], yvalues[1], yvalues[2]);

  const Eigen::FullPivLU<Eigen::Matrix3d> lu(M);
  if(!lu.isInvertible()) throw std::runtime_error("solve_parabola: the three x values are not distinct.");
  const Eigen::Vector3d coefficients = lu.solve(Z);
  return std::vector<double>(coefficients.data(), coefficients.data()+3);
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

void run(int argc, char **argv){

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

  const double mult_anomalous=p["mult_anomalous"];

  //general error check
  double beta_old=parms_old["BETA"];
  double mu_old=parms_old["MU"];
  double h_old=parms_old["dca.H"];
  double eta1_old=parms_old["dca.ETA1"];
  double eta2_old=parms_old["dca.ETA2"];
  double nsites_old=parms_old["dca.SITES"];
  int nfreq_old = parms_old["NMATSUBARA"];

  double beta_new=parms_new["BETA"];
  double mu_new=parms_new["MU"];
  double h_new=parms_new["dca.H_INIT"];
  double eta1_new=parms_new["dca.ETA1_INIT"];
  double eta2_new=parms_new["dca.ETA2_INIT"];
  int nfreq_new=parms_new["NMATSUBARA"];
  double nsites_new=parms_new["dca.SITES"];
  
  std::cout<<"adjusting parameters: "<<std::endl;
  std::cout<<"new beta: "<<beta_new<<"(old: "<<beta_old<<")"<<std::endl;
  std::cout<<"new mu: "<<mu_new<<"(old: "<<mu_old<<")"<<std::endl;
  std::cout<<"new eta: "<<eta1_new<<"(old: "<<eta1_old<<")"<<std::endl;

  //create old and new cluster transformers
  NambuClusterTransformer CT_old(parms_old);
  NambuClusterTransformer CT_new(parms_new);

  double_cell_matsubara_kspace_gf selfenergy_orig(
    alps::gf::omega_k_sigma1_sigma2_gf(
      alps::gf::matsubara_positive_mesh(beta_old,nfreq_old), CT_old.momentum_mesh(), alps::gf::index_mesh(2), alps::gf::index_mesh(2)));

  double_cell_matsubara_kspace_gf selfenergy_kspace_interp(
    alps::gf::omega_k_sigma1_sigma2_gf(
      alps::gf::matsubara_positive_mesh(beta_old,nfreq_old), CT_new.momentum_mesh(), alps::gf::index_mesh(2), alps::gf::index_mesh(2)));

  double_cell_matsubara_kspace_gf selfenergy_new(
    alps::gf::omega_k_sigma1_sigma2_gf(
      alps::gf::matsubara_positive_mesh(beta_new,nfreq_new), CT_new.momentum_mesh(), alps::gf::index_mesh(2), alps::gf::index_mesh(2)));

  double_cell_matsubara_kspace_gf green_new(
    alps::gf::omega_k_sigma1_sigma2_gf(
      alps::gf::matsubara_positive_mesh(beta_new,nfreq_new), CT_new.momentum_mesh(), alps::gf::index_mesh(2), alps::gf::index_mesh(2)));

  double_cell_matsubara_kspace_gf green0_new(
    alps::gf::omega_k_sigma1_sigma2_gf(
      alps::gf::matsubara_positive_mesh(beta_new,nfreq_new), CT_new.momentum_mesh(), alps::gf::index_mesh(2), alps::gf::index_mesh(2)));
  //cluster_nambu_itime_kspace_gf green0_itime(ntau, nsites_new, 2, diagonal);

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
  // The self-energy is given on the old frequencies and wanted on the new ones.
  if(selfenergy_orig.mesh1().extent()!=nfreq_old)
    throw std::invalid_argument("the self-energy has "+std::to_string(selfenergy_orig.mesh1().extent())
                                +" Matsubara frequencies, but NMATSUBARA of the old parameter file is "
                                +std::to_string(nfreq_old));
  if(std::abs(selfenergy_orig.mesh1().beta()-beta_old)>1e-10*beta_old)
    throw std::invalid_argument("the self-energy is given at beta="+std::to_string(selfenergy_orig.mesh1().beta())
                                +", but BETA of the old parameter file is "+std::to_string(beta_old));
  std::vector<double> matsubara_freqs_old(nfreq_old);
  std::vector<double> matsubara_freqs_new(nfreq_new);
  for(int i=0;i<nfreq_old;++i) matsubara_freqs_old[i]=(2.*i+1)*M_PI/beta_old;
  for(int i=0;i<nfreq_new;++i) matsubara_freqs_new[i]=(2.*i+1)*M_PI/beta_new;
  
  std::vector<double> sigma_real(nfreq_old);
  std::vector<double> sigma_imag(nfreq_old);
  for(alps::gf::momentum_index j(0);j<nsites_new;++j){
    for(alps::gf::index p(0);p<2;++p){
      for(alps::gf::index q(0);q<2;++q){
        for(alps::gf::matsubara_index i(0);i<nfreq_old;++i){
          sigma_real[i()]=selfenergy_kspace_interp(i,j,p,q).real();
          sigma_imag[i()]=selfenergy_kspace_interp(i,j,p,q).imag();
        }

        std::cout<<"interpolating vector of size: "<<sigma_real.size()<<std::endl;
        std::vector<double> sigma_real_interp=interpolate_vector(matsubara_freqs_old, sigma_real, matsubara_freqs_new);
        std::vector<double> sigma_imag_interp=interpolate_vector(matsubara_freqs_old, sigma_imag, matsubara_freqs_new);
        for (alps::gf::matsubara_index i(0); i<nfreq_new; ++i){
          selfenergy_new(i,j,p,q) = std::complex<double>(sigma_real_interp[i()], sigma_imag_interp[i()]);
          //std::cout<<"setting to: "<<std::complex<double>(sigma_real_interp[i()], sigma_imag_interp[i()])<<std::endl;
        }

        if(p() != q() && mult_anomalous!=1){
          for(alps::gf::matsubara_index i(0); i<nfreq_new; ++i){
            selfenergy_new(i,j,p,q)*=mult_anomalous;
          }
        }
      }
    }
  }

  {
    std::ofstream s_file("sigma_in.dat"); s_file<<selfenergy_new;
  }
  std::cout<<"...running self-consistency..."<<std::endl;
  CT_new.matsubara_k_sum(green_new, selfenergy_new, eta1_new, eta2_new, h_new);
  CT_new.dyson_get_green0(green0_new, selfenergy_new, green_new);

  {
    std::string filename_new=p["g0"];
    alps::hdf5::archive ar(filename_new, alps::hdf5::archive::WRITE);
    green0_new    .save(ar, "/G0");
    green_new     .save(ar, "/G");
    selfenergy_new.save(ar, "/Sigma");
  }
}
int main(int argc, char **argv){
  try {
    run(argc, argv);
    return 0;
  } catch (const std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    return 1;
  }
}
