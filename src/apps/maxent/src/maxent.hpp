/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1998-2026 ALPS Collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Legacy ALPS solver declarations from src/apps/maxent/src/maxent.hpp
// at d4533fd30c65656cf5ab7fa448dd6b921b39aa6f.

// Adapted from CQMP/Maxent src/maxent.hpp at
// e727650a1f0a427e95efebc409d907b4c84011ec.

#pragma once

#include <fstream>
#include <boost/function.hpp>
#include "maxent_matrix_def.hpp"
#include "maxent_params.hpp"
#include <random>

struct ofstream_ : std::ofstream{
    explicit ofstream_(std::streamsize precision=10){
	    this->precision(precision);
    }
};

class MaxEntHelper : private MaxEntParameters
{
public : 

  MaxEntHelper(alps::params& p);

  double omega_coord(const int i) const { return MaxEntParameters::omega_coord(i); }

  ///getter function for the discretized default model
  double Default(const int i) const { return def_[i]; }  
  ///getter function for the discretized default model
  const vector_type& Default() const { return def_; }

  //compute the 'free energy' Q (eq. 6.8 in Sebastian's thesis) according to equation D.8 in Sebastian's thesis
  double Q(const vector_type& u, const double alpha) const {
    vector_type A=transform_into_real_space(u);
    return 0.5*chi2(A)-alpha*entropy(A);
  }

  ///number of points for the Matsubara data
  int ndat() const { return MaxEntParameters::ndat(); }

  ///A->u
  vector_type transform_into_singular_space(vector_type A) const;
  ///u-A
  vector_type transform_into_real_space(vector_type u) const;
  ///A=A_i/\sigma_i; this removes \sigma_i
  vector_type get_spectrum(const vector_type& u) const;
  vector_type PrincipalValue(const vector_type &w,const vector_type &a) const;
  /// \Sigma*(V^T*RealSpace(u)*V)*\Sigma
  matrix_type left_side(const vector_type& u) const;
  /// \Sigma*U^T*(K*RealSpace(u)-y)
  vector_type right_side(const vector_type& u) const;
  ///  \delta \dot (V^T*RealSpace(u)*V)
  double step_length(const vector_type& delta, const vector_type& u) const;
  double convergence(const vector_type& u, const double alpha) const;
  double log_prob(const vector_type& u, const double alpha) const;
  double chi_scale_factor(vector_type A, const double chi_sq, const double alpha) const;
  double chi2(const vector_type& u) const;
  void print_chi2(const vector_type& u, std::ostream &os) const;
  /// compute S=\int d\omega(A(\omega)-D(\omega)-A(\omega)\ln[A(\omega)/D(\omega))
  double entropy(const vector_type& u) const;
  /// (back)continue A(\omega) to the imaginary axis; also writes to file
  void backcontinue(ofstream_ &os, const vector_type &A, const double norm, const std::string name,vector_type &ext_back);
  matrix_type constructGamma(const vector_type& A, const double alpha);
  void generateCovariantErr(const vector_type& A_in, const double alpha, ofstream_ &os);

private:
  std::mt19937::result_type bootstrap_seed_;
  ///discretized and normalized version of the default model.
  vector_type def_;
  ///check that the default model is non-zero
  void checkDefaultModel(const vector_type &D) const;
  vector_type generateGaussNoise(vector_type data, vector_type err,std::mt19937 &rng);
protected:
  bool text_output;
};




class MaxEntSimulation : private MaxEntHelper
{

public:

  ///setup of parameters
  MaxEntSimulation(alps::params parms, const std::string& outfile);
  bool run(boost::function<bool()> const& stop_callback);
  ///the maxent calculation
  void run();
  ///the evaluation and writing of files
  void evaluate();
  vector_type levenberg_marquardt(vector_type u, const double alpha) const;
  vector_type iteration(vector_type u, const double alpha, const double mu) const;
  ///define parameter defaults
  static alps::params& define_parameters(alps::params &p);

private:

  ///grid of alpha values
  vector_type alpha;
  const double norm;
  const int max_it;
  std::string name, output_file, Kernel_type;
  bool finished = false;
  bool verbose,self,make_back,gen_err;
  const int nfreq;

  vector_type lprob;
  vector_type chi_sq;
  std::vector<vector_type> spectra;
  vector_type u;
  ///averaged spectrum
  vector_type avspec;
  ///classic MaxEnt
  vector_type maxspec;
  ///historic MaxEnt
  vector_type chispec;
  ///grid of Omega points
  vector_type omegaGrid;
  ///averaged spectrum back-continued
  vector_type avspec_back;
  ///classic MaxEnt back-continued
  vector_type maxspec_back;
  ///historic MaxEnt back-continued
  vector_type chispec_back;
  ///grid of Omega points
  //posterior probability of the default model
  double postprobdef;
  ///vector of calculated Q values for each alpha iteration
  //this is one method of checking if a minimum was found
  vector_type qvec;

public:
  ///getter for avspec, the averaged spectrum
  const vector_type getAvspec() const{return avspec;}
   ///getter for chispec, the spectrum with the best chi^2
  const vector_type getChispec() const{return chispec;} 
  ///getter for maxspec, the most probable spectrum
  const vector_type getMaxspec() const{return maxspec;}
  ///getter for avspec_back, the averaged spectrum back-continued
  const vector_type getAvspecBack() const{return avspec_back;}
  ///getter for maxspec_back, the back-continued maxspec
  const vector_type getMaxspecBack() const{return maxspec_back;}
   ///getter for chispec_back, the spectrum chispec back-contonued
  const vector_type getChispecBack() const{return chispec_back;} 
  ///getter for the grid of omega points used in determining 
  ///the spectral function A(omega)
  const vector_type getOmegaGrid() const{return omegaGrid;}
  ///getter for the posterior probability of the default model
  double getPostProb() const{return postprobdef;}
  ///getter for alpha values
  const vector_type getAlphaGrid() const{return alpha;} 
  ///getter for vector of Q value per alpha iteration
  const vector_type getQvec() const{return qvec;}
}; 

///calculates the varience of a std::vector of eigen3 vectors
//note: mean,std_dev must be initialized to Zeros(nfreq())
//located in maxent_helper.cpp
void determineVariance(std::vector<vector_type> &in,vector_type &mean, vector_type &std_dev);
