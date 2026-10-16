/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1998-2018 ALPS Collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Legacy ALPS solver interface and output from src/apps/maxent/src/maxent_simulation.cpp
// at d4533fd30c65656cf5ab7fa448dd6b921b39aa6f; original notice:
/*****************************************************************************
 *
 * ALPS Project Applications
 *
 * Copyright (C) 2010 by Sebastian  Fuchs <fuchs@comp-phys.org>
 *                       Thomas Pruschke <pruschke@comp-phys.org>
 *                       Matthias Troyer <troyer@comp-phys.org>
 *               2012 by Emanuel Gull <gull@pks.mpg.de>
 *
 * ALPS Project: https://alps.comp-phys.org/
 * SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

// Adapted from CQMP/Maxent src/maxent_simulation.cpp at
// e727650a1f0a427e95efebc409d907b4c84011ec; see ../README.md.

#include "maxent.hpp"
#include <alps/hdf5/vector.hpp>
#include <cmath>
#include <Eigen/LU>
#include "eigen_hdf5.hpp"
#include <iomanip>
#include <alps/utilities/signal.hpp>

MaxEntSimulation::MaxEntSimulation(alps::params parms, const std::string& outfile)
: MaxEntHelper(define_parameters(parms))
, alpha((int)parms["N_ALPHA"])              //This is the # of \alpha parameters that should be tried.
, norm(parms["NORM"])                                             //The integral is normalized to NORM (use e.g. for self-energies
, max_it(parms["MAX_IT"])                                       //The number of iterations done in the root finding procedure
, output_file(outfile)
, Kernel_type(parms["KERNEL"].as<std::string>())
, verbose(parms["VERBOSE"])
, self(parms["SELF"])
, make_back(parms["BACKCONTINUE"])
, gen_err(parms["GENERATE_ERR"])
, nfreq(parms["NFREQ"].as<int>())
, qvec((int)parms["N_ALPHA"])
{
  // Preserve ALPS's explicit HDF5 destination and diagnostic filename prefix.
  name = outfile.size() >= 6 && outfile.substr(outfile.size() - 6) == "out.h5"
      ? outfile.substr(0, outfile.size() - 6) : outfile + ".";
  alps::signal::listen();

  if(norm != 1.) std::cerr<<"WARNING: Redefinition of parameter NORM: Input (and output) data are assumed to be normalized to NORM."<<std::endl;
  const double alpha_min = parms["ALPHA_MIN"];                                          //Smallest value of \alpha that is tried
  const double alpha_max = parms["ALPHA_MAX"];                                          //Largest  value of \alpha that is tried
  alpha[0] = alpha_max;
  for (Eigen::Index a=1; a<alpha.size(); ++a)                                            //These are all the alpa values on a log grid
    alpha[a] =  alpha[a-1] * std::pow(alpha_min/alpha_max, 1./double(alpha.size()-1));
}
namespace {
// Whether a parameter was given: supplied as text, or assigned a value.
bool given(const alps::params& p, const std::string& name) {
  return p.supplied(name) || p.exists(name);
}
// Upstream defines these without a default, since only some inputs and
// default models read them. ALPS rejects missing parameters before a run, so
// they are defined only when given.
void define_optional(alps::params& p, const std::string& name, const std::string& descr) {
  if (given(p, name)) p.define<double>(name, descr);
}
}

///define parameter defaults
alps::params& MaxEntSimulation::define_parameters(alps::params &p){
  p.description("Maxent - a utility for " 
    "performing analytic continuation \n \t using the method of Maximum Entropy\n");
  //---------------------------------
  //    General
  //---------------------------------
  p.define<bool>("DATA_IN_HDF5",false,"true if data is in HDF5 format");
  // ALPS: text output and backcontinuation stay opt-in, as in its earlier solver.
  p.define<bool>("TEXT_OUTPUT",false,"true if results should be output to text files");
  p.define<bool>("VERBOSE",false,"true to print verbose output");
  p.define<bool>("SELF",false,"input is a self energy");
  p.define<int>("MAX_IT",1000,"Maximum Iterations for the fitting routine");
  p.define<int>("N_ALPHA",60,"Number of alpha samples");
  p.define<double>("ALPHA_MIN",0.01,"Minimum alpha");
  p.define<double>("ALPHA_MAX",20,"Maximum alpha");
  p.define<double>("NORM",1.0,"NORM");
  p.define<bool>("BACKCONTINUE",false,"Output A(omega) back to imaginary axis");
  p.define<bool>("GENERATE_ERR",false,"Generate a bootstrap approximation for error bars");
  // ALPS: mcbase's SEED default of 42 instead of upstream's 0.
  p.define<unsigned int>("SEED",42,"Seed for reproducible bootstrap error bars");
  //*********************************
  // ALPS: T may be given instead of BETA.
  if (given(p, "T")) p.define<double>("T","temperature; BETA = 1/T unless BETA is given");
  if (given(p, "T") && !given(p, "BETA"))
    p.define<double>("BETA",1./p["T"].as<double>(),"beta, inverse temperature");
  else
    p.define<double>("BETA","beta, inverse temperature");
  p.define<int>("NDAT","# of input points");
  p.define<std::string>("DATA","","data file input");
  p.define<std::string>("COVARIANCE_MATRIX","","name of covariance matrix file");
  p.define<double>("CM_EIGENVALUE_CUTOFF",1e-10, "Cutoff for eigenvalues of the covariance matrix");
  p.define<std::string>("BASENAME","","Specified output name (generated if not given)");
  // ALPS: no MODEL_RUNS; the entry point runs one default model.
  define_optional(p, "X_0","G input for param file entry");
  define_optional(p, "SIGMA_0","G error input for param file entry");
  define_optional(p, "TAU_0","Used for input tau points");
  // ALPS: the run's time limit and its normalization constraint.
  p.define<int>("MAX_TIME",60,"time limit in seconds");
  p.define<bool>("ENFORCE_NORMALIZATION",false,"add the normalization as a data point with error SIGMA_NORMALIZATION");
  if (p["ENFORCE_NORMALIZATION"].as<bool>())
    p.define<double>("SIGMA_NORMALIZATION","error of the normalization data point");


  //---------------------------------
  //      Default Model
  //---------------------------------
  p.define<double>("OMEGA_MAX",10,"Maximum frequency for A(omega) grid");
  define_optional(p, "OMEGA_MIN","Minimum frequency, or =-OMEGA_MAX");
  p.define<std::string>("DEFAULT_MODEL","flat","Default model for entropy");
  p.define<double>("NORM1",0.5,"for Two Gaussians model");
  p.define<double>("SHIFT",0.0,"shift of a model");
  p.define<double>("SHIFT1",0.0,"for Two Gaussians model");
  define_optional(p, "SHIFT2","for Two Gaussians model");
  define_optional(p, "SIGMA1","for Two Gaussians model");
  define_optional(p, "SIGMA2","for Two Gaussians model");
  define_optional(p, "SIGMA","stddev - For Gaussian models");
  define_optional(p, "GAMMA","width of Lorentzian model"); 
  define_optional(p, "GAMMA1","for Two Lorentzian models");
  define_optional(p, "GAMMA2","for Two Lorentzian models");

  define_optional(p, "LAMBDA","for ___ExpDecay models");
  define_optional(p, "BOSE_NORM","General Double Gaussian Norm");
  
  //---------------------------------
  //    Grid
  //---------------------------------
  p.define<double>("CUT",0.01,"cut for lorentzian grids");
  p.define<double>("SPREAD",4,"spread for quadratic grid");
  p.define<double>("LOG_MIN",1.0e-4,"log_min for log grid");
  p.define<int>("NFREQ",1000,"Number of A(omega) frequencies");
  // ALPS: the default grid stays linear, as in its earlier solver.
  p.define<std::string>("FREQUENCY_GRID","linear","Type of frequency grid");

  //---------------------------------
  //    Kernel
  //---------------------------------
  p.define<std::string>("DATASPACE","time","Time or Frequency space");
  p.define<std::string>("KERNEL","fermionic","Type of kernel: Fermionic,Bosonic,TZero,Legendre"); 
  p.define<bool>("PARTICLE_HOLE_SYMMETRY",false,"Set =true if particle hole symmetric"); 

  // In ALPS HDF5 input, supplying this parameter selects /Covariance even
  // when its value is empty; an absent parameter selects /Error instead. Only
  // supplied text can tell the two apart once the default is in place, and
  // every call re-derives the marker, since define() re-reads that text.
  if (p["DATA_IN_HDF5"].as<bool>() && p.supplied("COVARIANCE_MATRIX")
      && p["COVARIANCE_MATRIX"].as<std::string>().empty())
    p["COVARIANCE_MATRIX"] = std::string("HDF5");
  // The upstream backcontinuation uses the error-scaled kernel. After a
  // covariance rotation its rows no longer label the original input grid.
  if (p["BACKCONTINUE"].as<bool>() && !p["COVARIANCE_MATRIX"].as<std::string>().empty())
    throw std::invalid_argument("BACKCONTINUE with COVARIANCE_MATRIX is not supported");
  return p;
}

// Preserve ALPS's stop-callback ordering around a complete alpha sweep.
bool MaxEntSimulation::run(boost::function<bool()> const& stop_callback) {
  bool stopped = false;
  while (!(stopped = stop_callback()) && !finished) {
    run();
    evaluate();
    finished = true;
  }
  return !stopped;
}

void MaxEntSimulation::run()
{
  lprob.resize(alpha.size());
  chi_sq.resize(alpha.size());
  spectra.resize(alpha.size());
  u = transform_into_singular_space(Default());

  ofstream_ spectral_function_file;
  ofstream_ fits_file;

  if (text_output) {
    spectral_function_file.open((name+"spex.dat").c_str());
    fits_file.open((name+"fits.dat").c_str());
  }
  //this loop is the 'core' of the maxent program: iterate over all alphas, compute the spectra, normalization, and probabilities
  //loop over all alpha values
  for (Eigen::Index a=0; a<alpha.size(); ++a) {
    std::cerr << "alpha it: " << a << "\t";
    //fitting procedure for 'u'
    u = levenberg_marquardt(u, alpha[a]);
    //computation of spectral function out of 'u'
    vector_type A = get_spectrum(u);
    //computation of normalization
    std::cerr << "norm: " << transform_into_real_space(u).sum() << "\t";
    if (text_output) {
      spectral_function_file<<"# alpha: "<<alpha[a]<<std::endl;
      for (Eigen::Index i=0; i<A.size(); ++i)
        spectral_function_file << omega_coord(i) << " " << A[i] << "\n";
      spectral_function_file << "\n";
    }
    //computation of probability
    lprob[a] = log_prob(u, alpha[a]);
    spectra[a] = A;
    //computation of chi2
    double chi_squared = chi2(transform_into_real_space(u));
    chi_sq[a] = chi_squared;
    if (verbose) std::cerr << "0.5*chi2  : " << 0.5*chi_squared;
    std::cerr << std::endl;
    if (text_output) print_chi2(transform_into_real_space(u), fits_file);
    qvec(a)=Q(u,alpha[a]); 
  }
    omegaGrid.resize(nfreq);
    for(Eigen::Index i=0;i<omegaGrid.size();i++)
	    omegaGrid(i)=omega_coord(i);
}
  //everything from here on down is evaluation.
void MaxEntSimulation::evaluate(){
  if (text_output) {
    ofstream_ chi_squared_file;
    chi_squared_file.open((name+"chi2.dat").c_str());
    for (Eigen::Index a=0; a<chi_sq.size(); ++a){
      chi_squared_file << alpha[a] << " " << chi_sq[a] << std::endl;
    }
  }
  int a_chi = 0;
  double diff = std::abs(chi_sq[0]-ndat());
  for (Eigen::Index a=1; a<chi_sq.size(); ++a) {
    double diff_new = std::abs(chi_sq[a]-ndat());
    if (diff_new < diff) {
      diff = diff_new;
      a_chi = a;
    }
  }

  vector_type def = get_spectrum(transform_into_singular_space(Default()));
  chispec = spectra[a_chi]*norm;
  if (text_output){
    ofstream_ chispec_file;
    chispec_file.open((name+"chispec.dat").c_str());
    for (Eigen::Index i=0; i<spectra[0].size(); ++i){
      chispec_file << omega_coord(i) << " " << spectra[a_chi][i]*norm << " " << def[i]*norm << std::endl;
    }
  }
  //const int max_a = max_lprob-lprob.begin();
  int max_a,nothing; double max_lprob;
  max_lprob=lprob.maxCoeff(&max_a,&nothing);
  const double factor = chi_scale_factor(spectra[max_a], chi_sq[max_a], alpha[max_a]);
  if (verbose) std::cerr << "chi scale factor: " << factor << std::endl;

	alps::hdf5::archive ar(output_file, "w");
	ar << alps::make_pvp("/alpha/values",alpha);

  vector_type om(spectra[0].size());
  for (Eigen::Index i=0;i<om.size();i++) om[i] = omega_coord(i);
  ar<<alps::make_pvp("/spectrum/omega",om);

  //output 'maximum' spectral function (classical maxent metod)
  if (text_output){
    ofstream_ maxspec_file;
    maxspec_file.open((name+"maxspec.dat").c_str());
    for (Eigen::Index i=0; i<spectra[0].size(); ++i)
      maxspec_file << omega_coord(i) << " " << spectra[max_a][i]*norm << " " << def[i]*norm << std::endl;
  }
	
  maxspec = spectra[max_a]*norm;
  vector_type specchi = spectra[a_chi]*norm;
	ar << alps::make_pvp("/spectrum/chi",specchi);
  ar << alps::make_pvp("/spectrum/maximum",maxspec);
  
	vector_type prob(lprob.size());
  for (Eigen::Index a=0; a<prob.size(); ++a)
    prob[a] = exp(lprob[a]-max_lprob);
  double probnorm = 0;
  for (Eigen::Index a=0; a<prob.size()-1; ++a)
    probnorm += 0.5*(prob[a]+prob[a+1])*(alpha[a]-alpha[a+1]);
  prob /= probnorm;
  ar << alps::make_pvp("/alpha/probability",prob);
  if (text_output){
    ofstream_ prob_str;
    prob_str.open((name+"prob.dat").c_str());
    for (Eigen::Index a=0; a<prob.size(); ++a) {
      prob_str << alpha[a] << "\t" << prob[a] << "\n";
    }
  }
  postprobdef = 0;
  for (Eigen::Index a=0; a<lprob.size()-1; ++a)
    postprobdef += 0.5*(exp(lprob[a])+exp(lprob[a+1]))*(alpha[a]-alpha[a+1]);
  std::cout << "posterior probability of the default model: " << postprobdef << std::endl;

  //compute 'average' spectral function (Brian's method)
  avspec.resize(spectra[0].size());
  for (Eigen::Index i=0; i<avspec.size(); ++i) {
    avspec[i] = 0.;
    for (Eigen::Index a=0; a<prob.size()-1; ++a)
      avspec[i] += 0.5*(prob[a]*spectra[a][i] +prob[a+1]*spectra[a+1][i])*(alpha[a]-alpha[a+1]);
  }
  //Estimate the variance for the spectrum
  vector_type varspec(spectra[0].size());
  for (Eigen::Index i=0; i<varspec.size(); ++i) {
    varspec[i] = 0.;
    for (Eigen::Index a=0; a<prob.size()-1; ++a)
      varspec[i] += 0.5*(prob[a]*(spectra[a][i]-avspec[i])*(spectra[a][i]-avspec[i]) + prob[a+1]*(spectra[a+1][i]-avspec[i])*(spectra[a+1][i]-avspec[i]))*(alpha[a]-alpha[a+1]);
  }
  avspec *= norm;
  varspec *= norm*norm;

  if (text_output){
    ofstream_ avspec_file;
    avspec_file.open((name+"avspec.dat").c_str());
    for (Eigen::Index i=0; i<avspec.size(); ++i)
      avspec_file << omega_coord(i) << " " << avspec[i] << " " << def[i]*norm << std::endl;
  }
  ar << alps::make_pvp("/spectrum/average",avspec);
  ar << alps::make_pvp("/spectrum/variance",varspec);
  
  if(gen_err){
    //Bootstrap errors
    ofstream_ boot_file;
    boot_file.open((name+"booterr.dat").c_str());
    generateCovariantErr(maxspec,alpha[max_a],boot_file);
  }
  

  if(Kernel_type=="anomalous"){ //for the anomalous function: use A(omega)=Im Sigma(omega)/(pi omega).
    ofstream_ maxspec_anom_str;maxspec_anom_str.open((name+"maxspec_anom.dat").c_str());
    ofstream_ avspec_anom_str; avspec_anom_str.open((name+"avspec_anom.dat").c_str());
    vector_type spec(avspec.size());
    for (Eigen::Index i=0; i<avspec.size(); ++i){
      //if(omega_coord(i)>=0.)
      spec[i] = avspec[i]*omega_coord(i)*M_PI;
      avspec_anom_str << omega_coord(i) << " " << spec[i]<<std::endl;
    }
    ar << alps::make_pvp("/spectrum/anomalous/average",spec);
    for (Eigen::Index i=0; i<spectra[0].size(); ++i){
      //if(omega_coord(i)>=0.)
      spec[i] = spectra[max_a][i]*norm*omega_coord(i)*M_PI;
      maxspec_anom_str << omega_coord(i) << " " << spec[i] << std::endl;
    }
    ar << alps::make_pvp("/spectrum/anomalous/maximum",spec);
  }
  if(Kernel_type=="bosonic"){ //for the anomalous function: use A(Omega_)=Im chi(Omega_)/(pi Omega_) (as for anomalous)
    vector_type spec(avspec.size());
    for (Eigen::Index i=0; i<avspec.size(); ++i){
      // ALPS reports Im chi = pi * omega * A.
      spec[i] = avspec[i]*omega_coord(i)*M_PI;
    }
    if (text_output) {
      ofstream_ avspec_anom_str;avspec_anom_str.open((name+"maxspec_bose.dat").c_str());
      for (Eigen::Index i=0; i<avspec.size(); ++i){
        //if(omega_coord(i)>=0.)
        avspec_anom_str << omega_coord(i) << " " << spec[i]<<std::endl;
      }
    }
    ar << alps::make_pvp("/spectrum/bosonic/average",spec);
    for (Eigen::Index i=0; i<spectra[0].size(); ++i){
      //if(omega_coord(i)>=0.)
      spec[i] = spectra[max_a][i]*norm*omega_coord(i)*M_PI;
    }
    if (text_output) {
      ofstream_ maxspec_anom_str;maxspec_anom_str.open((name+"avspec_bose.dat").c_str());
      for (Eigen::Index i=0; i<spectra[0].size(); ++i){
        maxspec_anom_str << omega_coord(i) << " " << spec[i] << std::endl;
      }
    }
    ar << alps::make_pvp("/spectrum/bosonic/maximum",spec);
  }

  //don't understand why this was commented out...
  if(self){
    // A quick word about normalization here. Usually we have G(iomega_n) = -1/pi \int_{-\infty}^\infty Im G(omega)/(omega_n - omega).
    // However, we are not interested in Im G but instead in A. In the case of the self-energy we have, analogously,
    // Sigma(i\omega_n) = -1/pi \int_{-\infty}^\infty Im \Sigma(omega)/(omega_n - omega); and we define A_\Sigma(omega) = -1/pi Sigma(omega). This makes
    // A_\Sigma be always positive, whereas Im Sigma(omega) is always negative.
    // here we compute Im Sigma out of A:
    //
    // for the self energy: use Im Sigma(omega)=-A(omega)*pi
    ofstream_ maxspec_self_str;maxspec_self_str.open((name+"maxspec_self.dat").c_str());
    ofstream_ avspec_self_str; avspec_self_str.open((name+"avspec_self.dat").c_str());
    ofstream_ chispec_self_str; chispec_self_str.open((name+"chispec_self.dat").c_str());
    for (Eigen::Index i=0; i<avspec.size(); ++i){
      avspec_self_str << omega_coord(i) << " " << -avspec[i]*M_PI<<std::endl;
    }
    for (Eigen::Index i=0; i<spectra[0].size(); ++i){
      maxspec_self_str << omega_coord(i) << " " << -spectra[max_a][i]*norm*M_PI << std::endl;
    }
    for (Eigen::Index i=0; i<specchi.size(); ++i){
      chispec_self_str << omega_coord(i) << " " << -specchi[i]*M_PI<< " " << -def[i]*norm*M_PI<<std::endl;
    }
    //for public facing variables
    avspec*=-M_PI;
    maxspec*=-M_PI;
    specchi*=-M_PI;
  }

  if(make_back){
    ofstream_ avspec_back_file,maxspec_back_file,chispec_back_file;

    if(text_output){
      avspec_back_file.open((name+"avspec_back.dat").c_str());
      maxspec_back_file.open((name+"maxspec_back.dat").c_str());
      chispec_back_file.open((name+"chispec_back.dat").c_str());
    }
    const std::string sp = "    ";
    double norm_fix = norm;
    //fix self-energy normalization 
    if(self){
      norm_fix *=-M_PI;
    }

    std::cerr << "spectra"<<sp<< " max backcont diff" <<sp<<  "chi^2 value " <<std::endl;
    std::cerr << "======="<< sp<<" ================="<< sp<< "=========== " <<std::endl;
    backcontinue(chispec_back_file,specchi,norm_fix,"chispec",chispec_back);
    backcontinue(avspec_back_file,avspec,norm_fix,"avspec ",avspec_back);
    backcontinue(maxspec_back_file,maxspec,norm_fix,"maxspec",maxspec_back);
  }
  ar.close();
}



//this is the levenberg marquardt fitting procedure. It minimizes the quantity Q = 1/2 chi^2 - \alpha S
// 
vector_type MaxEntSimulation::levenberg_marquardt(vector_type u, const double alpha) const
{

  double mu = 1e-18;
  const double nu = 1.3;
  double Q1=0.;
  int it = 0;
  int it2 = 0;
  for (; it<max_it; it++) {
    vector_type delta;
    if(std::isnan(Q1))
        throw std::logic_error("Q=NaN, something went wrong");
    for (it2=0; it2<max_it; ++it2) {
      //compute change vector delta to u
      delta = iteration(u, alpha, mu);
      //compute Q = 1/2 chi^2 - \alpha S
      Q1 = Q(u+delta, alpha);
      if (step_length(delta, u)<=0.02) {
        break;
      }
      else if (mu<1e20) {
        mu *= nu;
      }

    } 
    u += delta;
    if (convergence(u, alpha)<=1e-4)
      break;
  }
  if (it == max_it) std::cerr<<"WARNING: iteration reached max_it without converging, your minimizer is having problems. Please be careful!"<<std::endl;
  if (verbose) std::cerr <<"Iterations: " << it+1 << "\t";
  std::cerr << "Q = 0.5chi^2-\\alpha*entropy: " << Q1 << "\t";
  if (verbose) std::cerr << "entropy: "<<entropy(transform_into_real_space(u))<<"\talpha*entropy: "<<alpha*entropy(transform_into_real_space(u))<<"\t ";
  return u;
}



//this function computes the change delta to the vector 'u' 
//to be used in the Levenberg Marquardt fitting procedure
vector_type MaxEntSimulation::iteration(vector_type u, const double alpha, const double mu) const
{
  matrix_type M = left_side(u);
  for (Eigen::Index i=0; i<M.rows(); ++i)
    M(i,i) += alpha + mu;
  vector_type b = right_side(u) + alpha*u;
  matrix_type B(b.size(),1);
  for (Eigen::Index i=0; i<M.rows(); ++i)
    B(i,0) = -b[i];
  //NOTE: gesv uses LU decomp, but we can switch to a safe QR routine as well
  matrix_type Bp = M.lu().solve(B);
  //may need a transposeInPlace();
  return Bp;
}
