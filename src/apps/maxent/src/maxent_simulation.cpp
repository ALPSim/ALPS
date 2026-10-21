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

#include <utility> // Boost.Math 1.76 includes this inside a namespace.
#include <boost/math/constants/constants.hpp>
#include "maxent.hpp"
#include <alps/config.h> // needed to set up correct bindings
#include <alps/hdf5/ublas/vector.hpp>
#include <alps/utilities/signal.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/bindings/lapack/driver/gesv.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas.hpp>





void MaxEntSimulation::define_parameters(alps::params &p)
{
  //data
  p.define<int>("NDAT", "number of data points");
  p.define<double>("NORM", "normalization of the spectral function; the data are divided by it");
  p.define<std::string>("DATA", "", "file with the data; without it, the data are the parameters X_i and SIGMA_i");
  p.define<bool>("DATA_IN_HDF5", false, "read DATA (and COVARIANCE_MATRIX) from an HDF5 file");
  p.define<std::string>("COVARIANCE_MATRIX", "", "file with the covariance matrix of the data; replaces SIGMA_i");
  p.define<std::string>("DATASPACE", "time", "space of the data: time, frequency or legendre");
  p.define<std::string>("KERNEL", "fermionic", "kernel: fermionic, bosonic, anomalous or Boris");
  p.define<bool>("PARTICLE_HOLE_SYMMETRY", false, "the data are particle-hole symmetric");
  p.define<bool>("ENFORCE_NORMALIZATION", false, "add the normalization as a data point with error SIGMA_NORMALIZATION");
  if(p.exists("ENFORCE_NORMALIZATION") && p["ENFORCE_NORMALIZATION"].as<bool>())
    p.define<double>("SIGMA_NORMALIZATION", "error of the normalization data point");
  //data given as parameters, and the imaginary times of the data points
  if(p.exists("NDAT")) {
    for(int i=0; i<p["NDAT"].as<int>(); ++i) {
      std::string const index=boost::lexical_cast<std::string>(i);
      if(p.supplied("X_"+index)) p.define<double>("X_"+index, "data point i");
      if(p.supplied("SIGMA_"+index)) p.define<double>("SIGMA_"+index, "error of data point i");
      if(p.supplied("TAU_"+index)) p.define<double>("TAU_"+index, "imaginary time of data point i");
    }
  }
  //temperature: T, or 1/BETA
  if(p.supplied("T")) p.define<double>("T", "temperature");
  else {
    p.define<double>("BETA", "inverse temperature");
    if(p.exists("BETA")) p.define<double>("T", 1./p["BETA"].as<double>(), "temperature; 1/BETA by default");
  }
  //real-frequency grid
  p.define<int>("NFREQ", "number of real frequencies");
  p.define<double>("OMEGA_MAX", "upper end of the frequency grid");
  if(p.exists("OMEGA_MAX")) p.define<double>("OMEGA_MIN", -p["OMEGA_MAX"].as<double>(), "lower end of the frequency grid; -OMEGA_MAX by default");
  p.define<std::string>("FREQUENCY_GRID", "", "frequency grid: Lorentzian, half Lorentzian, quadratic, log, linear or lin_log");
  p.define<double>("CUT", 0.01, "cutoff of the Lorentzian grids");
  p.define<double>("SPREAD", 4., "spread of the quadratic grid");
  p.define<double>("LOG_MIN", 1.0e-4, "smallest frequency of the logarithmic grid");
  //default model, and the parameters of the chosen one
  p.define<std::string>("DEFAULT_MODEL", "flat", "default model: flat, gaussian, twogaussians, shifted gaussian, double gaussian, general double gaussian, linear rise exp decay, quadratic rise exp decay, or a file name");
  p.define<double>("BLOW_UP", 1., "factor by which the default model is widened");
  p.define<double>("SHIFT1", 0., "position of the first Gaussian of twogaussians");
  p.define<double>("NORM1", 0.5, "weight of the first Gaussian of twogaussians");
  std::string const model=p["DEFAULT_MODEL"].as<std::string>();
  bool const gaussian=(model=="gaussian" || model=="shifted gaussian" || model=="double gaussian" || model=="general double gaussian");
  if(gaussian) p.define<double>("SIGMA", "width of the Gaussian default model");
  if(gaussian && model!="gaussian") p.define<double>("SHIFT", "shift of the Gaussian default model");
  if(model=="general double gaussian") p.define<double>("BOSE_NORM", "bosonic weight of the general double Gaussian default model");
  if(model=="twogaussians") {
    p.define<double>("SIGMA1", "width of the first Gaussian of twogaussians");
    p.define<double>("SIGMA2", "width of the second Gaussian of twogaussians");
    p.define<double>("SHIFT2", "position of the second Gaussian of twogaussians");
  }
  if(model=="linear rise exp decay" || model=="quadratic rise exp decay")
    p.define<double>("LAMBDA", "decay rate of the default model");
  //iteration and output
  p.define<int>("N_ALPHA", "number of values of alpha");
  p.define<double>("ALPHA_MIN", "smallest alpha");
  p.define<double>("ALPHA_MAX", "largest alpha");
  p.define<int>("MAX_IT", 1000, "number of iterations of the root finding");
  p.define<int>("MAX_TIME", 60, "time limit in seconds");
  p.define<int>("SEED", 42, "random seed; accepted for compatibility, MaxEnt draws no random numbers");
  p.define<bool>("VERBOSE", false, "print extra information");
  p.define<bool>("TEXT_OUTPUT", false, "also write results as text files");
  p.define<bool>("SELF", false, "the data are a self-energy");
}

MaxEntSimulation::MaxEntSimulation(const alps::params &parms,const std::string &outfile)
: MaxEntHelper(parms)
, alpha((int)parms["N_ALPHA"])              //This is the # of \alpha parameters that should be tried.
, norm(parms["NORM"])                                                  //The integral is normalized to NORM (use e.g. for self-energies
, max_it(parms["MAX_IT"])                                       //The number of iterations done in the root finding procedure
, name(outfile,0,outfile.size()-6)
, Kernel_type(parms["KERNEL"].as<std::string>())
, finished(false)
, verbose(parms["VERBOSE"])
, text_output(parms["TEXT_OUTPUT"])
, self(parms["SELF"])
{
  // Preserve the seed conversion and archive signal setup previously supplied
  // by mcbase, even though deterministic MaxEnt does not consume random values.
  (void)static_cast<int>(parms["SEED"]);
  alps::signal::listen();
  if(norm != 1.) std::cerr<<"WARNING: Redefinition of parameter NORM: Input (and output) data are assumed to be normalized to NORM."<<std::endl;
  const double alpha_min = parms["ALPHA_MIN"];                                          //Smallest value of \alpha that is tried
  const double alpha_max = parms["ALPHA_MAX"];                                          //Largest  value of \alpha that is tried
  alpha[0] = alpha_max;
  for (std::size_t a=1; a<alpha.size(); ++a)                                            //These are all the alpa values on a log grid
    alpha[a] =  alpha[a-1] * std::pow(alpha_min/alpha_max, 1./double(alpha.size()-1));
}


MaxEntSimulation::~MaxEntSimulation() 
{
}

bool MaxEntSimulation::run(boost::function<bool()> const& stop_callback)
{
  // Keep mcbase's callback ordering: check before a step and once after it,
  // even when the step finished. One step computes the complete alpha sweep.
  bool stopped = false;
  while (!(stopped = stop_callback()) && !finished)
    dostep();
  return !stopped;
}



void MaxEntSimulation::dostep() 
{
  if (finished)
    return;
  vector_type lprob(alpha.size());
  vector_type chi_sq(alpha.size());
  std::vector<vector_type> spectra(alpha.size());
  vector_type u = transform_into_singular_space(Default());

  if (text_output) {
    spex_str.open((name+"spex.dat").c_str());
    chisq_str.open((name+"chi2.dat").c_str());
    avspec_str.open((name+"avspec.dat").c_str());
    maxspec_str.open((name+"maxspec.dat").c_str());
    chispec_str.open((name+"chispec.dat").c_str());
    fits_str.open((name+"fits.dat").c_str());
    prob_str.open((name+"prob.dat").c_str());
  }
  //this loop is the 'core' of the maxent program: iterate over all alphas, compute the spectra, normalization, and probabilities
  //loop over all alpha values
  for (std::size_t a=0; a<alpha.size(); ++a) {
    std::cerr << "alpha it: " << a << "\t";
    //fitting procedure for 'u'
    u = levenberg_marquardt(u, alpha[a]);
    //computation of spectral function out of 'u'
    vector_type A = get_spectrum(u);
    //computation of normalization
    std::cerr << "norm: " << boost::numeric::ublas::sum(transform_into_real_space(u)) << "\t";
    if (text_output) {
      spex_str<<"# alpha: "<<alpha[a]<<std::endl;
      for (std::size_t i=0; i<A.size(); ++i)
        spex_str << omega_coord(i) << " " << A[i] << "\n";
      spex_str << "\n";
    }
    //computation of probability
    lprob[a] = log_prob(u, alpha[a]);
    spectra[a] = A;
    //computation of chi2
    double chi_squared = chi2(transform_into_real_space(u));
    chi_sq[a] = chi_squared;
    if (verbose) std::cerr << "0.5*chi2  : " << 0.5*chi_squared;
    std::cerr << std::endl;
    if (text_output) print_chi2(transform_into_real_space(u), fits_str);
  }
  
  //everything from here on down is evaluation.
  if (text_output) {
    spex_str << "\n";
    for (std::size_t a=0; a<chi_sq.size(); ++a)
      chisq_str << alpha[a] << " " << chi_sq[a] << std::endl;
  }
  int a_chi = 0;
  double diff = std::abs(chi_sq[0]-ndat());
  for (std::size_t a=1; a<chi_sq.size(); ++a) {
    double diff_new = std::abs(chi_sq[a]-ndat());
    if (diff_new < diff) {
      diff = diff_new;
      a_chi = a;
    }
  }

  vector_type def = get_spectrum(transform_into_singular_space(Default()));
  if (text_output)
    for (std::size_t i=0; i<spectra[0].size(); ++i)
      chispec_str << omega_coord(i) << " " << spectra[a_chi][i]*norm << " " << def[i]*norm << std::endl;
  boost::numeric::ublas::vector<double>::const_iterator max_lprob = std::max_element(lprob.begin(), lprob.end());  
  const int max_a = max_lprob-lprob.begin();
  const double factor = chi_scale_factor(spectra[max_a], chi_sq[max_a], alpha[max_a]);
  if (verbose) std::cerr << "chi scale factor: " << factor << std::endl;
  
  alps::hdf5::archive ar(name+"out.h5", alps::hdf5::archive::WRITE);
  ar << alps::make_pvp("/alpha/values",alpha);
  
  vector_type om(spectra[0].size());
  for (int i=0;i<om.size();i++) om[i] = omega_coord(i);        
  ar<<alps::make_pvp("/spectrum/omega",om);
      
  //output 'maximum' spectral function (classical maxent metod)
  if (text_output) for (std::size_t i=0; i<spectra[0].size(); ++i)
    maxspec_str << omega_coord(i) << " " << spectra[max_a][i]*norm << " " << def[i]*norm << std::endl;
  {
    vector_type specmax = spectra[max_a]*norm,specchi = spectra[a_chi]*norm;
    ar << alps::make_pvp("/spectrum/chi",specchi);
    ar << alps::make_pvp("/spectrum/maximum",specmax);
  }
  vector_type prob(lprob.size());
  for (std::size_t a=0; a<prob.size(); ++a) 
    prob[a] = exp(lprob[a]-*max_lprob);
  double probnorm = 0;
  for (std::size_t a=0; a<prob.size()-1; ++a) 
    probnorm += 0.5*(prob[a]+prob[a+1])*(alpha[a]-alpha[a+1]);
  prob /= probnorm;
  ar << alps::make_pvp("/alpha/probability",prob);
  if (text_output) for (std::size_t a=0; a<prob.size(); ++a) {
    prob_str << alpha[a] << "\t" << prob[a] << "\n";
  }
  double postprobdef = 0;
  for (std::size_t a=0; a<lprob.size()-1; ++a) 
    postprobdef += 0.5*(exp(lprob[a])+exp(lprob[a+1]))*(alpha[a]-alpha[a+1]);
  std::cout << "posterior probability of the default model: " << postprobdef << std::endl;
  
  //compute 'average' spectral function (Brian's method)
  vector_type avspec(spectra[0].size());
  for (std::size_t i=0; i<avspec.size(); ++i) {
    avspec[i] = 0.;
    for (std::size_t a=0; a<prob.size()-1; ++a) 
      avspec[i] += 0.5*(prob[a]*spectra[a][i] +prob[a+1]*spectra[a+1][i])*(alpha[a]-alpha[a+1]);
  }
  //Estimate the variance for the spectrum
  vector_type varspec(spectra[0].size());
  for (std::size_t i=0; i<varspec.size(); ++i) {
    varspec[i] = 0.;
    for (std::size_t a=0; a<prob.size()-1; ++a)
    varspec[i] += 0.5*(prob[a]*(spectra[a][i]-avspec[i])*(spectra[a][i]-avspec[i]) + prob[a+1]*(spectra[a+1][i]-avspec[i])*(spectra[a+1][i]-avspec[i]))*(alpha[a]-alpha[a+1]);
  }
  avspec *= norm;
  varspec *= norm*norm;
  
  if (text_output) for (std::size_t  i=0; i<avspec.size(); ++i)
    avspec_str << omega_coord(i) << " " << avspec[i] << " " << def[i]*norm << std::endl;
  ar << alps::make_pvp("/spectrum/average",avspec);
  ar << alps::make_pvp("/spectrum/variance",varspec);
  
  if(Kernel_type=="anomalous"){ //for the anomalous function: use A(omega)=Im Sigma(omega)/(pi omega).
    std::ofstream maxspec_anom_str((name+"maxspec_anom.dat").c_str());
    std::ofstream avspec_anom_str (boost::filesystem::absolute(name+"avspec_anom.dat", dir).string().c_str());
    vector_type spec(avspec.size());
    for (std::size_t  i=0; i<avspec.size(); ++i){ 
      //if(omega_coord(i)>=0.)
      spec[i] = avspec[i]*omega_coord(i)*boost::math::constants::pi<double>();
      avspec_anom_str << omega_coord(i) << " " << avspec[i]*omega_coord(i)*boost::math::constants::pi<double>()<<std::endl;
    }
    ar << alps::make_pvp("/spectrum/anomalous/average",spec);
    for (std::size_t i=0; i<spectra[0].size(); ++i){
      //if(omega_coord(i)>=0.)
      spec[i] = spectra[max_a][i]*norm*omega_coord(i)*boost::math::constants::pi<double>();
      maxspec_anom_str << omega_coord(i) << " " << spectra[max_a][i]*norm*omega_coord(i)*boost::math::constants::pi<double>() << std::endl;
    }
    ar << alps::make_pvp("/spectrum/anomalous/maximum",spec);
  }
  if(Kernel_type=="bosonic"){ //for the anomalous function: use A(Omega)=Im chi(Omega)/(pi Omega) (as for anomalous)
    vector_type spec(avspec.size());
    for (std::size_t  i=0; i<avspec.size(); ++i){
      spec[i] = avspec[i]*omega_coord(i)*boost::math::constants::pi<double>();
    }
    if (text_output) {
      std::ofstream avspec_anom_str(boost::filesystem::absolute(name+"maxspec_bose.dat", dir).string().c_str());
      for (std::size_t  i=0; i<avspec.size(); ++i){
      //if(omega_coord(i)>=0.)
        avspec_anom_str << omega_coord(i) << " " << avspec[i]*omega_coord(i)*boost::math::constants::pi<double>()<<std::endl;
      }
    }
    ar << alps::make_pvp("/spectrum/bosonic/average",spec);
    for (std::size_t i=0; i<spectra[0].size(); ++i){
      //if(omega_coord(i)>=0.)
      spec[i] = spectra[max_a][i]*norm*omega_coord(i)*boost::math::constants::pi<double>();
    }
    if (text_output) {
      std::ofstream maxspec_anom_str (boost::filesystem::absolute(name+"avspec_bose.dat", dir).string().c_str());
      for (std::size_t i=0; i<spectra[0].size(); ++i){
        maxspec_anom_str << omega_coord(i) << " " << spectra[max_a][i]*norm*omega_coord(i)*boost::math::constants::pi<double>() << std::endl;
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
    std::ofstream maxspec_self_str(boost::filesystem::absolute(name+"maxspec_self.dat", dir).string().c_str());
    std::ofstream avspec_self_str (boost::filesystem::absolute(name+"avspec_self.dat", dir).string().c_str());
    for (std::size_t  i=0; i<avspec.size(); ++i){ 
      avspec_self_str << omega_coord(i) << " " << -avspec[i]*boost::math::constants::pi<double>()<<std::endl;
    }
    for (std::size_t i=0; i<spectra[0].size(); ++i){
      maxspec_self_str << omega_coord(i) << " " << -spectra[max_a][i]*norm*boost::math::constants::pi<double>() << std::endl;
    }
  }
 
  finished = true;
}



//this is the levenberg marquardt fitting procedure. It minimizes the quantity Q = 1/2 chi^2 - \alpha S
// 
MaxEntSimulation::vector_type MaxEntSimulation::levenberg_marquardt(vector_type u, const double alpha) const 
{
  using namespace boost::numeric;
  double mu = 1e-18;
  const double nu = 1.3;
  double Q1=0.;
  int it = 0;
  int it2 = 0;
  for (; it<max_it; it++) {
    vector_type delta;
    for (it2=0; it2<max_it; ++it2) {
      //compute change vector delta to u
      delta = iteration(u, alpha, mu);
      /*std::cout<<"delta is: "<<delta<<std::endl;
      vector_type z=transform_into_real_space(delta);
      for(int i=0;i<z.size();++i){
        std::cout<<omega_coord(i)<<" "<<z(i)<<std::endl;
      }*/
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
MaxEntSimulation::vector_type MaxEntSimulation::iteration(vector_type u, const double alpha, const double mu) const 
{
  using namespace boost::numeric;
  matrix_type M = left_side(u);
  for (std::size_t i=0; i<M.size1(); ++i) 
    M(i,i) += alpha + mu;
  vector_type b = right_side(u) + alpha*u;
  matrix_type B(b.size(),1);
  for (std::size_t i=0; i<M.size1(); ++i) 
    B(i,0) = -b[i];
  ublas::vector<fortran_int_t> ipiv(b.size());
  bindings::lapack::gesv(M, ipiv, B);
  return ublas::matrix_column<matrix_type>(B, 0);
}
