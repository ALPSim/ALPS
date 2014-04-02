/*****************************************************************************
 *
 * ALPS Project Applications
 *
 * Copyright (C) 2010 by Sebastian  Fuchs <fuchs@comp-phys.org>
 *                       Thomas Pruschke <pruschke@comp-phys.org>
 *                       Matthias Troyer <troyer@comp-phys.org>
 *               2012 by Emanuel Gull <gull@pks.mpg.de>
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

#include "maxent.hpp"
#include <alps/config.h> // needed to set up correct bindings
#include <boost/filesystem/operations.hpp>
#include <boost/numeric/bindings/lapack/driver/gesv.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas.hpp>





MaxEntSimulation::MaxEntSimulation(const alps::params &parms,const std::string &outfile)
: MaxEntHelper(parms)
, alps::mcbase(parms)
, alpha((int)parms["N_ALPHA"])              //This is the # of \alpha parameters that should be tried.
, norm(parms["NORM"]|1.0)                                             //The integral is normalized to NORM (use e.g. for self-energies
, max_it(parms["MAX_IT"]|1000)                                       //The number of iterations done in the root finding procedure
, name(outfile,0,outfile.size()-6)
, Kernel_type(parms["KERNEL"]|"")
, finished(false)
, verbose(parms["VERBOSE"]|false)
, text_output(parms["TEXT_OUTPUT"]|false)
, self(parms["SELF"]|false)
{
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
      spec[i] = avspec[i]*omega_coord(i)*M_PI;
      avspec_anom_str << omega_coord(i) << " " << avspec[i]*omega_coord(i)*M_PI<<std::endl;
    }
    ar << alps::make_pvp("/spectrum/anomalous/average",spec);
    for (std::size_t i=0; i<spectra[0].size(); ++i){
      //if(omega_coord(i)>=0.)
      spec[i] = spectra[max_a][i]*norm*omega_coord(i)*M_PI;
      maxspec_anom_str << omega_coord(i) << " " << spectra[max_a][i]*norm*omega_coord(i)*M_PI << std::endl;
    }
    ar << alps::make_pvp("/spectrum/anomalous/maximum",spec);
  }
  if(Kernel_type=="bosonic"){ //for the anomalous function: use A(Omega)=Im chi(Omega)/(pi Omega) (as for anomalous)
    vector_type spec(avspec.size());
    for (std::size_t  i=0; i<avspec.size(); ++i){
      spec[i] = avspec[i]*omega_coord(i)*M_PI;
    }
    if (text_output) {
      std::ofstream avspec_anom_str(boost::filesystem::absolute(name+"maxspec_bose.dat", dir).string().c_str());
      for (std::size_t  i=0; i<avspec.size(); ++i){
      //if(omega_coord(i)>=0.)
        avspec_anom_str << omega_coord(i) << " " << avspec[i]*omega_coord(i)*M_PI<<std::endl;
      }
    }
    ar << alps::make_pvp("/spectrum/bosonic/average",spec);
    for (std::size_t i=0; i<spectra[0].size(); ++i){
      //if(omega_coord(i)>=0.)
      spec[i] = spectra[max_a][i]*norm*omega_coord(i)*M_PI;
    }
    if (text_output) {
      std::ofstream maxspec_anom_str (boost::filesystem::absolute(name+"avspec_bose.dat", dir).string().c_str());
      for (std::size_t i=0; i<spectra[0].size(); ++i){
        maxspec_anom_str << omega_coord(i) << " " << spectra[max_a][i]*norm*omega_coord(i)*M_PI << std::endl;
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
      avspec_self_str << omega_coord(i) << " " << -avspec[i]*M_PI<<std::endl;
    }
    for (std::size_t i=0; i<spectra[0].size(); ++i){
      maxspec_self_str << omega_coord(i) << " " << -spectra[max_a][i]*norm*M_PI << std::endl;
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



//this function is nonsensical. Why do we need it? It has zero content!
void MaxEntSimulation::write_xml_body(alps::oxstream& out, const boost::filesystem::path&, bool write_all_xml) const
{
  if (write_all_xml) {
    out << alps::start_tag("AVERAGES");
    out << alps::start_tag("SCALAR_AVERAGE") << alps::attribute("name","Zeug") << alps::no_linebreak
    << alps::start_tag("MEAN") << 42 << alps::end_tag("MEAN")
    << alps::end_tag("SCALAR_AVERAGE");
    out << alps::end_tag("AVERAGES");
  }
}
