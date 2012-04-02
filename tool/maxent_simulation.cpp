/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2010 by Sebastian  Fuchs <fuchs@comp-phys.org>
*                       Thomas Pruschke <pruschke@comp-phys.org>
*                       Matthias Troyer <troyer@comp-phys.org>
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





MaxEntSimulation::MaxEntSimulation(const alps::ProcessList& w, const boost::filesystem::path& fn) 
 : alps::scheduler::Task(w,fn)
 , MaxEntHelper(parms)
 , alpha(parms["N_ALPHA"])                                                              //This is the # of \alpha parameters that should be tried.
 , norm(parms.value_or_default("NORM", 1.))                                             //The integral is normalized to NORM (use e.g. for self-energies
 , max_it(parms.value_or_default("MAX_IT", 1000))                                       //The number of iterations done in the root finding procedure
 , name(fn.filename().string(),0,fn.filename().string().size()-6)
 , dir(fn.branch_path())
 , spex_str(boost::filesystem::absolute(name+"spex.dat", dir).string().c_str())
 , chisq_str(boost::filesystem::absolute(name+"chi2.dat", dir).string().c_str())
 , avspec_str(boost::filesystem::absolute(name+"avspec.dat", dir).string().c_str())
 , maxspec_str(boost::filesystem::absolute(name+"maxspec.dat", dir).string().c_str())
 , chispec_str(boost::filesystem::absolute(name+"chispec.dat", dir).string().c_str())
 , prob_str(boost::filesystem::absolute(name+"prob.dat", dir).string().c_str())
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
  if (finished())
    return;
  vector_type lprob(alpha.size());
  vector_type chi_sq(alpha.size());
  std::vector<vector_type> spectra(alpha.size());
  vector_type u = transform_into_singular_space(Default());
  for (std::size_t a=0; a<alpha.size(); ++a) {
    std::cerr << "alpha it: " << a << "\t";
    u = levenberg_marquardt(u, alpha[a]);
    vector_type A = get_spectrum(u);
    std::cerr << "norm: " << boost::numeric::ublas::sum(transform_into_real_space(u)) << "\t";
    for (std::size_t i=0; i<A.size(); ++i) 
      spex_str << alpha[a] << " " << omega_coord(i) << " " << A[i] << "\n";
    spex_str << "\n";
    lprob[a] = log_prob(u, alpha[a]);
    spectra[a] = A;
    double chi_squared = chi2(transform_into_real_space(u));
    chi_sq[a] = chi_squared;
    std::cerr << "chi2  : " << chi_squared << std::endl;
  }
  spex_str << "\n";
  for (std::size_t a=0; a<chi_sq.size(); ++a) 
    chisq_str << alpha[a] << " " << chi_sq[a] << std::endl;
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
  for (std::size_t i=0; i<spectra[0].size(); ++i) 
    chispec_str << omega_coord(i) << " " << spectra[a_chi][i]*norm << " " << def[i]*norm << std::endl;
  boost::numeric::ublas::vector<double>::const_iterator max_lprob = std::max_element(lprob.begin(), lprob.end());  
  const int max_a = max_lprob-lprob.begin();
  const double factor = chi_scale_factor(spectra[max_a], chi_sq[max_a], alpha[max_a]);
  std::cerr << "chi scale factor: " << factor << std::endl;

  //output 'maximum' spectral function (classical maxent metod)
  for (std::size_t i=0; i<spectra[0].size(); ++i) 
    maxspec_str << omega_coord(i) << " " << spectra[max_a][i]*norm << " " << def[i]*norm << std::endl;
  vector_type prob(lprob.size());
  for (std::size_t a=0; a<prob.size(); ++a) 
    prob[a] = exp(lprob[a]-*max_lprob);
  double probnorm = 0;
  for (std::size_t a=0; a<prob.size()-1; ++a) 
    probnorm += 0.5*(prob[a]+prob[a+1])*(alpha[a]-alpha[a+1]);
  prob /= probnorm;
  for (std::size_t a=0; a<prob.size(); ++a) {
    prob_str << alpha[a] << "\t" << prob[a] << "\n";
  }
  double postprobdef = 0;
  for (std::size_t a=0; a<lprob.size()-1; ++a) 
    postprobdef += 0.5*(exp(lprob[a])+exp(lprob[a+1]))*(alpha[a]-alpha[a+1]);
  std::cout << "posterior probability of the default model: " << postprobdef << std::endl;
 
  //compute 'average' spectral function (Brian's metod)
  vector_type avspec(spectra[0].size());
  for (std::size_t i=0; i<avspec.size(); ++i) {
    avspec[i] = 0.;
    for (std::size_t a=0; a<prob.size()-1; ++a) 
      avspec[i] += 0.5*(prob[a]*spectra[a][i] +prob[a+1]*spectra[a+1][i])*(alpha[a]-alpha[a+1]);
  }
  for (std::size_t  i=0; i<avspec.size(); ++i) 
    avspec_str << omega_coord(i) << " " << avspec[i]*norm << " " << def[i]*norm << std::endl;
  if(parms["KERNEL"]=="anomalous"){ //for the anomalous function: use A(omega)=Im Sigma(omega)/pi omega. 
    std::ofstream maxspec_anom_str(boost::filesystem::absolute(name+"maxspec_anom.dat", dir).string().c_str());
    std::ofstream avspec_anom_str (boost::filesystem::absolute(name+"avspec_anom.dat", dir).string().c_str());
    for (std::size_t  i=0; i<avspec.size(); ++i){ 
      //if(omega_coord(i)>=0.)
        avspec_anom_str << omega_coord(i) << " " << avspec[i]*norm*omega_coord(i)*M_PI<<std::endl;
    }
    for (std::size_t i=0; i<spectra[0].size(); ++i){
      //if(omega_coord(i)>=0.)
        maxspec_anom_str << omega_coord(i) << " " << spectra[max_a][i]*norm*omega_coord(i)*M_PI << std::endl;
    }
  }
  if(parms.defined("SELF")){ //for the self energy: use Im Sigma(omega)=-A(omega)*pi
    std::ofstream maxspec_self_str(boost::filesystem::absolute(name+"maxspec_self.dat", dir).string().c_str());
    std::ofstream avspec_self_str (boost::filesystem::absolute(name+"avspec_self.dat", dir).string().c_str());
    for (std::size_t  i=0; i<avspec.size(); ++i){ 
        avspec_self_str << omega_coord(i) << " " << -avspec[i]*norm*M_PI<<std::endl;
    }
    for (std::size_t i=0; i<spectra[0].size(); ++i){
        maxspec_self_str << omega_coord(i) << " " << -spectra[max_a][i]*norm*M_PI << std::endl;
    }
  }
  finish();
}




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
      delta = iteration(u, alpha, mu);
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
  std::cerr <<"Iterations: " << it+1 << "\t"
            << "Q: " << Q1 << "\t";
  return u;
}




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



//this function is crap. Why do we need it? It has zero content!
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
