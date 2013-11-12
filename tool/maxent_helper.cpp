/*****************************************************************************
 *
 * ALPS Project Applications
 *
 * Copyright (C) 2010 by Sebastian  Fuchs <fuchs@comp-phys.org>
 *                       Thomas Pruschke <pruschke@comp-phys.org>
 *                       Matthias Troyer <troyer@comp-phys.org>
 *               2011 by Emanuel Gull  <gull@phys.columbia.edu>
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
#include <boost/numeric/bindings/lapack/driver/posv.hpp>
#include <boost/numeric/bindings/lapack/computational/potrf.hpp>
#include <boost/numeric/bindings/lapack/driver/syev.hpp>
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/upper.hpp>
#include <boost/numeric/bindings/lower.hpp>

/*
MaxEntHelper::MaxEntHelper(const alps::Parameters& p) : 
MaxEntParameters(p) , def_(nfreq())
{
  for (int i=0; i<nfreq(); ++i) 
    def_[i] = MaxEntParameters::Default().D(omega_coord(i)) * delta_omega(i); 
  //std::cout<<"sum of def: "<<sum(def_)<<std::endl;
  def_ /= sum(def_);
  std::ofstream out;
  out.open("deltaOmega.dat");
  for (int i=0; i<nfreq(); ++i)
    out << i<<" "<<delta_omega(i) << " "<<def_[i]<<std::endl;
}
 */

MaxEntHelper::MaxEntHelper(const alps::params& p) :
MaxEntParameters(p) , def_(nfreq())
{
    for (int i=0; i<nfreq(); ++i)
        def_[i] = MaxEntParameters::Default().D(omega_coord(i)) * delta_omega(i);
    //std::cout<<"sum of def: "<<sum(def_)<<std::endl;
    def_ /= sum(def_);
    std::ofstream out;
    out.open("deltaOmega.dat");
    for (int i=0; i<nfreq(); ++i)
        out << i<<" "<<delta_omega(i) << " "<<def_[i]<<std::endl;
}


//the opposite of 'transform_into_real_space'; takes a vector 'A' and makes a vector 'u' out of it:
// u = V^T* log(A/Default)
MaxEntHelper::vector_type MaxEntHelper::transform_into_singular_space(vector_type A) const
{
  for (unsigned int i=0; i<A.size(); ++i) {
    A[i] /= Default(i);
    A[i] = A[i]==0. ? 0. : log(A[i]);
  }
  return prec_prod(Vt(), A);
}



//returns exp(V^T*u)*Default(i). This quantity is then usually called 'A'
MaxEntHelper::vector_type MaxEntHelper::transform_into_real_space(vector_type u) const
{
  u = prec_prod(trans(Vt()), u);
  for (unsigned int i=0; i<u.size(); ++i) {
    u[i] = exp(u[i]);
    u[i] *= Default(i);
  }
  return u; 
}




MaxEntHelper::vector_type MaxEntHelper::get_spectrum(const vector_type& u) const
{
  vector_type A = transform_into_real_space(u);
  for (unsigned int i=0; i<A.size(); ++i) 
    A[i] /= delta_omega(i);
  return A;
}



//'left side' is defined as Sigma*(V^T*RealSpace(u)*V)*Sigma
//see Bryan's paper near Eq. 11 
MaxEntHelper::matrix_type MaxEntHelper::left_side(const vector_type& u) const 
{
  vector_type A = transform_into_real_space(u);
  matrix_type M = trans(Vt());
  for (unsigned int i=0; i<M.size1(); ++i) 
    for (unsigned int j=0; j<M.size2(); ++j) 
      M(i,j) *= A[i];
  M = prec_prod(Vt(), M);
  M = prec_prod(Sigma() ,M);
  M = prec_prod(Sigma(), M);
  M *= 2./ndat();
  return M;
}


//this function computes
//Sigma*U^T*(K*RealSpace(u)-y)
//up to a factor of 2./ndat(). Compare this to Eq. D.12 in Sebastian's thesis
MaxEntHelper::vector_type MaxEntHelper::right_side(const vector_type& u) const 
{
  vector_type b = 2./ndat()*(prec_prod(K(), transform_into_real_space(u)) - y());
  b = prec_prod(trans(U()), b);
  b = prec_prod(Sigma(), b);
  return b;
}



//this function constructs delta \dot (V^T*RealSpace(u)*V)
double MaxEntHelper::step_length(const vector_type& delta, const vector_type& u) const 
{
  vector_type A = transform_into_real_space(u);
  matrix_type L = trans(Vt());
  for (unsigned int i=0; i<L.size1(); ++i) 
    for (unsigned int j=0; j<L.size2(); ++j) 
      L(i,j) *= A[i];
  L = prec_prod(Vt(), L);
  return inner_prod(delta, prec_prod(L, delta));
}




double MaxEntHelper::convergence(const vector_type& u, const double alpha) const 
{
  using namespace boost::numeric::ublas;
  vector_type A = transform_into_real_space(u);
  matrix_type L = trans(Vt());
  for (unsigned int i=0; i<L.size1(); ++i) 
    for (unsigned int j=0; j<L.size2(); ++j) 
      L(i,j) *= A[i];
  L = prec_prod(Vt(), L);
  vector_type alpha_dSdu = -alpha*prec_prod(L, u);
  vector_type dLdu = prec_prod(L, right_side(u));
  vector_type diff = alpha_dSdu - dLdu;
  double denom = norm_2(alpha_dSdu) + norm_2(dLdu);
  denom = denom*denom;
  return 2*inner_prod(diff, diff)/denom;
}



double MaxEntHelper::log_prob(const vector_type& u, const double alpha) const
{
  matrix_type L = prec_prod_trans(K(), K());
  const vector_type A = transform_into_real_space(u);
  for (unsigned int i=0; i<L.size1(); ++i)
    for (unsigned int j=0; j<L.size2(); ++j)
      L(i,j) *= sqrt(A[i])*sqrt(A[j]);
  for (unsigned int i=0; i<L.size1(); ++i)
    L(i,i) += alpha;
  boost::numeric::bindings::lapack::potrf(boost::numeric::bindings::lower(L));
  double log_det = 0.;
  for (unsigned int i=0; i<L.size1(); ++i) 
    log_det  += log(L(i,i)*L(i,i));
  return 0.5*( (nfreq())*log(alpha) - log_det ) - Q(u, alpha);
}



double MaxEntHelper::chi_scale_factor(vector_type A, const double chi_sq, const double alpha) const
{
  for (unsigned int i=0; i<A.size(); ++i) 
    A[i] *= delta_omega(i);
  using namespace boost::numeric;
  matrix_type L = prec_prod_trans(K(), K());
  for (unsigned int i=0; i<L.size1(); ++i)
    for (unsigned int j=0; j<L.size2(); ++j)
      L(i,j) *= sqrt(A[i])*sqrt(A[j]);
  vector_type lambda(L.size1());
  bindings::lapack::syev('N', bindings::upper(L) , lambda, bindings::lapack::optimal_workspace());
  double Ng = 0.;
  for (unsigned int i=0; i<lambda.size(); ++i) {
    if (lambda[i]>=0) 
      Ng += lambda[i]/(lambda[i]+alpha);
  }
  std::cerr << "Ng: " << Ng << std::endl;
  std::cerr << "chi2 max: " << chi_sq << std::endl;
  return sqrt(chi_sq/(ndat()-Ng));
}


//This function computes chi^2 as in equation D.6 in Sebastian's thesis
double MaxEntHelper::chi2(const vector_type& A) const 
{
  vector_type del_G = prec_prod(K(), A) - y();
  
  /*std::cout<<"in computation of chi2:"<<std::endl;
   for(int i=0;i<y().size();++i){
   std::cout<<i<<" "<<prec_prod(K(), A)[i]<<" "<<y()[i]<<std::endl;
   }*/
  
  double c = 0;
  for (unsigned int i=0; i<del_G.size(); ++i) 
    c += del_G[i]*del_G[i];
  return c;
}
void MaxEntHelper::print_chi2(const vector_type& A, std::ostream &os) const 
{
  vector_type backcont=prec_prod(K(), A);
  vector_type defaultm=prec_prod(K(), Default());
  os<<"#first column: index (Matsubara frequency). second column: fitted function. third: input data. fourth: default model."<<std::endl;
  for(int i=0;i<y().size();++i){
    os<<i<<" "<<backcont[i]<<" "<<y()[i]<<" "<<defaultm[i]<<std::endl;
  }
  os<<std::endl;
}


//This function computes the entropy as in Eq. D.7 in Sebastian's thesis
double MaxEntHelper::entropy(const vector_type& A) const 
{
  double S = 0;
  for (unsigned int i=0; i<A.size(); ++i) {
    double lg = A[i]==0. ? 0. : log(A[i]/Default(i));
    S += A[i] - Default(i) - A[i]*lg;
  }
  return S;
}

// Generic principal-value integration. Assumes that spectral function does not vary too
// strongly within [x_{i+1},x_i]
MaxEntHelper::vector_type MaxEntHelper::PrincipalValue(const vector_type &w,const vector_type &a) const
{
    int N = w.size();
    vector_type r(N);
    for (int i=2;i<N-2;i++) {
      double scr = -a[i]*std::log(std::abs((w[i+1]-w[i])/(w[i]-w[i-1])));
      for (int j=0;j<i-1;j++) scr -= 0.5*(a[j]+a[j+1])*std::log(std::abs((w[i]-w[j+1])/(w[i]-w[j])));
      for (int j=i+1;j<N-2;j++) scr -= 0.5*(a[j]+a[j+1])*std::log(std::abs((w[j+1]-w[i])/(w[j]-w[i])));
      r[i] = scr;
    }
    r[0] = w[2]*r[2]/w[0];
    r[1] = w[2]*r[2]/w[1];
    r[N-2] = w[N-3]*r[N-3]/w[N-2];
    r[N-1] = w[N-3]*r[N-3]/w[N-1];
    return r;
}


