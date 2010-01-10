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
#include <boost/numeric/bindings/lapack/posv.hpp>
#include <boost/numeric/bindings/lapack/syev.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>


MaxEntHelper::MaxEntHelper(const alps::Parameters& p) : 
  MaxEntParameters(p) , def_(nfreq())
{
  for (int i=0; i<nfreq(); ++i) 
    def_[i] = MaxEntParameters::Default().D(omega_coord(i)) * delta_omega(i); 
  def_ /= sum(def_);
  std::ofstream out;
  out.open("deltaOmega.dat");
  for (int i=0; i<nfreq(); ++i)
    out << delta_omega(i) << std::endl;
}



MaxEntHelper::vector_type MaxEntHelper::transform_into_singular_space(vector_type A) const
{
  for (int i=0; i<A.size(); ++i) {
    A[i] /= Default(i);
    A[i] = A[i]==0. ? 0. : log(A[i]);
  }
  return boost::numeric::ublas::prec_prod(Vt(), A);
}




MaxEntHelper::vector_type MaxEntHelper::transform_into_real_space(vector_type u) const
{
  using namespace boost::numeric::ublas;
  u = prec_prod(trans(Vt()), u);
  for (int i=0; i<u.size(); ++i) {
    u[i] = exp(u[i]);
    u[i] *= Default(i);
  }
  return u; 
}




MaxEntHelper::vector_type MaxEntHelper::get_spectrum(const vector_type& u) const
{
  vector_type A = transform_into_real_space(u);
  for (int i=0; i<A.size(); ++i) 
    A[i] /= delta_omega(i);
  return A;
}




MaxEntHelper::matrix_type MaxEntHelper::left_side(const vector_type& u) const 
{
  using namespace boost::numeric::ublas;
  vector_type A = transform_into_real_space(u);
  matrix_type M = trans(Vt());
  for (int i=0; i<M.size1(); ++i) 
    for (int j=0; j<M.size2(); ++j) 
      M(i,j) *= A[i];
  M = prec_prod(Vt(), M);
  M = prec_prod(Sigma() ,M);
  M = prec_prod(Sigma(), M);
  M *= 2./ndat();
  return M;
}




MaxEntHelper::vector_type MaxEntHelper::right_side(const vector_type& u) const 
{
  using namespace boost::numeric::ublas;
  vector_type b = 2./ndat()*(prec_prod(K(), transform_into_real_space(u)) - y());
  b = prec_prod(trans(U()), b);
  b = prec_prod(Sigma(), b);
  return b;
}




double MaxEntHelper::step_length(const vector_type& delta, const vector_type& u) const 
{
  using namespace boost::numeric::ublas;
  vector_type A = transform_into_real_space(u);
  matrix_type L = trans(Vt());
  for (int i=0; i<L.size1(); ++i) 
    for (int j=0; j<L.size2(); ++j) 
      L(i,j) *= A[i];
  L = prec_prod(Vt(), L);
  return inner_prod(delta, prec_prod(L, delta));
}




double MaxEntHelper::convergence(const vector_type& u, const double alpha) const 
{
  using namespace boost::numeric::ublas;
  vector_type A = transform_into_real_space(u);
  matrix_type L = trans(Vt());
  for (int i=0; i<L.size1(); ++i) 
    for (int j=0; j<L.size2(); ++j) 
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
  using namespace boost::numeric;
  matrix_type L = ublas::prec_prod(ublas::trans(K()), K());
  const vector_type A = transform_into_real_space(u);
  for (int i=0; i<L.size1(); ++i)
    for (int j=0; j<L.size2(); ++j)
      L(i,j) *= 2./ndat()*sqrt(A[i])*sqrt(A[j]);
  for (int i=0; i<L.size1(); ++i)
    L(i,i) += alpha;
  bindings::lapack::potrf('L', L);
  double log_det = 0.;
  for (int i=0; i<L.size1(); ++i) 
    log_det  += log(L(i,i)*L(i,i));
  double lprob = 0.5*( (nfreq())*log(alpha) - log_det ) - Q(u, alpha);
  /*double gamma = 100;
  for (int g=0; g<50; ++g) {
    gamma *= std::pow(0.01/100, 1./double(50-1));
    double lp = 0.5*nfreq()*log(gamma) + 0.5*( (nfreq())*log(alpha) - log_det ) - Q(u, alpha, gamma);
    std::cout << alpha << "\t" << gamma << "\t" << lp << "\n";
    }*/
  return lprob;
}



double MaxEntHelper::chi_scale_factor(vector_type A, const double chi_sq, const double alpha) const
{
  for (int i=0; i<A.size(); ++i) 
    A[i] *= delta_omega(i);
  using namespace boost::numeric;
  matrix_type L = ublas::prec_prod(ublas::trans(K()), K());
  for (int i=0; i<L.size1(); ++i)
    for (int j=0; j<L.size2(); ++j)
    L(i,j) *= 2./ndat()*sqrt(A[i])*sqrt(A[j]);
  vector_type lambda(L.size1());
  bindings::lapack::syev('N', 'U', L , lambda, bindings::lapack::optimal_workspace());
  double Ng = 0.;
  for (int i=0; i<lambda.size(); ++i) {
    if (lambda[i]>=0) 
      Ng += lambda[i]/(lambda[i]+alpha);
  }
  Ng /= ndat();
  std::cerr << "Ng: " << Ng << std::endl;
  std::cerr << "chi2 max: " << chi_sq << std::endl;
  return sqrt(chi_sq/(1-Ng));
}



double MaxEntHelper::chi2(const vector_type& A) const 
{
  vector_type del_G = prec_prod(K(), A) - y();
  double c = 0;
  for (int i=0; i<del_G.size(); ++i) 
    c += del_G[i]*del_G[i];
  c /= ndat();
  return c;
}



double MaxEntHelper::entropy(const vector_type& A) const 
{
  double S = 0;
  for (int i=0; i<A.size(); ++i) {
    double lg = A[i]==0. ? 0. : log(A[i]/Default(i));
    S += A[i] - Default(i) - A[i]*lg;
  }
  return S;
}



