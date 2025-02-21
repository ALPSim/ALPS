/****************************************************************************
*
* ALPS Project Applications: Directed Worm Algorithm  
*
* Copyright (C) 2013 by Matthias Troyer  <troyer@phys.ethz.ch> ,
*                       Lode Pollet      <pollet@phys.ethz.ch> ,
*                       Ping Nang Ma     <pingnang@phys.ethz.ch> 
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef BANDSTRUCTURE_HPP
#define BANDSTRUCTURE_HPP

#include <vector>
#include <boost/multi_array.hpp>
#include <alps/numeric/vector_functions.hpp>
#include <alps/ngs/boost_python.hpp>
#include <boost/bind.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <alps/python/numpy_array.hpp>


// LAPACK Library: dsteqr -- description: diagonalize a real tridiagonal matrix
//                        -- arguments:   compz, n, d, e, z, ldz, work, info
extern "C" void dsteqr_(char*,int*,double*,double*,double*,int*,double*,int*); // (compz,n,d,e,z,ldz,work,info) --diagonalize a tridiagonal matrix


class bandstructure
{
public:
  bandstructure(double V0_, double lambda_, double a_, double m_, unsigned int L_, int Mmax_=10);
  bandstructure(boost::python::object V0_, boost::python::object lambda_, double a_, double m_, unsigned int L_, int Mmax_=10); 

  std::vector<double> get_t()    {  if(!is_evaluated)  evaluate();  return t;     }
  double              get_U()    {  if(!is_evaluated)  evaluate();  return U;     }
  std::vector<double> get_Ut()   {  if(!is_evaluated)  evaluate();  using alps::numeric::operator/; return U/t;  }
  std::vector<double> get_norm() {  if(!is_evaluated)  evaluate();  return norm;  }
  std::vector<double> get_q  (unsigned int i)  {  if(!is_evaluated)  evaluate();  return q[i];    }
  std::vector<double> get_wk2(unsigned int i)  {  if(!is_evaluated)  evaluate();  return wk2[i];  }
  double              get_wk2_c()  {  if(!is_evaluated)  evaluate();  return wk2_c;  }
  double              get_wk2_d()  {  if(!is_evaluated)  evaluate();  return wk2_d;  }

  void output(std::ostream & out);
  std::string representation()   { std::ostringstream oss; output(oss); return oss.str(); }

private:
  void evaluate();

  static const double pi; 
  static const double hbar;
  static const double k;  
  static const double amu;
  static const double a0; 

  bool is_evaluated;

  std::vector<double>  V0;        // in recoil energies (note that they are different in every direction)
  std::vector<double>  lambda;    // in nanometer
  std::vector<double>  Er2nK;
  unsigned int         L;
  double               g;         // in nK  

  std::vector<double>  t;
  double               U;

  std::vector<double>  norm;

  std::vector<std::vector<double> > q;
  std::vector<std::vector<double> > wk2;  

  double wk2_c;
  double wk2_d;   

  int Mmax;
};

const double bandstructure::pi   = 3.141592654;
const double bandstructure::hbar = 1.05457148;
const double bandstructure::k    = 1.3806503;
const double bandstructure::amu  = 1.66053886;
const double bandstructure::a0   = 0.052917720859;


// Constructor
bandstructure::bandstructure(double V0_, double lambda_, double a_, double m_, unsigned int L_, int Mmax_)
  : is_evaluated (false)
  , V0     (std::vector<double>(3, V0_))
  , lambda (std::vector<double>(3, lambda_))
  , L      (L_)
  , Mmax   (Mmax_)
{
  if (a_ == 0. || m_ == 0. || L == 0)   
    boost::throw_exception(std::runtime_error("Illegal initialization parameters for bandstructure class"));

  for (int i=0; i<3; ++i) 
    if (V0[i] == 0. || lambda[i] == 0.)
      boost::throw_exception(std::runtime_error("Illegal initialization parameters for bandstructure class"));

  using boost::numeric::operators::operator*;
  using alps::numeric::operator/;
  Er2nK = (2e9 * pi*pi * hbar*hbar) / (m_ * amu * k * lambda * lambda);

  g = ((32e9 * pi * a0 * hbar*hbar * a_) / (amu * k * m_ * lambda[0]*lambda[1]*lambda[2]));

  t.resize(3,0.);
  U = g;

  norm.resize(3,0.);

  q  .resize(3);
  wk2.resize(3);

  wk2_c = 1.;
  wk2_d = 1.;
}

bandstructure::bandstructure(boost::python::object V0_, boost::python::object lambda_, double a_, double m_, unsigned int L_, int Mmax_)
  : is_evaluated (false)
  , L            (L_)
  , Mmax         (Mmax_)
{
  alps::python::numpy::convert(V0_, V0);
  alps::python::numpy::convert(lambda_, lambda);

  if (V0.empty() || lambda.empty() || a_ == 0. || m_ == 0. || L == 0) 
    boost::throw_exception(std::runtime_error("Illegal initialization parameters for bandstructure class"));

  V0.resize    (3, V0.back());
  lambda.resize(3, lambda.back());

  for (int i=0; i<3; ++i) 
    if (V0[i] == 0. || lambda[i] == 0.)
      boost::throw_exception(std::runtime_error("Illegal initialization parameters for bandstructure class"));

  using boost::numeric::operators::operator*;
  using alps::numeric::operator/;
  Er2nK = (2e9 * pi*pi * hbar*hbar) / (m_ * amu * k * lambda * lambda); 

  g = ((32e9 * pi * a0 * hbar*hbar * a_) / (amu * k * m_ * lambda[0]*lambda[1]*lambda[2]));

  t.resize(3,0.);
  U = g;

  norm.resize(3,0.);

  q  .resize(3);
  wk2.resize(3);

  wk2_c = 1.;
  wk2_d = 1.;
}


void bandstructure::output(std::ostream & out)
{
  using alps::numeric::operator<<;
  out << "\nOptical lattice:\n"
      << "================\n"
      << "V0    [Er] = " << V0 << "\n"
      << "lamda [nm] = " << lambda << "\n"
      << "Er2nK      = " << Er2nK  << "\n"
      << "L          = " << L << "\n"
      << "g          = " << g << "\n"
      ;

  out << "\nBand structure:\n"
      << "===============\n"
      << "t [nK] : " << get_t()  << "\n"
      << "U [nK] : " << get_U()  << "\n"
      << "U/t    : " << get_Ut() << "\n\n"
      << "wk2[0 ,0 ,0 ] : " << get_wk2_c()  << "\n"
      << "wk2[pi,pi,pi] : " << get_wk2_d()  << "\n"
      ;
}


void bandstructure::evaluate()
{
  for (int alpha=0; alpha<3; ++alpha)   // in each dimension: alpha=x,y,z
  {
    int size = 2*Mmax+1;

    std::vector<double> k;
    std::vector<double> energy;
    boost::multi_array<double, 2> c(boost::extents[L][size]); 

    for (int idx=0; idx<L; ++idx)
      k.push_back(static_cast<double>(idx)/L);

    // diagonalize for all values of k
    for (int idx=0; idx<L; ++idx)
    {
      std::vector<double> diagonal;
      for (int m=-Mmax; m<=Mmax; ++m)
        diagonal.push_back(4*alps::numeric::sq(m + k[idx]) + V0[alpha]/2.);
      std::vector<double> offdiagonal(2*Mmax, -V0[alpha]/4.);

      char                 type = 'I';
      int                  info;
      std::vector<double>  tmp(2*size-2);

      boost::multi_array<double, 2> c_tmp(boost::extents[size][size]);

      dsteqr_(&type, &size, &diagonal[0], &offdiagonal[0], &c_tmp[0][0], &size, &tmp[0], &info);

      energy.push_back(diagonal[0]);

      if (c_tmp[0][Mmax] >= 0.)
        std::copy(&c_tmp[0][0], &c_tmp[0][size], &c[idx][0]);
      else
        std::transform(&c_tmp[0][0], &c_tmp[0][size], &c[idx][0], std::negate<double>());
    }

    // evaluate t
    double this_t=0.;
    for (int idx=0; idx<L; ++idx)
      this_t += energy[idx] * std::cos(2*pi*k[idx]);
    this_t /= (L);
    t[alpha] = -(this_t*Er2nK[alpha]);

    // collapse c as much as possible 
    std::vector<double> c_tmp(c.shape()[1], 0.);
    for (int midx=0; midx<c_tmp.size(); ++midx)
      for (int idx=0; idx<L; ++idx)
        c_tmp[midx] += c[idx][midx];
    int M = Mmax;
    while (c_tmp[Mmax-M] < 1e-3 && c_tmp[Mmax+M] < 1e-3)
      --M;

    typedef boost::multi_array_types::index_range range;
    typedef range::index index;
    c = c[boost::indices[range()][index(Mmax-M) <= range() <= index(Mmax+M)]];
    Mmax = M;

    // convenient values
    std::vector<double> m;
    for (int midx=0; midx<=2*Mmax; ++midx)
      m.push_back(midx-Mmax);

    boost::multi_array<double,2> k_eff(boost::extents[L][m.size()]);
    for (int idx=0; idx<L; ++idx)
    for (int midx=0; midx<m.size(); ++midx)
      k_eff[idx][midx] = k[idx] + m[midx];

    // calculate wannier, norm and U
    std::vector<double> this_wannier(1000, 0.);  // x=0, 0.01, 0.02, ..., 10.00
    for (int i=0; i<1000; ++i) {
      double x = static_cast<double>(i)/100.;
      for (int idx=0; idx<L; ++idx)
      for (int midx=0; midx<m.size(); ++midx)
        this_wannier[i] += c[idx][midx] * std::cos(2. * pi * k_eff[idx][midx] * x);
      this_wannier[i] /= L;
    }
    std::vector<double> this_wannier2 = alps::numeric::sq(this_wannier);
    std::vector<double> this_wannier4 = alps::numeric::sq(this_wannier2);
   
    norm[alpha] = 0.01 * (2. * std::accumulate(this_wannier2.begin(), this_wannier2.end(), 0.) - this_wannier2[0] - this_wannier2.back());
    U *= 0.01 * (2. * std::accumulate(this_wannier4.begin(), this_wannier4.end(), 0.) - this_wannier4[0] - this_wannier4.back());

    // calculate wk2
    std::map<double, double> wk;
    for (int idx=0; idx<L; ++idx)
    for (int midx=0; midx<m.size(); ++midx)
      wk.insert(std::make_pair(k_eff[idx][midx], c[idx][midx]));

    for (std::map<double, double>::iterator it=wk.begin(); it!=wk.end(); ++it) {
      q  [alpha].push_back(it->first);
      wk2[alpha].push_back(it->second);
    } 
    using boost::numeric::operators::operator/;
    wk2[alpha] = alps::numeric::sq(wk2[alpha] / std::sqrt(static_cast<double>(L)));

    wk2_c *= alps::numeric::sq(wk[0.]) / (static_cast<double>(L));
    wk2_d *= alps::numeric::sq(wk[0.5]) / (static_cast<double>(L));
  }

  is_evaluated = true;
}


#endif
