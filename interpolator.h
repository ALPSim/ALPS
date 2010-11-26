/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2009 by Sebastian Fuchs <fuchs@comp-phys.org>
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

#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H


#include "cluster.h"
#include "alps/parameter.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>



class Interpolator 
{

public:

  typedef ClusterHelper::vector_type vector_t;

  Interpolator(const alps::Parameters& parms);
  virtual ~Interpolator() {}
  virtual void local_green_function() = 0;
  void interpolate();

protected:

  boost::shared_ptr<DCATransformer> transform_ptr;
  const double mu;
  const double delta;
  const double delta_green;
  const std::string name;
  const unsigned int n_zone;
  const unsigned int n_k;
  
  virtual void interpolate_spectrum(const unsigned int counter, const unsigned int, 
				    const std::vector<vector_t>& kvec,
				    const double offset, const int sign) = 0;

  void bare_dispersion(const unsigned int counter, const unsigned int, 
		       const std::vector<vector_t>& kvec,
		       const double offset, const int sign);

  virtual void interpolate_1D(const unsigned int) = 0;
  virtual void interpolate_2D(const unsigned int) = 0;
  virtual void interpolate_3D(const unsigned int) = 0;
 
};




class ParaInterpolator : public Interpolator 
{
public:
  
  typedef std::map<double, std::vector<std::complex<double> > > selfenergy_t;

  ParaInterpolator(const alps::Parameters& parms);

  virtual void local_green_function();
  
private:
  
  virtual void interpolate_spectrum(const unsigned int counter, const unsigned int zone, 
				    const std::vector<vector_t>& kvec,
				    const double offset, const int sign);  

  virtual void dmft_spectrum(const unsigned int counter, const unsigned int zone,
                             const std::vector<vector_t>& kvec,
                             const double offset, const int sign);
  
  virtual void interpolate_1D(const unsigned int);
  virtual void interpolate_2D(const unsigned int);
  virtual void interpolate_3D(const unsigned int);
  
  std::vector<selfenergy_t> selfenergies;
  std::vector<std::map<double, std::complex<double> > > averaged_selfenergy;
};


      

class AntiferroInterpolator : public Interpolator 
{
public:

  typedef ClusterHelper::Matrix<std::complex<double> >::Type matrix_t;
  typedef std::map<double, std::vector<matrix_t> > selfenergy_t;
  
  AntiferroInterpolator(const alps::Parameters& parms);
  
  virtual void local_green_function() { ; }

private:
  
  virtual void interpolate_spectrum(const unsigned int counter, const unsigned int zone, 
				    const std::vector<vector_t>& kvec,
				    const double offset, const int sign);  

  virtual void interpolate_1D(const unsigned int);
  virtual void interpolate_2D(const unsigned int);
  virtual void interpolate_3D(const unsigned int);
  
  std::vector<selfenergy_t> selfenergies;
};


template<typename T>
inline T ClusterHelper::linear_interpolation(vector_type k, const std::vector<T>& values, const bool antiferromagnet) const 
{
  k = round(matrix_vector_mult(reciprocal_basis_inv_, k));
  std::vector<vector_type> neighbor;
  std::vector<unsigned int> neighbor_index;
  switch (dimension()) {
    case 1:
      neighbor.resize(2);
      neighbor_index.resize(2);
      for (int i=0; i<2; ++i) 
        neighbor[i].resize(1);
      neighbor[0][0] = floor(k[0]);
      neighbor[1][0] = ceil(k[0]);
      break;
    case 2:
      neighbor.resize(4);
      neighbor_index.resize(4);
      for (int i=0; i<4; ++i) 
        neighbor[i].resize(2);
      for (int x=0; x<=1; ++x) {
        double xt = (x==0) ? floor(k[0]) : ceil(k[0]);
        for (int y=0; y<=1; ++y) {
          double yt = (y==0) ? floor(k[1]) : ceil(k[1]);
          neighbor[y*2+x][0] = xt;
          neighbor[y*2+x][1] = yt;
        }
      }
      break;
    case 3:
      neighbor.resize(8);
      neighbor_index.resize(8);
      for (int i=0; i<8; ++i) 
        neighbor[i].resize(3);
      for (int x=0; x<=1; ++x) {
        double xt = (x==0) ? floor(k[0]) : ceil(k[0]);
        for (int y=0; y<=1; ++y) {
          double yt = (y==0) ? floor(k[1]) : ceil(k[1]);
          for (int z=0; z<=1; ++z) {
            double zt = (z==0) ? floor(k[2]) : ceil(k[2]);
            neighbor[z*2*2+y*2+x][0] = xt;
            neighbor[z*2*2+y*2+x][1] = yt;
            neighbor[z*2*2+y*2+x][2] = zt;
          }
        }
      }
      break;
  }
  for (int i=0; i<neighbor.size(); ++i) {
    vector_type nv = matrix_vector_mult(reciprocal_basis_, neighbor[i]);
    nv = map_back_into_brillouin_zone(nv);
    if (antiferromagnet) {
      if (!inside_reduced_brillouin_zone(nv)) {
        for (int d=0; d<dimension(); ++d) 
          nv[d] -= M_PI;
        nv = map_back_into_brillouin_zone(nv);
      }
      if (!inside_reduced_brillouin_zone(nv)) {
        boost::throw_exception(std::logic_error("something doesn't work in ClusterHelper::interpolate)"));
      }
    }
    neighbor_index[i] = identify_vertex(nv, cluster_momenta_);
  }
  std::vector<T> av(neighbor_index.size());
  for (int i=0; i<neighbor_index.size(); ++i) 
    av[i] = values[neighbor_index[i]];
  for (int d=0; d<dimension(); ++d) {
    double c = k[d] - floor(k[d]) ;
    std::vector<T> t(av.size()/2);
    for (int i=0; i<av.size(); i+=2) {
      t[i/2] = (1-c)*av[i] + c*av[i+1];
    }
    av = t;
  }
  return av[0];
}




template<typename T>
inline T ClusterHelper::interpolation(vector_type k, const std::vector<T>& values, const bool antiferromagnet) const 
{
  k = round(matrix_vector_mult(reciprocal_basis_inv_, k));
  return interpolator(k, dimension()-1, values, antiferromagnet);
}



template<typename T>
inline T ClusterHelper::interpolator(const vector_type k, const int d, const std::vector<T>& values,
                             const bool antiferromagnet) const
{
  if (d<0) {
    unsigned int index = identify_vertex(map_interpolation_vertex(k, antiferromagnet), cluster_momenta_);;
    return values[index];
  }
  else {
    std::map<double, T> pairs;
    vector_type kf = k;
    kf[d] = ceil(kf[d]);
    vector_type k1 = map_interpolation_vertex(kf, antiferromagnet);

    //double c = k[d] - floor(k[d]) ;
    //std::vector<T> t(av.size()/2);
    //for (int i=0; i<av.size(); i+=2) {
    //  t[i/2] = (1-c)*av[i] + c*av[i+1];
    //}
    //av = t;

    do {
      pairs[kf[d]] = interpolator(kf, d-1, values, antiferromagnet);
      kf[d] = kf[d]+1;
    } while(!almost_equal(k1, map_interpolation_vertex(kf, antiferromagnet)));
    return interpolate(k[d], pairs);
  }
}



template<>
inline double ClusterHelper::interpolate(const double k, const std::map<double, double>& values) const
{
  
  if (values.size() == 1) {
    return values.begin()->second;
  }
  else {
    const int N = values.size();
    gsl_interp_accel* acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_akima_periodic, 2*N+1);
    double x[2*N+1];
    double y[2*N+1];
    std::map<double, double>::const_iterator it = values.begin();
    x[2*N] = (it->first)+N;
    y[2*N] = it->second;
    for(unsigned int i=0; it!=values.end(); ++it, i++) {
      x[i] = (it->first)-N;
      y[i] = it->second;
      x[i+N] = it->first;
      y[i+N] = it->second;
    }
    gsl_spline_init(spline, x, y, 2*N+1);
    double val = gsl_spline_eval(spline, k, acc);
    gsl_interp_accel_free (acc);
    gsl_spline_free (spline);
    return val;
  }
}  



template<>
inline std::complex<double> ClusterHelper::interpolate(const double k, const std::map<double, std::complex<double> >& values) const
{
  std::map<double, double> re;
  std::map<double, double> im;
  for (std::map<double, std::complex<double> >::const_iterator it=values.begin(); it!=values.end(); it++) {
    re[it->first] = it->second.real();
    im[it->first] = it->second.imag();
  }
  double dim = interpolate(k, im);
  return std::complex<double>(interpolate(k, re), dim);
}



template<>
inline ClusterHelper::Matrix<std::complex<double> >::Type 
ClusterHelper::interpolate(const double k,
                           const std::map<double, Matrix<std::complex<double> >::Type>& values) const
{
  std::map<double, std::complex<double> > m1;
  std::map<double, std::complex<double> > m2;
  std::map<double, std::complex<double> > m3;
  std::map<double, std::complex<double> > m4;
  for (std::map<double, Matrix<std::complex<double> >::Type >::const_iterator it=values.begin(); it!=values.end(); it++) {
    m1[it->first] = it->second(0,0);
    m2[it->first] = it->second(1,0);
    m3[it->first] = it->second(0,1);
    m4[it->first] = it->second(1,1);
  }
  Matrix<std::complex<double> >::Type  M(2,2);
  M(0,0) = interpolate(k, m1);
  M(1,0) = interpolate(k, m2);
  M(0,1) = interpolate(k, m3);
  M(1,1) = interpolate(k, m4);
  return M;
}
    
#endif
