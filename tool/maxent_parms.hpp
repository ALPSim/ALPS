/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2010 by Sebastian Fuchs <fuchs@comp-phys.org>
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

#ifndef ALPS_TOOL_MAXENT_PARMS_HPP
#define ALPS_TOOL_MAXENT_PARMS_HPP

#include <alps/config.h> // needed to set up correct bindings
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <alps/osiris.h>
#include <alps/alea.h>
#include <alps/ngs.hpp>
#include <alps/mcbase.hpp>
#include <alps/ngs/signal.hpp>
//#include "maxent.hpp"
#include "default_model.hpp"


class ContiParameters {

public:
  
  typedef boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> matrix_type;
  typedef boost::numeric::ublas::vector<double> vector_type;
  typedef boost::numeric::ublas::vector<std::complex<double> > complex_vector_type;
  typedef std::pair<vector_type, complex_vector_type> omega_complex_type;

//  ContiParameters(const alps::Parameters& p);
  ContiParameters(const alps::params& p);
  
  const DefaultModel& Default() const { return *Default_; }  
  double y(const int i) const { return y_[i]; }
  double cov(const int i,const int j) const { return cov_(i,j); }
  double sigma(const int i) const { return sigma_[i]; }
  double x(const int i) const { return x_[i]; }
  double T() const { return T_; }
  int ndat() const { return ndat_; }
  double K(const int i, const int j) const {
    return K_(i,j);
  }

private:
  
  boost::shared_ptr<DefaultModel> Default_;
  const double T_;
  int ndat_;

protected:

  void setup_kernel(const alps::params& p, const int ntab, const vector_type& freq);
  const int nfreq_;
  vector_type y_,sigma_,x_;
  matrix_type K_,cov_;
  vector_type t_array_;
};




class MaxEntParameters : public ContiParameters
{
public:
  
//  MaxEntParameters(const alps::Parameters& p);
  MaxEntParameters(const alps::params& p);
  
  const vector_type& y() const { return y_; }
  const matrix_type& K() const { return K_; }
  const matrix_type& U() const { return U_; }
  const matrix_type& Vt() const { return Vt_; }
  const matrix_type& Sigma() const { return Sigma_; }
  double omega_coord(const int i) const { return omega_coord_[i]; }
  const vector_type& omega_coord() const { return omega_coord_; }
  double delta_omega(const int i) const { return delta_omega_[i]; }
  int ns() const { return ns_; }
  int nfreq() const { return nfreq_; }

private:

  matrix_type U_;
  matrix_type Vt_;
  matrix_type Sigma_;
  vector_type omega_coord_;
  vector_type delta_omega_;
  int ns_;
};

    
#endif
