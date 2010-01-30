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

/* $Id: factory.h 1101 2004-08-18 18:11:43Z troyer $ */

#ifndef ALPS_TOOL_MAXENT_HPP
#define ALPS_TOOL_MAXENT_HPP

#include <alps/scheduler.h>
#include "maxent_parms.hpp"




class MaxEntFactory : public alps::scheduler::Factory 
{
  alps::scheduler::Task* make_task(const alps::ProcessList&, const boost::filesystem::path&) const;  
  void print_copyright(std::ostream& out) const;
};




class MaxEntHelper : private MaxEntParameters
{
public : 
  
  typedef MaxEntParameters::matrix_type matrix_type;
  typedef MaxEntParameters::vector_type vector_type;
  typedef MaxEntParameters::omega_complex_type omega_complex_type;

  MaxEntHelper(const alps::Parameters& p);

  double omega_coord(const int i) const { return MaxEntParameters::omega_coord(i); }

  double Default(const int i) const { return def_[i]; }  
  const vector_type& Default() const { return def_; }

  double Q(const vector_type& u, const double alpha) const {
    using namespace boost::numeric::ublas;
    vector_type A = transform_into_real_space(u);
    return 0.5*chi2(A)-alpha*entropy(A);
  }
    
  int ndat() const { return MaxEntParameters::ndat(); }

  vector_type transform_into_singular_space(vector_type A) const;
  vector_type transform_into_real_space(vector_type u) const;
  vector_type get_spectrum(const vector_type& u) const;
  matrix_type left_side(const vector_type& u) const;
  vector_type right_side(const vector_type& u) const;
  double step_length(const vector_type& delta, const vector_type& u) const;
  double convergence(const vector_type& u, const double alpha) const;
  double log_prob(const vector_type& u, const double alpha) const;
  double chi_scale_factor(vector_type A, const double chi_sq, const double alpha) const;
  double chi2(const vector_type& u) const;
  double entropy(const vector_type& u) const;
  
  private:

  vector_type def_;
};
  



class MaxEntSimulation : public alps::scheduler::Task, private MaxEntHelper
{

public:
  
  MaxEntSimulation(const alps::ProcessList&, const boost::filesystem::path&);
  ~MaxEntSimulation();
  void write_xml_body(alps::oxstream&, const boost::filesystem::path&, bool writeall) const;
  void dostep();
  vector_type levenberg_marquardt(vector_type u, const double alpha) const;
  vector_type iteration(vector_type u, const double alpha, const double mu) const;

private:

  vector_type alpha;
  const double norm;
  const double hartree;
  std::string name;
  boost::filesystem::path dir;
  std::ofstream spex_str;
  std::ofstream chisq_str;
  std::ofstream avspec_str;
  std::ofstream maxspec_str;
  std::ofstream chispec_str;
  std::ofstream prob_str;
}; 


#endif
