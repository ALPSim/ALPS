/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2006 by Fabien Alet <alet@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
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

/* $Id$ */

#ifndef ALPS_APPLICATIONS_QMC_H
#define ALPS_APPLICATIONS_QMC_H

#include <alps/scheduler/montecarlo.h>
#include <alps/scheduler/measurement_operators.h>
#include <alps/alea.h>

#include <boost/optional.hpp>
#include <boost/assert.hpp>

template <class G=typename alps::scheduler::LatticeModelMCRun<>::graph_type, class StateType = boost::uint8_t>
class QMCRun : public alps::scheduler::LatticeModelMCRun<G> 
             , public alps::MeasurementOperators
{
public : 
  typedef StateType state_type;
  typedef alps::scheduler::LatticeModelMCRun<G> super_type;
  typedef typename super_type::site_iterator site_iterator;
  
  QMCRun(const alps::ProcessList&,const alps::Parameters&,int,bool=false);
protected: 
  double beta; 
  bool is_signed_;
  bool is_spin_model_;
  bool is_charge_model_;
  bool measure_local_density_;
  bool measure_local_magnetization_;
  bool measure_site_type_density_;
  bool measure_local_compressibility_;
  bool measure_site_compressibility_;
  bool measure_correlations_;
  bool measure_structure_factor_;
  bool measure_green_function_;
  bool measure_bond_type_stiffness_;
  bool measure_CBS_order_;
  boost::optional<typename super_type::site_descriptor> measurement_origin_;
  std::vector<unsigned int> distance_mult;
  std::size_t num_site_types_;
  std::size_t num_bond_types_;
  std::valarray<double> site_type_density_;
  std::valarray<double> local;

  unsigned maximum_sitetype;
  state_type maximum_number_of_states;
  std::map<int,int> number_states_for_site_type_;
  std::vector<int>        site_site_type;
  std::vector<state_type> site_number_of_states; 

  std::map<std::string,std::vector<std::vector<double> > > diagonal_matrix_element;
  
  boost::optional<int> restricted_particle_number;
  boost::optional<int> restricted_magnetization;
  
  
  void create_common_observables();
  bool do_common_measurements(double sign, const std::vector<state_type>& local, 
            const std::valarray<double>& local_int=std::valarray<double>());
  void initialize_site_states();
private:
  bool build_diagonal_operator(std::string const& op);
};  

template <class G, class StateType>
QMCRun<G,StateType>::QMCRun(const alps::ProcessList& w, const alps::Parameters& p,int n,bool issymbolic)
  : super_type(w,p,n,issymbolic )
  , alps::MeasurementOperators(p)
  , beta(p.defined("Beta") ? alps::evaluate<double>(p["Beta"],p)
      :  (p.defined("beta") ? alps::evaluate<double>(p["beta"],p)
        : (p.defined("BETA") ? alps::evaluate<double>(p["BETA"],p)
           : (p.defined("T") ? 1./alps::evaluate<double>(p["T"],p)
              : (p.defined("TEMPERATURE") ? 1./alps::evaluate<double>(p["TEMPERATURE"],p)
                 : 1./alps::evaluate<double>(p["temperature"],p)))))),
    is_signed_(this->has_sign_problem()),
    is_spin_model_(false),
    is_charge_model_(false),
    measure_local_density_(false),
    measure_local_magnetization_(false),
    measure_site_type_density_(p.value_or_default("MEASURE[Site Type Density]",false)),
    measure_local_compressibility_(false),
    measure_site_compressibility_(false),
    measure_correlations_(p.value_or_default("MEASURE[Correlations]",false)),
    measure_structure_factor_(p.value_or_default("MEASURE[Structure Factor]",false)),
    measure_green_function_(p.value_or_default("MEASURE[Green Function]",false)),
    measure_bond_type_stiffness_(p.value_or_default("MEASURE[Bond Type Stiffness]",false)),
    measure_CBS_order_(p.value_or_default("MEASURE[CBS Order]",false)),
    num_site_types_(alps::maximum_vertex_type(this->graph())+1),
    num_bond_types_(alps::maximum_edge_type(this->graph())+1)
    
{
  if (p.defined("INITIAL_SITE"))
    measurement_origin_=static_cast<int>(p["INITIAL_SITE"]);
  if (p.defined("RESTRICT_MEASUREMENTS[N]"))
    restricted_particle_number.reset(static_cast<int>(p["RESTRICT_MEASUREMENTS[N]"]));
  if (p.defined("RESTRICT_MEASUREMENTS[Sz]"))
    restricted_magnetization.reset(static_cast<int>(p["RESTRICT_MEASUREMENTS[Sz]"]));
    

  if(beta<0)
    boost::throw_exception(
         std::out_of_range("negative inverse temperature beta is illegal"));
         
}

template <class G, class StateType>
bool QMCRun<G,StateType>::build_diagonal_operator(std::string const& name)
{
  if (diagonal_matrix_element.find(name) != diagonal_matrix_element.end())
    return true;
  std::vector<std::vector<double> > vec(maximum_sitetype + 1);

  std::map<int,int>::const_iterator it = number_states_for_site_type_.begin();
  for (; it != number_states_for_site_type_.end(); ++it) {
    int sitetype = it->first;
    int number_states = it->second;
    if (number_states) {
      alps::SiteOperator term(name);
      boost::multi_array<alps::Expression,2> matrix_symbolic =  
          alps::get_matrix(alps::Expression(),term,this->model().basis().site_basis(sitetype),this->parms);
      for (int i=0;i<number_states;++i)
        for (int j=0;j<number_states;++j) {
          if (!matrix_symbolic[i][j].can_evaluate())
            return false;
          else if (i!=j && alps::evaluate<double>(matrix_symbolic[i][j]))
            return false;
          else if (i==j)
            vec[sitetype].push_back(alps::evaluate<double>(matrix_symbolic[i][j]));
        }
    }
  }

  diagonal_matrix_element[name]=vec;
  return true;
}



template <class G, class StateType>
void QMCRun<G,StateType>::initialize_site_states()
{  
  maximum_number_of_states=0;
  maximum_sitetype = 0;

  for (site_iterator it=this->sites().first; it!=this->sites().second;++it) {
    unsigned int sitetype=alps::scheduler::LatticeModelMCRun<G>::site_type(*it);
    if (sitetype > maximum_sitetype)
      maximum_sitetype = sitetype;
    site_site_type.push_back(sitetype);
    std::map<int,int>::const_iterator found = number_states_for_site_type_.find(sitetype);
    if (found != number_states_for_site_type_.end())
      site_number_of_states.push_back(found->second);
    else {
      // get site basis
      int num_states = this->model().basis().site_basis(sitetype).num_states();
    
      number_states_for_site_type_[sitetype]=num_states;
      if (num_states>maximum_number_of_states)
        maximum_number_of_states = num_states;
      site_number_of_states.push_back(num_states);
    }
  }
  
// #error get rid of the is_charge_model_ and is_spin_model_ below and replace it by the general measurements
  is_charge_model_ = build_diagonal_operator("n");
  is_spin_model_ = build_diagonal_operator("Sz");
   
  // set up the matrices for the general measurements
  std::map<std::string,std::string>::iterator it = this-> average_expressions.begin();
  while ( it != this-> average_expressions.end() ) {
    std::map<std::string,std::string>::iterator next = it;
    ++next;
    if (!build_diagonal_operator(it->second)) {
      std::clog << "Will not measure \"" << it->first << "\" since it is off-diagonal or not a site operator\n";
      this-> average_expressions.erase(it);
    }
    it = next;
  }

  it = this->local_expressions.begin();
  while ( it != this-> local_expressions.end() ) {
    std::map<std::string,std::string>::iterator next = it;
    ++next;
    if (!build_diagonal_operator(it->second)) {
      std::clog << "Will not measure \"" << it->first << "\" since it is off-diagonal or not a site operator\n";
      this-> local_expressions.erase(it);
    }
    it = next;
  }

  std::map<std::string,std::pair<std::string,std::string> >::iterator itp = this-> correlation_expressions.begin();
  while ( itp != this-> correlation_expressions.end() ) {
    std::map<std::string,std::pair<std::string,std::string> >::iterator next = itp;
    ++next;
    if (!build_diagonal_operator(itp->second.first) || !build_diagonal_operator(itp->second.second)) {
      std::clog << "Will not measure \"" << itp->first << "\" since it is off-diagonal\n";
      this-> correlation_expressions.erase(itp);
    }
    itp = next;
  }

  itp = this-> structurefactor_expressions.begin();
  while ( itp != this-> structurefactor_expressions.end() ) {
    std::map<std::string,std::pair<std::string,std::string> >::iterator next = itp;
    ++next;
    if (!build_diagonal_operator(itp->second.first) || !build_diagonal_operator(itp->second.second)) {
      std::clog << "Will not measure \"" << itp->first << "\" since it is off-diagonal\n";
      this-> structurefactor_expressions.erase(itp);
    }
    itp = next;
  }

}


template <class G, class StateType>
void QMCRun<G,StateType>::create_common_observables()
{
  local.resize(this->num_sites());
  if (is_signed_)
    this->measurements << alps::RealObservable("Sign");

  // Initialize this->measurements
  this->measurements << alps::make_observable(alps::RealObservable("Energy"),is_signed_);
  this->measurements << alps::make_observable(alps::RealObservable("Energy Density"),is_signed_);

  if (is_charge_model_) {
    measure_local_compressibility_=this->parms.value_or_default(
               "MEASURE[Local Compressibility]",false);
    measure_site_compressibility_=this->parms.value_or_default(
               "MEASURE[Site Compressibility]",false);
    measure_local_density_=this->parms.value_or_default(
                "MEASURE[Local Density]",false) 
          || measure_local_compressibility_ || measure_site_compressibility_;
  }

   if (is_spin_model_) {
     measure_local_magnetization_=this->parms.value_or_default("MEASURE[Local Magnetization]",false);
  }
 

  alps::RealVectorObservable::label_type sitelabels;
  if (measure_local_density_ || measure_local_compressibility_ || measure_local_magnetization_)
    sitelabels =  this->site_labels();

  alps::RealVectorObservable::label_type sitetypelabels;
  if (measure_site_type_density_)
    for (std::size_t i=0; i< num_site_types_;++i)
      sitetypelabels.push_back(boost::lexical_cast<std::string>(i));

  alps::RealVectorObservable::label_type bondtypelabels;
  alps::RealVectorObservable::label_type bondtypecorrlabels;
  if (measure_bond_type_stiffness_) {
    for (std::size_t i=0; i< num_bond_types_;++i) {
      bondtypelabels.push_back(boost::lexical_cast<std::string>(i));
      for (std::size_t j=0; j< num_bond_types_;++j)
        bondtypecorrlabels.push_back(boost::lexical_cast<std::string>(i) + " -- " 
                                   + boost::lexical_cast<std::string>(j));
    }
  }

  alps::RealVectorObservable::label_type corrlabels;
  if (measure_correlations_ || measure_green_function_ || !this->correlation_expressions.empty()) {
    if (measurement_origin_) {
      corrlabels.clear();
      for (unsigned int j=0;j<this->num_sites();++j)
        corrlabels.push_back(this->coordinate_string(measurement_origin_.get())+" -- " + this->coordinate_string(j));
    }
    else {
      corrlabels = this->distance_labels();
      distance_mult = this->distance_multiplicities();
    }
  }

  
  alps::RealVectorObservable::label_type momentalabels;
  if(measure_structure_factor_ || !this->structurefactor_expressions.empty())
    momentalabels = this->momenta_labels();
    
  this->measurements << alps::make_observable(alps::RealObservable("Stiffness"),is_signed_);
  if (measure_bond_type_stiffness_) {
    this->measurements << alps::make_observable(alps::RealVectorObservable("Bond Type Stiffness",bondtypelabels),is_signed_);
    this->measurements << alps::make_observable(alps::RealVectorObservable("Bond Type Stiffness Correlations",bondtypecorrlabels),is_signed_);
  }
  if (is_spin_model_) {
    this->measurements << alps::make_observable(alps::RealObservable("Magnetization"),is_signed_);
    this->measurements << alps::make_observable(alps::RealObservable("Magnetization Density"),is_signed_);
    this->measurements << alps::make_observable(alps::RealObservable("|Magnetization|"),is_signed_);
    this->measurements << alps::make_observable(alps::RealObservable("|Magnetization Density|"),is_signed_);
    this->measurements << alps::make_observable(alps::RealObservable("Magnetization^2"),is_signed_);
    this->measurements << alps::make_observable(alps::RealObservable("Magnetization Density^2"),is_signed_);
    this->measurements << alps::make_observable(alps::RealObservable("Magnetization^4"),is_signed_);
    this->measurements << alps::make_observable(alps::RealObservable("Magnetization Density^4"),is_signed_);
    this->measurements << alps::make_observable(alps::RealObservable("Susceptibility"),is_signed_);
    if (this->is_bipartite()) {
      this->measurements << alps::make_observable(alps::RealObservable("Staggered Magnetization"),is_signed_);
      this->measurements << alps::make_observable(alps::RealObservable("Staggered Magnetization Density"),is_signed_);
      this->measurements << alps::make_observable(alps::RealObservable("Staggered Magnetization^2"),is_signed_);
      this->measurements << alps::make_observable(alps::RealObservable("Staggered Magnetization Density^2"),is_signed_);
      this->measurements << alps::make_observable(alps::RealObservable("Staggered Magnetization^4"),is_signed_);
      this->measurements << alps::make_observable(alps::RealObservable("Staggered Magnetization Density^4"),is_signed_);
    }
    if(measure_correlations_)
      this->measurements << alps::make_observable(alps::RealVectorObservable("Spin Correlations",corrlabels),is_signed_);
    if(measure_structure_factor_)
      this->measurements << alps::make_observable(alps::RealVectorObservable("Spin Structure Factor",momentalabels),is_signed_);
    if (measure_local_magnetization_)
      this->measurements << alps::make_observable(alps::RealVectorObservable("Local Magnetization",sitelabels),is_signed_);
  }

if(measure_CBS_order_)
{
  this->measurements << alps::make_observable(alps::RealObservable("CBS Order"),is_signed_);
  this->measurements << alps::make_observable(alps::RealObservable("CBS Order^2"),is_signed_);
  this->measurements << alps::make_observable(alps::RealObservable("CBS Order^4"),is_signed_);
}
  
  if (is_charge_model_) {
    this->measurements << alps::make_observable(alps::RealObservable("Density"),is_signed_);
    this->measurements << alps::make_observable(alps::RealObservable("Density^2"),is_signed_);
    if (this->is_bipartite()) {
      this->measurements << alps::make_observable(alps::RealObservable("Checkerboard Charge Order"),is_signed_);
      this->measurements << alps::make_observable(alps::RealObservable("Checkerboard Charge Order^2"),is_signed_);
      this->measurements << alps::make_observable(alps::RealObservable("Checkerboard Charge Order^4"),is_signed_);
    }
    if (measure_local_density_)
      this->measurements << alps::make_observable(alps::RealVectorObservable("Local Density",sitelabels),is_signed_);
    if (measure_site_type_density_)
      this->measurements << alps::make_observable(alps::RealVectorObservable("Site Type Density",sitetypelabels),is_signed_);
    if (measure_site_compressibility_)
      this->measurements << alps::make_observable(alps::RealVectorObservable("Integrated Local Density Correlations",sitelabels),is_signed_);
    if (measure_local_compressibility_)
      this->measurements << alps::make_observable(alps::RealVectorObservable("Local Density * Global Density",sitelabels),is_signed_);
    if(measure_correlations_)
      this->measurements << alps::make_observable(alps::RealVectorObservable("Density Correlations",corrlabels),is_signed_);
    if(measure_structure_factor_)
      this->measurements << alps::make_observable(alps::RealVectorObservable("Density Structure Factor",momentalabels),is_signed_);
  }

  if (measure_green_function_)
    this->measurements << alps::make_observable(alps::RealVectorObservable("Green's Function",corrlabels),is_signed_);

  // set up the matrices for the general measurements
  std::map<std::string,std::string>::iterator it = this-> average_expressions.begin();
  while ( it != this-> average_expressions.end() ) {
    std::map<std::string,std::string>::iterator next = it;
    ++next;
    if (this->measurements.has(it->first)) {
      std::clog << "Will not measure custom average measurement \"" << it->first << "\" since a hard-coded measurement already exists\n";
      this-> average_expressions.erase(it);
    }
    else
      this->measurements << alps::make_observable(alps::RealObservable(it->first),is_signed_);
    it = next;
  }

  alps::RealVectorObservable::label_type bondlabels;
  it = this->local_expressions.begin();
  while ( it != this-> local_expressions.end() ) {
    std::map<std::string,std::string>::iterator next = it;
    ++next;
    if (this->measurements.has(it->first)) {
      std::clog << "Will not measure custom local measurement \"" << it->first << "\" since a hard-coded measurement already exists\n";
      this-> local_expressions.erase(it);
    }
    else if (this->has_bond_operator(it->second)) {
      if (bondlabels.empty())
        bondlabels = this->bond_labels();
      this->measurements << alps::make_observable(alps::RealVectorObservable(it->first,bondlabels),is_signed_);
    }
    else {
      if (sitelabels.empty())
        sitelabels = this->site_labels();
      this->measurements << alps::make_observable(alps::RealVectorObservable(it->first,sitelabels),is_signed_);
    }
    it = next;
  }

  std::map<std::string,std::pair<std::string,std::string> >::iterator itp = this-> correlation_expressions.begin();
  while ( itp != this-> correlation_expressions.end() ) {
   std::map<std::string,std::pair<std::string,std::string> >::iterator next = itp;
    ++next;
    if (this->measurements.has(itp->first)) {
      std::clog << "Will not measure custom correlation measurements \"" << itp->first << "\" since a hard-coded measurement already exists\n";
      this-> correlation_expressions.erase(itp);
    }
    else
      this->measurements << alps::make_observable(alps::RealVectorObservable(itp->first,corrlabels),is_signed_);
    itp = next;
  }

  itp = this-> structurefactor_expressions.begin();
  while ( itp != this-> structurefactor_expressions.end() ) {
   std::map<std::string,std::pair<std::string,std::string> >::iterator next = itp;
    ++next;
    if (this->measurements.has(itp->first)) {
      std::clog << "Will not measure custom correlation measurements \"" << itp->first << "\" since a hard-coded measurement already exists\n";
      this-> structurefactor_expressions.erase(itp);
    }
    else
      this->measurements << alps::make_observable(alps::RealVectorObservable(itp->first,momentalabels),is_signed_);
    itp = next;
  }
}


template <class G, class StateType>
bool QMCRun<G,StateType>::do_common_measurements(double sign, const std::vector<state_type>& state, const std::valarray<double>& localint)
{
  std::vector<std::vector<double> > const& matrix_element_Sz = diagonal_matrix_element["Sz"];
  std::vector<std::vector<double> > const& matrix_element_n = diagonal_matrix_element["n"];
  
  double NbSites=this->num_sites();

  if (is_signed_)
    this->measurements["Sign"] << sign;

  if (is_spin_model_) {
    double sz=0.;
    double ssz=0.;
    for (unsigned int i=0;i<this->num_sites();++i) {
      double x = local[i] = matrix_element_Sz[this->site_type(i)][state[i]];
      sz+=x;        
      if (this->is_bipartite())
        ssz += this->parity(i)*x;
    }
    
    if (restricted_magnetization && *restricted_magnetization != sz)
      return false;
    
    this->measurements["Magnetization"] << sz*sign;
    this->measurements["Magnetization Density"] << sz/NbSites*sign;
    this->measurements["|Magnetization|"] << std::abs(sz)*sign;
    this->measurements["|Magnetization Density|"] << std::abs(sz)/NbSites*sign;
    this->measurements["Magnetization^2"] << sz*sz*sign;
    this->measurements["Magnetization Density^2"] << sz*sz/NbSites/NbSites*sign;
    this->measurements["Magnetization^4"] << sz*sz*sz*sz*sign;
    this->measurements["Magnetization Density^4"] << sz*sz*sz*sz/NbSites/NbSites/NbSites/NbSites*sign;
    this->measurements["Susceptibility"] << sz*sz*beta/NbSites*sign;
    if (this->is_bipartite()) {
      this->measurements["Staggered Magnetization"] << ssz*sign;
      this->measurements["Staggered Magnetization Density"] << ssz/NbSites*sign;
      this->measurements["Staggered Magnetization^2"] << ssz*ssz*sign;
      this->measurements["Staggered Magnetization Density^2"] << ssz*ssz/NbSites/NbSites*sign;
      this->measurements["Staggered Magnetization^4"] << ssz*ssz*ssz*ssz*sign;
      this->measurements["Staggered Magnetization Density^4"] << ssz*ssz*ssz*ssz/NbSites/NbSites/NbSites/NbSites*sign;
    }
    if (measure_local_magnetization_) {
      std::valarray<double> local_sz(local);
      local_sz *= sign;
      this->measurements["Local Magnetization"] << local_sz;
    }
  } 
  
  if (is_charge_model_) {
    if(measure_site_type_density_) {
      site_type_density_.resize(num_site_types_);
      for (int i=0;i<num_site_types_;++i)
        site_type_density_[i]=0.;
    }
    double n=0.;
    double sn=0.;
    for (unsigned int i=0;i<this->num_sites();++i) {
      double nloc = local[i] = matrix_element_n[this->site_type(i)][state[i]];
      n+=nloc;
      if(measure_site_type_density_)
        site_type_density_[this->site_type(i)] += nloc;
      if (this->is_bipartite())
        sn += this->parity(i)*nloc;
      
    }
    
    if (restricted_particle_number && *restricted_particle_number != n)
      return false;

    if (this->is_bipartite()) {
      this->measurements["Checkerboard Charge Order"] << sn/NbSites/NbSites*sign;
      this->measurements["Checkerboard Charge Order^2"] << sn*sn/NbSites/NbSites*sign;
      this->measurements["Checkerboard Charge Order^4"] << sn*sn*sn*sn/NbSites/NbSites/NbSites/NbSites*sign;
    }
    if(measure_site_type_density_) {
      site_type_density_ *= sign/this->num_sites();
      this->measurements["Site Type Density"] << site_type_density_;
    }
    if (measure_local_density_) {
      std::valarray<double> local_density(local);
      local_density *= sign;
      this->measurements["Local Density"] << local_density;
      if (measure_local_compressibility_) {
        local_density *= n/NbSites;
        this->measurements["Local Density * Global Density"] << local_density;
      }
      if (measure_site_compressibility_) {
        BOOST_ASSERT (local_density.size() == localint.size());
        local_density=localint*sign;
        this->measurements["Integrated Local Density Correlations"] << local_density;
      }
    }
    this->measurements["Density"] << n/NbSites*sign;
    this->measurements["Density^2"] << n*n/NbSites*sign;
  }

  // custom average measurements
  
  typedef std::pair<std::string,std::string> string_pair_type;
  BOOST_FOREACH(string_pair_type const& x, this->average_expressions) {
    std::vector<std::vector<double> > const& matrix_element = diagonal_matrix_element[x.second];
    double ave=0.;
    for (unsigned int i=0;i<this->num_sites();++i)
      ave += matrix_element[this->site_type(i)][state[i]];
    this->measurements[x.first] << ave/NbSites*sign;
  }

  // custom local measurements
  BOOST_FOREACH(string_pair_type const& x, this->local_expressions) {
    std::vector<std::vector<double> > const& matrix_element = diagonal_matrix_element[x.second];
    std::valarray<double> local(this->num_sites());
    for (unsigned int i=0;i<this->num_sites();++i)
      local[i] = matrix_element[this->site_type(i)][state[i]];
    local *= sign;
    this->measurements[x.first] << local;
  }

  // custom correlation measurements
  typedef std::pair<std::string,std::pair<std::string,std::string> > string_pair_pair_type;
  BOOST_FOREACH(string_pair_pair_type const& x, this->correlation_expressions) {
    std::vector<std::vector<double> > const& matrix_element_a = diagonal_matrix_element[x.second.first];
    std::vector<std::vector<double> > const& matrix_element_b = diagonal_matrix_element[x.second.second];
    std::valarray<double> corr(this->num_distances());
    std::vector<double> local_a(this->num_sites());
    std::vector<double> local_b(this->num_sites());
    for (int i=0;i<this->num_sites();++i) {
      local_a[i] = matrix_element_a[this->site_type(i)][state[i]];
      local_b[i] = matrix_element_b[this->site_type(i)][state[i]];
    }
    corr=0.;
    for (int i1=0;i1<this->num_sites();++i1)
      for (int i2=0;i2<this->num_sites();++i2) {
        typename super_type::size_type d = this->distance(i1,i2);
        BOOST_ASSERT(d>=0 && d<corr.size());
        corr[this->distance(i1,i2)] += local_a[i1]*local_b[i2];
      }
    for (int i=0;i<corr.size();++i)
      corr[i]*=sign/distance_mult[i];
    this->measurements[x.first] << corr;
  }
  
  // Correlation measurements
  
  if(measure_correlations_) {
    std::valarray<double> corr;
    if(measurement_origin_) {
      corr.resize(this->num_sites());
      int i1=measurement_origin_.get();
      for (int i2=0;i2<this->num_sites();++i2)
        corr[i2] += local[i1]*local[i2];
    }
    else {
      corr.resize(this->num_distances());
      corr=0.;
      for (int i1=0;i1<this->num_sites();++i1)
        for (int i2=0;i2<this->num_sites();++i2) {
          typename super_type::size_type d = this->distance(i1,i2);
          BOOST_ASSERT(d>=0 && d< corr.size());
          corr[this->distance(i1,i2)] += local[i1]*local[i2];
        }
      for (int i=0;i<corr.size();++i)
        corr[i]/=distance_mult[i];
    }
    this->measurements[is_charge_model_ ? "Density Correlations" : "Spin Correlations"] << corr;
  }

  // Structure factor measurements
  
  // custom structure factor measurements
  BOOST_FOREACH(string_pair_pair_type const& x, this->structurefactor_expressions) {
    std::vector<std::vector<double> > const& matrix_element_a = diagonal_matrix_element[x.second.first];
    std::vector<std::vector<double> > const& matrix_element_b = diagonal_matrix_element[x.second.second];
    std::vector<double> str;
    str.clear();
    std::vector<double> local_a(this->num_sites());
    std::vector<double> local_b(this->num_sites());
    for (int i=0;i<this->num_sites();++i) {
      local_a[i] = matrix_element_a[this->site_type(i)][state[i]];
      local_b[i] = matrix_element_b[this->site_type(i)][state[i]];
    }
    for (typename super_type:: momentum_iterator mit=this->momenta().first; mit != this->momenta().second; ++mit) {
      std::complex<double> vala, valb;
      for (typename super_type::site_iterator sit=this->sites().first; sit!=this->sites().second;++sit) {
        double phase = alps::numeric::scalar_product(this->momentum(*mit), alps::scheduler::LatticeModelMCRun<G>::coordinate(*sit));
        std::complex<double> cphase(std::cos(phase), std::sin(phase));
        vala += local_a[*sit] * cphase;
        valb += local_b[*sit] * cphase;
      }
      str.push_back(std::real(std::conj(vala)*valb));
    }
    std::valarray<double> strv(str.size());
    for (int i=0;i<str.size();++i)
      strv[i]=str[i]*sign/this->num_sites();
    this->measurements[x.first] << str;
  }

  if(measure_structure_factor_) {
    std::vector<double> str;
    for (typename super_type:: momentum_iterator mit=this->momenta().first; mit != this->momenta().second; ++mit) {
      std::complex<double> val;
      for (typename super_type::site_iterator sit=this->sites().first; sit!=this->sites().second;++sit) {
        double phase = alps::numeric::scalar_product(this->momentum(*mit), alps::scheduler::LatticeModelMCRun<G>::coordinate(*sit));
        val += local[*sit] * std::complex<double>(std::cos(phase), std::sin(phase));
      }
      str.push_back(std::norm(val));
    }
    std::valarray<double> strv(str.size());
    for (int i=0;i<str.size();++i)
      strv[i]=str[i]*sign/this->num_sites();
    this->measurements[is_charge_model_ ? "Density Structure Factor" : "Spin Structure Factor"] << strv;
  }
  if (measure_CBS_order_)
  {
        std::complex<double> val;
        for (typename super_type::site_iterator sit=this->sites().first; sit!=this->sites().second;++sit)
        {
            double phase = M_PI * (alps::scheduler::LatticeModelMCRun<G>::coordinate(*sit)[0] + alps::scheduler::LatticeModelMCRun<G>::coordinate(*sit)[1]);
            val += local[*sit] * (std::complex<double> ( std::cos(phase), std::sin(phase) ));            
        }
        double CBS_order = (std::abs(val) * std::abs(val) * sign) / this->num_sites();
        this->measurements["CBS Order"] << CBS_order;
        this->measurements["CBS Order^2"] << CBS_order * CBS_order;
        this->measurements["CBS Order^4"] << CBS_order * CBS_order * CBS_order * CBS_order;
  }
  
  return true;
}

#endif
