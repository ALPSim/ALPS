/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1999-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Fabian Stoeckli <fabstoec@student.ethz.ch>
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

#ifndef ALPS_APPLICATIONS_MC_SPIN_ABSTRACTSPINMC_H_
#define ALPS_APPLICATIONS_MC_SPIN_ABSTRACTSPINMC_H_

#include <alps/scheduler.h>
#include <alps/alea.h>
#include <boost/property_map/vector_property_map.hpp>

#include <alps/model/sign.h>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include "matrices.h"

typedef alps::uint64_t uint64_t;

struct update_info_type 
{
  std::size_t clustersize;
  double m2_;
};

 
template<class MAT>
class AbstractSpinSim : public alps::scheduler::LatticeMCRun<> {
public : 
  AbstractSpinSim(const alps::ProcessList&,const alps::Parameters&,int);
  ~AbstractSpinSim();  
 
  bool change_parameter(const std::string& name, const alps::StringValue& value);
  void dostep();
  void save(alps::ODump& dump) const;
  void load(alps::IDump& dump);
  bool is_thermalized() const;
  double work_done() const;

  virtual update_info_type do_update()=0;
  virtual void do_measurements(update_info_type)=0;

protected: 
  double beta_;
  uint64_t sweeps_;
  TinyVector<double,MAT::dim> h_;
  TinyVector<double,MAT::dim> h_normalized;
  bool has_magnetic_field_;
  uint64_t sweeps_done_;
  uint64_t thermalization_;
  double thermalization_fraction_;
  uint64_t thermalization_sweeps_;
  
  bool cluster_updates_;
  bool general_case_;
  typedef boost::property_map<graph_type,alps::bond_type_t>::type  bond_type_map;
  typedef boost::property_map<graph_type,alps::site_type_t>::type  site_type_map;
  alps::property_map<alps::bond_type_t,graph_type,int>::type bond_type;

  boost::vector_property_map<double,site_type_map> spinfactor_; /* S */
  boost::vector_property_map<MAT,site_type_map> selfinteraction_; /* D */
  boost::vector_property_map<MAT,bond_type_map> couplings_; /* J */ 
  boost::vector_property_map<double,bond_type_map> couplings_det; /* |J| */
 
  bool use_error_limit;
  bool ferromagnetic_;
  bool antiferromagnetic_;
  bool quantum_convention_;
  double g_;
  int print_sweeps_;
};  

template<class MAT>
AbstractSpinSim<MAT>::AbstractSpinSim(const alps::ProcessList& w, 
        const alps::Parameters& myparms,int n)
  : alps::scheduler::LatticeMCRun<>(w,myparms,n),
  beta_(parms.defined("beta") ? double(parms["beta"]) : 1./double(parms.required_value("T"))),
  sweeps_(static_cast<uint64_t>(parms.required_value("SWEEPS"))),
  h_(parms.defined("h") ? TinyVector<double,MAT::dim>(parms["h"],parms) : TinyVector<double,MAT::dim>(0.0)),
  sweeps_done_(0),
  thermalization_(parms.value_or_default("THERMALIZATION",sweeps_/10)),
  thermalization_fraction_(0.),
  thermalization_sweeps_(0),
  general_case_(false),
  bond_type(alps::get_or_default(alps::bond_type_t(),graph(),0)),
  spinfactor_(boost::get(alps::site_type_t(),graph())),
  selfinteraction_(boost::get(alps::site_type_t(),graph())),
  couplings_(boost::get(alps::bond_type_t(),graph())),
  couplings_det(boost::get(alps::bond_type_t(),graph())),
  ferromagnetic_(true),
  antiferromagnetic_(true),
  quantum_convention_(parms.defined("CONVENTION") && parms["CONVENTION"]=="quantum"),
  g_(parms.defined("g") ? alps::evaluate<double>(parms["g"],parms) : 1.),
  print_sweeps_(parms.value_or_default("PRINT_SWEEPS",0))
{
  // 2009-01-27: wistaria@comp-phys.org
  // moved from the inilializaton list to avoid 'internal compiler error' in gcc 4.3.x
  h_normalized = h_*(1./std::sqrt(h_.get_length2()));
  has_magnetic_field_ = h_normalized.get_length2()>1.e-10;

  if (inhomogeneous())
    boost::throw_exception(std::runtime_error("Disordered lattices not supported by the classical Monte Carlo program.\n"));

  std::map <int,bool> bond_type_visited;
  std::map <int,bool> site_type_visited;

  if (parms.defined("ERROR_VARIABLE"))
    if (parms.defined("ERROR_LIMIT"))
      use_error_limit = true;
    else {
      use_error_limit = false;
    }
  else
    use_error_limit = false; 

  if (has_magnetic_field_ && quantum_convention_) 
    boost::throw_exception(
      std::runtime_error("Cannot use quantum convention with nonzero magnetic field ."));


  if (parms["UPDATE"]=="cluster" && has_magnetic_field_) 
    {boost::throw_exception(std::runtime_error("Illegal update type " 
               + std::string(parms["UPDATE"]) 
               + " with magnetic field not set to zero."));}

  std::string J_string;
  if (parms.defined("J"))
    J_string = parms["J"];
  else
    J_string = "1.0";

  for (site_iterator si = sites().first; si!=sites().second; si++) {
    int type = site_type(*si);
    if (!site_type_visited[type]) {
      // set the values for the spin factors
      double spinvalue;
      std::string spin = "S"+boost::lexical_cast<std::string, int>(type);
      if (!parms.defined(spin))
        spin="S";
      if (!parms.defined(spin))
        spinvalue = quantum_convention_ ? 0.5 : 1.;
      else 
        spinvalue = alps::evaluate<double>(parms[spin],parms);
      spinfactor_[*si]=std::sqrt(spinvalue*(spinvalue+(quantum_convention_ ? 1. : 0.)));

      // set the values for the diagonal matrix terms
      site_type_visited[type]=true;
      std::string diag="D"+boost::lexical_cast<std::string, int>(type);
      std::string diag_string;
      if (!parms.defined(diag)) {
        diag="D";
      }
      if (!parms.defined(diag))
        selfinteraction_[*si].setMatrix(0.0);
      else {
        // parse input string to obtain values for the diagonal term 
        diag_string = parms[diag];
        selfinteraction_[*si].setMatrix(diag_string,parms);
        general_case_=true;
      }
      selfinteraction_[*si] = selfinteraction_[*si]*spinvalue*spinvalue;
      
      if (!selfinteraction_[*si].is_symmetric())
         boost::throw_exception(std::runtime_error("Invalid matrix for D"));
    }
  }

  typedef boost::tuple<int,int,int> bond_tuple_type;
  std::map<bond_tuple_type,int> new_bond_type;
  int num_bond_types=0;


  for (bond_iterator b = bonds().first; b!=bonds().second ;++b) {
    bond_tuple_type this_bond(bond_type[*b],site_type(source(*b)),site_type(target(*b)));
    std::map<bond_tuple_type,int>::const_iterator found = new_bond_type.find(this_bond);
    if (found!=new_bond_type.end())
      bond_type[*b] = found->second;
    else {
      new_bond_type[this_bond] = num_bond_types;
      std::string J_name = parms.value_or_default("J"+
               boost::lexical_cast<std::string,int>(bond_type[*b]),J_string);
      MAT this_J;
      this_J.setMatrix(J_name,parms);
      bond_type[*b] = num_bond_types;
      this_J = this_J * (quantum_convention_ ? 1. : -1.);
      couplings_[*b]=this_J;
      double det = this_J.det();
      couplings_det[*b] = det;
      if (det>0.)
        ferromagnetic_=false;
      else if (det<0.)
        antiferromagnetic_=false;
      num_bond_types++;
    }
  }

  if (parms["UPDATE"]=="cluster") {
    if (MAT::allow_cluster_update)
      cluster_updates_=true;
    else {
      std::cerr << "Invalid input: \n"
                << "Cluster updates are not allowed with the given matrix "
                << "types.\nCluster updates are only allowed for matrices "
                << "reducing to a\nsingle element (these are diagonal matrices "
                << "with always the same\nelement on the diagonal). For this "
                << "case, give one single value\nin the parameter file.\n";
      boost::throw_exception(std::runtime_error("wrong parameter submitted"));
    }
  }
  else if (parms["UPDATE"]=="local")
    cluster_updates_=false;
  else if (parms["UPDATE"]=="")
    cluster_updates_ = !alps::is_frustrated(graph(),couplings_det);
  else
    boost::throw_exception(std::runtime_error("Illegal update type " + std::string(parms["UPDATE"])));

  cluster_updates_ = cluster_updates_ && MAT::allow_cluster_update;

  if ( cluster_updates_) std::cerr << "Using cluster updates\n";
  // check whether frustrated 
  measurements << alps::RealObservable("Energy");
  measurements << alps::RealVectorObservable("Bond-type Energy");
  measurements << alps::RealObservable("Energy Density");
  measurements << alps::RealVectorObservable("Bond-type Energy Density");
  measurements << alps::RealObservable("Energy^2");
  measurements << alps::RealObservable("beta * Energy / sqrt(N)");
  measurements << alps::RealObservable("(beta * Energy)^2 / N");
  measurements << alps::RealObservable("|Magnetization|");
  measurements << alps::RealObservable("Magnetization along Field");
  measurements << alps::RealObservable("Magnetization^2");
  measurements << alps::RealObservable("E.Magnetization^2");
  measurements << alps::RealObservable("Magnetization^4");
  measurements << alps::RealObservable("E.Magnetization^4");
  measurements << alps::RealObservable("Susceptibility");
  if(cluster_updates_)
    measurements << alps::RealObservable("Cluster size");
  if (is_bipartite()) {
    measurements << alps::RealObservable("|Staggered Magnetization|");
    measurements << alps::RealObservable("Staggered Magnetization^2");
  }
  
}



template<class MAT>
AbstractSpinSim<MAT>::~AbstractSpinSim() 
{
// nothing to be done
}  
 
// Check if simulation is finished
// take the measurement that is closer to be finished
// ie. the amount of work done is max(work_done(sweeps),work_done(error_limit)).
template<class MAT>
double AbstractSpinSim<MAT>::work_done() const {
  if (!is_thermalized())
    return 0.;
  double wd_error = 0.0;
  double wd_sweeps = 0.0;
  if (use_error_limit) {
    alps::RealObservable myObs;
#ifndef BOOST_NO_EXCEPTIONS
    try {
#endif
      myObs = measurements.template get<alps::RealObservable>(parms["ERROR_VARIABLE"].c_str());
#ifndef BOOST_NO_EXCEPTIONS
    }
    catch (std::exception& e) {
      std::cerr << "Caught exception: " << e.what() << "\n";
      std::exit(-5);
   }
#endif
    if (myObs.count() == 0)
      return 0.;

    double current_error = myObs.error();
    double error_limit = parms["ERROR_LIMIT"];
    double tmp = error_limit/current_error;
    wd_error = tmp*tmp;
  }
  double good_sweeps = sweeps_done_-thermalization_sweeps_;
  wd_sweeps = (good_sweeps)/double(sweeps_);

  if (wd_sweeps >= 1.0)
    return wd_sweeps;

  if (wd_error > wd_sweeps) return wd_error; 
  else return wd_sweeps;
}

// It is desirable to be able to increase the number of SWEEPS
template<class MAT>
bool AbstractSpinSim<MAT>::change_parameter(const std::string& name, const alps::StringValue& value) {
  uint64_t new_sweeps = 0;
  if(name=="SWEEPS")
    new_sweeps = uint64_t(value);
  // Is it sensible to do it ?
  if(new_sweeps > 0) {
    sweeps_ = new_sweeps;
    return true;
   }
  // Otherwise we cannot do it
  return false;
}

// Do a MC step - and do measurements (if thermalized)
// Do a MC step only sweeps_ times.
template<class MAT>
void AbstractSpinSim<MAT>::dostep() {
  ++sweeps_done_;
  update_info_type s= do_update();
  if (!is_thermalized()) {
    ++thermalization_sweeps_;
    thermalization_fraction_+=double(s.clustersize)/double(num_sites());
  }
  else {
// Is there any good reason to discard sweeps beyond the required ones ?
// if (sweeps_done_ < sweeps_) {
    do_measurements(s);
  }
}
  
// Check if simulation is thermalized
template<class MAT>
bool AbstractSpinSim<MAT>::is_thermalized() const {
  return (thermalization_fraction_ >= thermalization_); 
}


// save simulation
template<class MAT>
void AbstractSpinSim<MAT>::save(alps::ODump& dump) const {
  dump << sweeps_done_ << thermalization_fraction_ << thermalization_sweeps_;
}

// load simulation
template<class MAT>
void AbstractSpinSim<MAT>::load(alps::IDump& dump) {
  dump >> sweeps_done_ >> thermalization_fraction_ >> thermalization_sweeps_;
  if(where.empty())
    measurements.compact();
}

#endif
