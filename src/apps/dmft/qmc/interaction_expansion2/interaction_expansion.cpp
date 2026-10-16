/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#include <utility> // Boost.Math 1.76 includes this inside a namespace.
#include <boost/math/constants/constants.hpp>
#include "interaction_expansion.hpp"
#include <ctime>
#include "xml.h"
#include "alps/ngs/make_deprecated_parameters.hpp"

//global variables

frequency_t c_or_cdagger::nm_;
bool c_or_cdagger::use_static_exp_;
unsigned int c_or_cdagger::ntau_;
double c_or_cdagger::beta_;
double *c_or_cdagger::omegan_;
std::complex<double> *c_or_cdagger::exp_iomegan_tau_;



namespace {
//defines a parameter that has an older or a derived alternative only if the input supplies it
template<typename T> void define_if_supplied(alps::params &p, std::string const &name, std::string const &description){
  if(p.supplied(name)) p.define<T>(name, description);
}
}

void InteractionExpansion::define_parameters(alps::params &p)
{
  alps::mcbase::define_parameters(p); //SEED
  p.define<long>("SWEEPS", "number of measurement steps");
  p.define<long>("THERMALIZATION", "number of thermalization steps");
  p.define<long>("MAX_TIME", "time limit in seconds");
  p.define<double>("BETA", "inverse temperature");
  p.define<double>("U", "on-site interaction");
  p.define<double>("ALPHA", "shift of the auxiliary field");
  p.define<int>("MAX_ORDER", 2048, "maximum perturbation order");
  p.define<int>("SITES", 1, "number of sites");
  p.define<int>("N_ORBITALS", 2, "number of orbitals; the number of flavors unless FLAVORS is set");
  p.define<int>("N_MEAS", 200, "number of updates between measurements unless MEASUREMENT_PERIOD is set");
  p.define<int>("RECALC_PERIOD", 5000, "number of updates between recomputations of the inverse Green's function matrix");
  p.define<double>("ALMOSTZERO", 1.e-16, "roundoff threshold");
  p.define<bool>("HISTOGRAM_MEASUREMENT", false, "measure the self-energy in imaginary time with histograms instead of on Matsubara frequencies");
  p.define<bool>("ATOMIC", false, "start from the atomic limit instead of the bare Green's function in INFILE");
  p.define<bool>("DMFT_FRAMEWORK", false, "read INFILE in the layout of the DMFT framework");
  p.define<int>("TAU_DISCRETIZATION_FOR_EXP", 0, "number of imaginary-time points at which exp(i omega_n tau) is tabulated; 0 computes it on the fly");
  if(!p["ATOMIC"].as<bool>()) p.define<std::string>("INFILE", "HDF5 file with the bare Green's function G0");
  else define_if_supplied<std::string>(p, "INFILE", "HDF5 file with the bare Green's function G0");
  //sizes with older names (which take precedence) or derived defaults; see ctint_parameters
  define_if_supplied<int>(p, "FLAVORS", "number of flavors; N_ORBITALS by default");
  define_if_supplied<int>(p, "NMATSUBARA", "older name of N_MATSUBARA");
  if(p.supplied("N_MATSUBARA") || !p.supplied("NMATSUBARA")) p.define<int>("N_MATSUBARA", "number of Matsubara frequencies");
  define_if_supplied<int>(p, "NMATSUBARA_MEASUREMENTS", "number of Matsubara frequencies measured; N_MATSUBARA by default");
  define_if_supplied<int>(p, "N", "older name of N_TAU");
  if(p.supplied("N_TAU") || !p.supplied("N")) p.define<int>("N_TAU", "number of imaginary-time slices");
  define_if_supplied<int>(p, "NSELF", "number of imaginary-time points of the self-energy; 10*N_TAU by default");
  define_if_supplied<int>(p, "MEASUREMENT_PERIOD", "number of updates between measurements; N_MEAS by default");
  define_if_supplied<int>(p, "CONVERGENCE_CHECK_PERIOD", "number of updates between convergence checks; RECALC_PERIOD by default");
  //read by the Fourier transformer and the interaction matrix through alps::Parameters
  define_if_supplied<double>(p, "MU", "chemical potential");
  define_if_supplied<double>(p, "H", "magnetic field");
  define_if_supplied<double>(p, "J", "Hund's coupling");
  define_if_supplied<double>(p, "U'", "interaction between different orbitals");
  define_if_supplied<std::string>(p, "U_MATRIX", "interaction matrix");
  std::vector<std::string> flavor_energies;
  for(auto const &supplied : p.supplied_text())
    if(supplied.first.compare(0, 4, "EPS_")==0 || supplied.first.compare(0, 6, "EPSSQ_")==0)
      flavor_energies.push_back(supplied.first);
  for(auto const &name : flavor_energies)
    p.define<double>(name, "EPS_f: energy of flavor f; EPSSQ_f: its second moment");
}

namespace ctint_parameters {
unsigned int n_flavors(const alps::params &p){
  return p.defined("FLAVORS") ? (int)p["FLAVORS"] : (int)p["N_ORBITALS"];
}
unsigned int n_matsubara(const alps::params &p){
  return p.defined("NMATSUBARA") ? (int)p["NMATSUBARA"] : (int)p["N_MATSUBARA"];
}
unsigned int n_matsubara_measurements(const alps::params &p){
  return p.defined("NMATSUBARA_MEASUREMENTS") ? (int)p["NMATSUBARA_MEASUREMENTS"] : (int)n_matsubara(p);
}
unsigned int n_tau(const alps::params &p){
  return p.defined("N") ? (int)p["N"] : (int)p["N_TAU"];
}
unsigned int n_self(const alps::params &p){
  return p.defined("NSELF") ? (int)p["NSELF"] : (int)(10*n_tau(p));
}
unsigned int measurement_period(const alps::params &p){
  return p.defined("MEASUREMENT_PERIOD") ? (int)p["MEASUREMENT_PERIOD"] : (int)p["N_MEAS"];
}
unsigned int convergence_check_period(const alps::params &p){
  return p.defined("CONVERGENCE_CHECK_PERIOD") ? (int)p["CONVERGENCE_CHECK_PERIOD"] : (int)p["RECALC_PERIOD"];
}
}

InteractionExpansion::InteractionExpansion(const alps::params &parms, int node)
: alps::mcbase(parms,node),
max_order(parms["MAX_ORDER"]),
n_flavors(ctint_parameters::n_flavors(parms)),
n_site(parms["SITES"]),
n_matsubara(ctint_parameters::n_matsubara(parms)),
n_matsubara_measurements(ctint_parameters::n_matsubara_measurements(parms)),
n_tau(ctint_parameters::n_tau(parms)),
n_tau_inv(1./n_tau),
n_self(ctint_parameters::n_self(parms)),
mc_steps((boost::uint64_t)parms["SWEEPS"]),
therm_steps((unsigned int)parms["THERMALIZATION"]),        
max_time_in_seconds((int)parms["MAX_TIME"]),
beta((double)parms["BETA"]),                        
temperature(1./beta),
onsite_U((double)parms["U"]),                        
alpha((double)parms["ALPHA"]),
U(alps::make_deprecated_parameters(parms)),                         
recalc_period(parms["RECALC_PERIOD"]),
measurement_period(ctint_parameters::measurement_period(parms)),
convergence_check_period(ctint_parameters::convergence_check_period(parms)),
almost_zero(parms["ALMOSTZERO"]),
green_matsubara(n_matsubara, n_site, n_flavors),
bare_green_matsubara(n_matsubara,n_site, n_flavors), 
bare_green_itime(n_tau+1, n_site, n_flavors),
green_itime(n_tau+1, n_site, n_flavors),
pert_hist(max_order)
{
  //initialize measurement method
  if (static_cast<bool>(parms["HISTOGRAM_MEASUREMENT"]))
    measurement_method=selfenergy_measurement_itime_rs;
  else
    measurement_method=selfenergy_measurement_matsubara;
  for(unsigned int i=0;i<n_flavors;++i)
    g0.push_back(green_matrix(n_tau, 20));
  //other parameters
  weight=0;
  sign=1;
  step=0;
  start_time=time(NULL);
  measurement_time=0;
  update_time=0;
  thermalized=therm_steps==0?true:false;
  if(!parms["ATOMIC"].as<bool>()) {
    alps::hdf5::archive ar(parms["INFILE"].as<std::string>(),"r");
    if(parms["DMFT_FRAMEWORK"].as<bool>()){
      //read in as green_function
//      std::cerr << "Reading G0 ...";
      bare_green_matsubara.read_hdf5(ar,"/G0");
//      std::cerr << " done.\n";
      
    } else { //plain hdf5
      std::vector<std::complex<double> > tmp(n_matsubara);
//      std::cerr << "Reading G0 ...";
      for(std::size_t j=0; j<n_flavors; j++){
        std::stringstream path; path<<"/G0_"<<j;
        ar>>alps::make_pvp(path.str(),tmp);
        for(std::size_t i=0; i<n_matsubara; i++)
          bare_green_matsubara(i,0,0,j)=tmp[i];
      }
//      std::cerr << " done.\n";
      tmp.clear();
    }
    
    FourierTransformer::generate_transformer(alps::make_deprecated_parameters(parms), fourier_ptr);
    fourier_ptr->backward_ft(bare_green_itime, bare_green_matsubara);
  } else {
    for(spin_t flavor=0; flavor<n_flavors; ++flavor) 
      for(site_t site1=0; site1<n_site; ++site1) 
        for(site_t site2=0; site2<n_site; ++site2) 
          for(itime_index_t t=0; t<=n_tau; ++t)
            bare_green_itime(t, site1, site2, flavor)=-0.5;
    for(spin_t flavor=0; flavor<n_flavors; ++flavor) 
      for(site_t site1=0; site1<n_site; ++site1) 
        for(site_t site2=0; site2<n_site; ++site2) 
          for(unsigned int k=0; k<n_matsubara; ++k) 
            bare_green_matsubara(k, site1, site2, flavor)=std::complex<double>(0, -beta/((2*k+1)*boost::math::constants::pi<double>()));
    //fourier transform of -1/2
  }
  //initialize the simulation variables
  initialize_simulation(parms);
  if(node==0) {print(std::cout);}
  vertex_histograms=new simple_hist *[n_flavors*n_flavors];
  vertex_histogram_size=100;
  for(unsigned int i=0;i<n_flavors*n_flavors;++i){
    vertex_histograms[i]=new simple_hist(vertex_histogram_size);
  }
  c_or_cdagger::initialize_simulation(parms);
  
  if(n_site !=1) throw std::invalid_argument("you're trying to run this code for more than one site. Do you know what you're doing?!?");
}



void InteractionExpansion::update()
{
  for(std::size_t i=0;i<measurement_period;++i){
    step++;
    interaction_expansion_step();                
    if(vertices.size()<max_order)
      pert_hist[vertices.size()]++;
    if(step % recalc_period ==0)
      reset_perturbation_series();
  }
}
void InteractionExpansion::measure(){
  measure_observables();
}



double InteractionExpansion::fraction_completed() const{
  //check for error convergence
  if (!thermalized) 
    return 0.;
  if(time(NULL)-start_time> max_time_in_seconds){
    std::cout<<"we ran out of time!"<<std::endl;
    return 1;
  }
  return ((step-therm_steps) / (double) mc_steps);
}



///do all the setup that has to be done before running the simulation.
void InteractionExpansion::initialize_simulation(const alps::params &parms)
{
  weight=0;
  sign=1;
  //set the right dimensions:
  for(spin_t flavor=0;flavor<n_flavors;++flavor)
    M.push_back(inverse_m_matrix());
  vertices.clear();
  pert_hist.clear();
  //initialize ALPS observables
  initialize_observables();
  green_matsubara=bare_green_matsubara;
  green_itime=bare_green_itime;
}



void c_or_cdagger::initialize_simulation(const alps::params &p)
{
  beta_=p["BETA"];
  nm_=ctint_parameters::n_matsubara_measurements(p);
  omegan_ = new double[nm_];
  for(unsigned int i=0;i<nm_;++i) {
    omegan_[i]=(2.*i+1.)*boost::math::constants::pi<double>()/beta_;
  }
  if((int)p["TAU_DISCRETIZATION_FOR_EXP"]!=0) {
    ntau_=p["TAU_DISCRETIZATION_FOR_EXP"];
    use_static_exp_=true;
    exp_iomegan_tau_=new std::complex<double> [2*nm_*ntau_];
    if(exp_iomegan_tau_==0){throw std::runtime_error("not enough memory for computing exp!"); }
    std::cout<<"starting computation of exp values for measurement"<<std::endl;
    for(unsigned int i=0;i<ntau_;++i){
      double tau=i*beta_/(double)ntau_;
      for(unsigned int o=0;o<nm_;++o)
        exp_iomegan_tau_[2*nm_*i + o] = std::complex<double>(cos(omegan_[o]*tau), sin(omegan_[o]*tau));
      for(unsigned int o=0;o<nm_;++o)
        exp_iomegan_tau_[2*nm_*i + nm_ + o] = std::complex<double>(cos(omegan_[o]*tau), -sin(omegan_[o]*tau));
    }
    std::cout<<"done exp computation."<<std::endl;
  } else {
    use_static_exp_=false;
  }
}
