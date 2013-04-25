/****************************************************************************
*
* ALPS Project Applications: Directed Worm Algorithm  
*
* Copyright (C) 2012 by Lode Pollet      <pollet@phys.ethz.ch>  
*                       Ping Nang Ma     <pingnang@phys.ethz.ch> 
*                       Matthias Troyer  <troyer@phys.ethz.ch>    
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
*
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

#ifndef DIRECTED_WORM_ALGORITHM_NGS_HPP
#define DIRECTED_WORM_ALGORITHM_NGS_HPP


#include <cassert>
#include <ctime>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <set>
#include <vector>
#include <utility>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iterator>

#include <boost/cstdint.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/tuple/tuple.hpp>

#include <alps/numeric/vector_functions.hpp>
#include <alps/numeric/vector_valarray_conversion.hpp>

#include "worldlines.hpp"

#include "../qmc.ngs.h"


class directed_worm_algorithm  
  : public qmcbase<>     
{
public:
  typedef boost::uint64_t                count_type;
  typedef worldlines::location_type      location_type;

  directed_worm_algorithm(alps::hdf5::archive & ar);

  static void print_copyright  (std::ostream & out);
  void        print_simulation (std::ostream & out);

  // i/o
  void save(alps::hdf5::archive & ar) const;
  void load(alps::hdf5::archive & ar);


private:

  /// private member functions 

  // regarding simulation backbone (ESSENTIAL)
  void  update();
  void  measure() {}  // I don't measure every sweep, and I do a checkpoint rather then only measure.      
  void  docheckpoint(bool const & measure_);

  // regarding simulation interprocess (ESSENTIAL)
  bool   wormhead_propagates_till_collision_with_wormtail (unsigned short wormpair_state, std::pair<neighbor_iterator, neighbor_iterator> const & neighbors_, bool const & measure_);
  void   insert_jump_or_bounce        (double onsite_energy_relative_, std::pair<neighbor_iterator, neighbor_iterator> const & neighbors_);
  void   delete_relink_jump_or_bounce (std::pair<neighbor_iterator, neighbor_iterator> const & neighbors_);

  // regarding simulation performance 
  double fraction_completed()          const  {  return static_cast<double>(_sweep_counter)/double(_total_sweeps);  }
  double probability_worm_insertion()  const  {  return (1. - static_cast<double>(_sweep_failure_counter)/_sweep_counter);         }
  double probability_bounce()          const  {  return (static_cast<double>(_propagation_failure_counter)/_propagation_counter);  }

  // regarding model
  void initialize_hamiltonian();
  void print_hamiltonian(std::ostream & out);

  std::vector<double> onsite_hamiltonian(unsigned int site_type);
  inline boost::multi_array<double,4>  
                      bond_hamiltonian(const bond_descriptor & bond);
  boost::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double> > 
                      bond_ladder_hamiltonians(const bond_descriptor & bond);

  inline double onsite_energy (unsigned int site_, unsigned short state_) const  { return onsite_matrix[site_type(site_)][state_]; }     
  double onsite_energy_relative (unsigned int site_, unsigned short state_, bool forward_, bool creation_) const;
  double hopping_energy (bond_descriptor const & bond_, unsigned short targetstate_, bool increasing_) const;

  std::vector<double> onsite_energies (std::vector<unsigned short> const & states_) const
  {
    std::vector<double> _onsite_energies;
    _onsite_energies.reserve(num_sites());
    for (unsigned int site=0; site<num_sites(); ++site)
      _onsite_energies.push_back(onsite_energy(site,states_[site]));
    return _onsite_energies;
  }

  // regarding lookups
  void initialize_lookups();
  void print_lookups(std::ostream & out) const;

  // regarding lattice
  std::vector<double> position(unsigned int site)                   const  { return position_lookup[site]; }
  double              position(unsigned int site, unsigned int dim) const  { return position_lookup[site][dim]; } 

  std::vector<double> component_momenta()                  const  { return component_momenta_lookup; }
  double              component_momentum(unsigned int idx) const  { return component_momenta_lookup[idx]; }
  unsigned int        num_component_momenta()              const  { return num_component_momenta_; }

  // regarding experiment
  double phase(unsigned int site) const  { return finite_tof ? phase_lookup[site] : 0.; }  

  // regarding worldline
  void print_worldline(std::ostream & out) const; 

  // regarding measurements
  void initialize_measurements(); 
  void perform_diagonal_measurements();

  // regarding on fly measurements
  void reinitialize_on_fly_measurements();


  // private member objects

  // regarding MC simulation 
  count_type  _sweep_counter;
  count_type  _sweep_failure_counter;
  count_type  _propagation_counter;
  count_type  _propagation_failure_counter;

  count_type  _total_sweeps;
  count_type  _sweep_per_measurement;

  std::pair<std::time_t, std::time_t>  _simulation_timer;  

  // regarding model
  std::vector<unsigned short> state_minimum;         // arg: site site type
  std::vector<unsigned short> state_maximum;         // arg: site site type

  std::vector<std::vector<double> > onsite_matrix;               // onsite matrix      - arg1: sitetype, arg2: state  
  std::vector<std::vector<double> > site_oneup_matrix;           // site 1-up matrix   - arg1: sitetype, arg2: state
  std::vector<std::vector<double> > site_onedown_matrix;         // site 1-down matrix - arg1: sitetype, arg2: state

  std::vector<double> bond_strength_matrix;                      // bond strength matrix - arg: bondtype

  // regarding lookups
  std::vector<std::vector<double> > position_lookup;
  std::vector<double>               component_momenta_lookup;
  std::vector<double>               phase_lookup;

  // regarding lattice
  bool is_periodic_;
  int  num_component_momenta_;

  std::vector<std::vector<double> > lattice_vector_;

  // regarding experiment
  bool finite_tof;

  // regarding worldline
  worldlines wl;
  wormpair   worm;

  // regarding measurements
  bool measure_only_simulation_speed_;
  bool measure_winding_number_;
  bool measure_local_density_;
  bool measure_local_density2_;
  bool measure_local_magnetization_;
  bool measure_local_magnetization2_;
  bool measure_local_energy_;
  bool measure_green_function_;
  bool measure_momentum_density_;
  bool measure_tof_image_;

  // regarding on-fly measurements
  double  green_onsite;
  double  green_neighbors;
  double  zero_momentum_density;
  double  zero_momentum_density_tof;

  std::vector<double>  green;
  std::vector<double>  green_tof;
  std::vector<double>  momentum_density;
  std::vector<double>  momentum_density_tof;

  // regarding caches (for optimization purpose)
#ifdef HEATBATH_ALGORITHM
  std::vector<location_type>  _neighborlocations_cache;
  std::vector<double>         _cummulative_weights_cache; 
#endif
};
#endif
