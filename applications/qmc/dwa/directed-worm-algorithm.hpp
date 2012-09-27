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

#ifndef ALPS_APPLICATIONS_DIRECTED_WORM_ALGORITHM_HPP
#define ALPS_APPLICATIONS_DIRECTED_WORM_ALGORITHM_HPP


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

#include <alps/alea.h>
#include <alps/osiris/comm.h>

#include "numeric.hpp"
#include "worldline.hpp"

#include "../qmc.h"


namespace alps {
namespace applications {


class directed_worm_algorithm  
  : public QMCRun<>     
{
public:
  typedef  boost::uint64_t                count_type;
  typedef  unsigned short                 State;
  typedef  double                         Time;
  typedef  site_descriptor                Site;
  typedef  Worldline::Kink                Kink;
  typedef  Worldline::Line                Line;
  typedef  Worldline::Lines               Lines;
  typedef  Worldline::LineIterator        LineIterator;
  typedef  Worldline::LineConstIterator   LineConstIterator;
  typedef  Worldline::LinesIterator       LinesIterator;
  typedef  Worldline::LinesConstIterator  LinesConstIterator;

  directed_worm_algorithm
    ( const alps::ProcessList &  processes
    , const alps::Parameters  &  parameters
    , int                        processnode
    );

  static void print_copyright  (std::ostream & out);
  void        print_simulation (std::ostream & out);

private:
  /// private member functions 

  // I/O
  ////void  save(alps::ODump & dump) const { worldline.save(dump); } 
  ////void  load(alps::IDump & dump)       { worldline.load(dump); }

  // regarding simulation backbone (ESSENTIAL)
  void  start();          
  void  dostep();                         // Check this out to see how the current diagonal configuration is being updated.

  // regarding simulation interprocess (ESSENTIAL)
  bool   wormhead_propagation_until_halted_by_wormtail (State const & wormpair_state);
  void   insert_jump_or_bounce        (double const onsite_energy_relative_);
  void   delete_relink_jump_or_bounce ();

  // regarding simulation performance 
  bool   is_thermalized()  const         {  return _sweep_counter >= _total_thermalization_sweeps;  }
  double work_done()       const         {  return (is_thermalized() ? (_sweep_counter - _total_thermalization_sweeps)/double(_total_sweeps) : 0.);  }
  double probability_worm_insertion()  const  {  return (1. - static_cast<double>(_sweep_failure_counter)/_sweep_counter);         }
  double probability_bounce()          const  {  return (static_cast<double>(_propagation_failure_counter)/_propagation_counter);  }

  // regarding model
  void initialize_hamiltonian();
  void print_hamiltonian(std::ostream & out);

  std::vector<double> onsite_hamiltonian(const site_descriptor & site);
  inline boost::multi_array<double,4>  
                      bond_hamiltonian(const bond_descriptor & bond);
  boost::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double> > 
                      bond_ladder_hamiltonians(const bond_descriptor & bond);

  inline double onsite_energy (Site const site_, State const state_) const  { return onsite_matrix[site_type(site_)][state_]; }     
  inline double onsite_energy_relative (Kink const & kink_, bool const forward_) const;
  inline double hopping_energy (bond_descriptor const & bond_, State const targetstate_, bool const increasing_) const;

  std::vector<double> onsite_energies (std::vector<int> const & states_) const
  {
    std::vector<double> _onsite_energies;
    _onsite_energies.reserve(num_sites());
    for (Site site=0; site<num_sites(); ++site)
      _onsite_energies.push_back(onsite_energy(site,states_[site]));
    return _onsite_energies;
  }

  // regarding lattice
  std::vector<double> position(const site_descriptor & site) 
  {
    using alps::numeric::operator-;
    using boost::numeric::operators::operator-;
    using boost::numeric::operators::operator/;
      
    return coordinate(site) - alps::applications::numeric::vector_cast<double>(lattice().extent() - 1)/2.;
  }

  int idx(const std::vector<int> & idx_vector) const
  {
    int _index = idx_vector[0];
    for (int i=0; i<dimension()-1; ++i)
      _index += (idx_vector[i+1]*lattice().extent()[i]);
    return _index;
  } 

  std::vector<int> nearest_neighbor_idx(const site_descriptor & site) const 
  {
    std::set<std::vector<int> > _unit_vectors;
    for (neighbor_bond_iterator it=neighbor_bonds(site).first; it != neighbor_bonds(site).second; ++it)
      _unit_vectors.insert(alps::applications::numeric::iround(alps::numeric::abs(bond_vector(*it))));
    std::vector<int> _nearest_neighbor_idx;
    _nearest_neighbor_idx.reserve(_unit_vectors.size());
    for (std::set<std::vector<int> >::iterator it=_unit_vectors.begin(); it!=_unit_vectors.end(); ++it)
      _nearest_neighbor_idx.push_back(idx(*it));
    std::sort(_nearest_neighbor_idx.begin(), _nearest_neighbor_idx.end());
    return _nearest_neighbor_idx;
  }

  std::vector<int> periodicize (const std::vector<int> & idx_vector_) const
  {
    std::vector<int> _vector = idx_vector_;
    for (int i=0; i<dimension(); ++i)
    {
      if (_vector[i] < 0)
        _vector[i] = lattice().extent()[i] + _vector[i];
    }
    return _vector;
  }

  std::vector<int> winding (const std::vector<int> & idx_vector_)  const
  {
    std::vector<int> _vector = idx_vector_;
    for (int i=0; i<dimension(); ++i)
    {
      if (_vector[i] == lattice().extent()[i]-1)
        _vector[i] = -1;
      else if (_vector[i] == -(lattice().extent()[i]-1))
        _vector[i] = 1;
    }
    return _vector;
  }

  // regarding worldline
  void print_worldline(std::ostream & out) const; 

  // regarding measurements
  void initialize_measurements(); 
  void perform_diagonal_measurements();
  void print_measurements(std::ostream & out);

  // regarding on fly measurements
  void reinitialize_on_fly_measurements();

  // regarding random numbers
  bool  randombool()   {  return random() < 0.5;          }  
  Site  randomsite()   {  return random() * num_sites();  }
  Time  randomtime()   {  return random();                }


  // private member objects

  // regarding MC simulation 
  count_type  _sweep_counter;
  count_type  _sweep_failure_counter;
  count_type  _propagation_counter;
  count_type  _propagation_failure_counter;

  count_type  _total_thermalization_sweeps;
  count_type  _total_sweeps;
  count_type  _sweep_per_measurement;

  // regarding model
  std::vector<State> state_minimum;         // arg: site site type
  std::vector<State> state_maximum;         // arg: site site type

  std::vector<std::vector<double> > onsite_matrix;               // onsite matrix      - arg1: sitetype, arg2: state  
  std::vector<std::vector<double> > site_oneup_matrix;           // site 1-up matrix   - arg1: sitetype, arg2: state
  std::vector<std::vector<double> > site_onedown_matrix;         // site 1-down matrix - arg1: sitetype, arg2: state

  std::vector<double> bond_strength_matrix;                      // bond strength matrix - arg: bondtype

  // regarding disorder
  std::vector<double> percentage_error_hopping_t;
  std::vector<double> percentage_error_onsite_U;
  std::vector<double> percentage_error_onsite_mu;

  // regarding experiment
  bool is_experiment;

  static const double hbar = 1.05457148;
  static const double amu  = 1.66053886;

  std::vector<double> phase;

  // regarding worldline
  Worldline            worldline;
  std::pair<Kink,Kink> wormpair;

  // regarding measurements
  bool _measure_now;

  // regarding on-fly measurements
  bool _on_fly_measurements_reinitialized;

  std::vector<double>  green_tof;
  std::vector<double>  green;
  std::vector<int>     winding_number;
};


} // ending namespace applications
} // ending namespace alps


#endif
