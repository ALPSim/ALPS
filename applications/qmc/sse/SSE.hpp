/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2009 by Fabien Alet <alet@comp-phys.org>,
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

#ifndef SSE_HPP
#define SSE_HPP

// Here do define
// #define AUTO_TIMES

#include "../qmc.h"
#include "SSE.Classes.hpp"
#include <boost/cstdint.hpp>
#include <boost/multi_array.hpp>
#include <boost/property_map/vector_property_map.hpp>
#include <vector>
#include <algorithm>
#include <utility>
#include <numeric>
#include <algorithm>

using namespace std;
using namespace alps;

class SSE : public QMCRun<> {
public : 
  typedef Vertex<state_type> vertex_type;
  static void print_copyright(std::ostream&);
  
  SSE(const ProcessList&,const Parameters&,int);
  void save(ODump&) const;
  void load(IDump&);
  void dostep();
  bool is_thermalized() const;
  double work_done() const;
  bool change_parameter(const std::string&, const StringValue&);


private : 

  typedef QMCRun<> super_type;
  /************* Worm weights ******************/
  double worm_plus(state_type s,state_type t) 
  { 
    return s< matrix_element_raise_[t].size() ? 
      matrix_element_raise_[t][s] : 0.;
  }

  double worm_minus(state_type s,state_type t) 
  { 
    return s< matrix_element_lower_[t].size() ? 
      matrix_element_lower_[t][s] : 0.;
  }
  
  /******** Simulation parameters ***********************
   * It may be desirable to do more than ~4.3e9 sweeps. *
   * Hence use data type uint64_t for nb_steps and     *
   * steps_done_total                                   *
   *****************************************************/

  double epsilon;
  double total_shift;
  alps::uint32_t cutoff_L; 
  alps::uint32_t number_of_worms_per_sweep;
  alps::uint64_t nb_steps; 
  alps::uint32_t nb_thermalisation_steps; 
  alps::uint32_t each_measurement;
  alps::uint64_t steps_done_total; 
  alps::uint32_t measurements_done;
  double worm_size;
  std::string WHICH_LOOP_TYPE;
  bool NO_WORMWEIGHT;
  int nb_worms_thermalization;
  int count_worms_thermalization;

  int number_of_bond_types; 
  vector<state_type> site_state; 
  vector<state_type> site_state_copy; 

  vector<std::pair<state_type,state_type> > site_type_for_bond_type_;
  vector<alps::uint32_t> number_of_bonds_for_bond_type;

  std::vector<std::vector<double> > matrix_element_raise_;
  std::vector<std::vector<double> > matrix_element_lower_;

  std::map<int,boost::multi_array<double,4> > matrix_element;
  boost::multi_array<double, 3> proba_diagonal;
  boost::multi_array<double, 7> proba_worm;
 
  vector<double> energy_offset;
  vector<vertex_type> operator_string;
  vector<vertex_type> operator_string_copy;
  std::vector<alps::uint32_t> op_indices;
  alps::uint32_t current_number_of_non_identity;
  std::valarray<double> green;
  double worm_weight;
  double matelsum;
  unsigned int initial_site;
  unsigned int initial_level;

  //boost::vector_property_map<int> original_bond_type;
  property_map<alps::boundary_crossing_t,graph_type,alps::boundary_crossing>::type boundary_crossing;
  property_map<alps::bond_type_t,graph_type,unsigned int>::type bond_type;
  
  // sign problem support
  double get_sign();
  boost::multi_array<bool, 5>  matrix_sign;

  
  // measurements support
  std::valarray<double> localint; // integrated local density correlations
  std::valarray<double> localint2; // integrated local density correlations

  /******** SSE.Initialization.cpp ********/

  boost::multi_array<double,4> bond_hamiltonian(const bond_descriptor&);
  void initialize_hamiltonian();
  void print_arrays();
  void determine_bonds_offset();
  void initialize_diagonal_update_probabilities();
  void initialize_simulation();

  /******** SSE.Update.cpp ********/
  void print_operator_string();
  void diagonal_update();
  void link_legs();
  std::pair<alps::uint32_t,state_type> initial_vertex();
  void worm_update();
  inline state_type return_exit_leg(alps::uint32_t,state_type);
  void increase_cutoff_L(alps::uint32_t);
  void do_update();

  /******** SSE.Directed.cpp ********/
  double the_sol[4][4];
  bool check_out(double[4]);
  //bool find_sol(double[4],double[4],double[4],double[4]);
  bool find_minbounce(double[4],double[4],double[4],double[4]);
  bool find_locopt(double[4],double[4]);
  double return_weight_after_flip(alps::uint32_t, state_type[4],state_type,bool);
  double return_weight_before_flip(alps::uint32_t, state_type[4],state_type,bool);
  double return_worm_weight_before(alps::uint32_t,state_type[4],state_type,bool);
  double return_worm_weight_after(alps::uint32_t,state_type[4],state_type,bool);
  //bool calculate_scattering_matrix(alps::uint32_t,state_type[4],bool);
  bool calculate_minbounce(alps::uint32_t,state_type[4],bool);
  bool calculate_locopt(alps::uint32_t,state_type[4],bool);
  bool calculate_all_scattering_matrices(alps::uint32_t,bool);            
  void put_in(alps::uint32_t,state_type[4],state_type);
  void calculate_proba_worm();
  void calculate_heat_bath_matrix(alps::uint32_t,state_type[4],bool);

  /******** SSE.Measurements.cpp ********/
  void create_observables();
  void do_measurements();
  
  std::vector<unsigned int> original_bond_type;
  double worm_abort;
};  

#endif
