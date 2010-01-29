/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>
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

#ifndef ALPS_APPLICATIONS_WORM_RUN_H
#define ALPS_APPLICATIONS_WORM_RUN_H

#define P_REMOVE 0.5

#include <alps/scheduler.h>
#include <alps/lattice.h>
#include <alps/expression.h>

#include "WKink.h"
#include "random.h"
#include "../qmc.h"


#include <algorithm>

#include <boost/multi_array.hpp>
#include <boost/property_map/vector_property_map.hpp>

//#define Print_config
//#define Print_steps
//#define Print_spins
//#define Print_update
//#define Print_detail
//#define Print_weight

//#define CHECK_OFTEN 

//#define SIMPLE

// non-local interactions?

// data structure for kinks, default is list
// #define USE_VECTOR
// #define USE_SET

struct subinterval_info{ time_struct start_time;
                         time_struct end_time;
                         double delta_t;
                         double delta_e;
                         double integrated_time;
                         double integrated_weight; };

class WRun : public QMCRun<> {
public:

  static void print_copyright(std::ostream&);
  WRun(const alps::ProcessList&,const alps::Parameters&,int);

  void save(alps::ODump&) const;       
  void load(alps::IDump&);       
  void dostep();
  bool is_thermalized() const;
  double work_done() const;

private:

  //- types --------------------------------------------------

  // model specification
  int number_of_bond_types; 
  std::vector<uint8_t> site_state; 

  std::vector<std::pair<uint8_t,uint8_t> > site_type_for_bond_type_;
  boost::vector_property_map<int>          original_bond_type;
  
  std::vector< std::vector<double> > matrix_element_raise_;
  std::vector< std::vector<double> > matrix_element_lower_;

  std::vector< std::vector<double> > site_matrix;
  std::vector<double>                hopping_matrix;

  std::map<int, boost::multi_array<double,2> > diagonal_matrix_element;
  
  // state type, quantum nunmbers
  typedef uint8_t state_type;

  // kinks, wormheads, iterators
  typedef Kink<state_type> kink_type;  
#if defined ( USE_VECTOR )
  typedef std::vector<kink_type> kinklist_type;
#elif defined ( USE_SET )
  typedef std::set<kink_type> kinklist_type;
#else
  typedef std::list<kink_type> kinklist_type;
#endif


  typedef std::vector<kinklist_type> kinkvector_type;
  typedef ::cyclic_iterator<kinklist_type> cyclic_iterator;
  typedef cyclic_iterator::base_iterator iterator;

  typedef Wormhead<kinklist_type,graph_type> wormhead_type;

  //- member functions ---------------------------------------

  // worm creation/annihilation
  inline bool create_worm();
  inline bool annihilate_worm();
  int64_t make_worm();

  // update methods
  void shift_kink(wormhead_type&,wormhead_type&);
  void insert_jump(wormhead_type&,wormhead_type&,int,int);
  void remove_jump(wormhead_type&,wormhead_type&,int,int);

  // unperturbed energy of site states, includes neighbor terms
  inline double H0(const site_descriptor&);

  // unperturbed energy differences
  inline double Delta_H0(const site_descriptor&, const state_type&, const state_type&, std::vector<state_type>&);
  inline double Delta_H0(const site_descriptor&, const state_type&, const state_type&, std::vector<state_type>&,
                         const site_descriptor&, const state_type&, const state_type&, std::vector<state_type>&, const bond_descriptor&);

  inline void Update_Delta_H0(double& delta_e, const state_type& state1, const state_type& state2,  
                              const state_type& nb_state1, const state_type& nb_state2, 
                              const bond_descriptor& b) {
    if (nonlocal) {
      delta_e -= neighbor_energy(state1, nb_state1, b);
      delta_e += neighbor_energy(state2, nb_state1, b);
      delta_e += neighbor_energy(state1, nb_state2, b);
      delta_e -= neighbor_energy(state2, nb_state2, b);
    }
  }   // Update_Delta_H0

  // calculate integrated weight taking into account all subintervals
  inline double integrated_weight(const double& lambda, const double& time)
    {
#ifdef SIMPLE
      if(lambda!=0.)
        return 0.;
      else
        return time;
#else
      if(lambda==0. || fabs(lambda*time)<1e-10)
        return time;
      if ((-lambda*time > log_numeric_limits_double ) && is_thermalized())
        boost::throw_exception(std::logic_error("Exceeding std::numeric_limits<double>::max() in integrated weight"));        
      return 1./lambda*(1.-std::exp(-lambda*time));
#endif
    }   // WRun::integrated_weight

  // explore all subintervals due to non-local interaction
  void traverse_subintervals(wormhead_type&);
  void traverse_subintervals(const time_struct&, const time_struct&, const double&, 
                             const site_descriptor&,
                             const state_type&, const state_type&, double);

  // find time steps where neighbor states change
  std::pair<time_struct, time_struct> adjacent_subinterval(const site_descriptor&, const time_struct&,
                                                                      std::vector<state_type>&);

  // calculate worm creation probability taking into account all subintervals
  double worm_creation_probability(const time_struct&, const double&, const site_descriptor&,
                                   const state_type&, const state_type&, double);
                                    
  // measurements
  void make_meas();
  void measure_green();
  inline void measure() {
    if(--measurements_done==0) {
      measurements_done=skip_measurements;
      make_meas();
    }
  }

  // adjustment
  int get_particle_number();
  void adjustment();
  bool canonical;
  int number_of_bosons;
  double correction;
  int corrections_upwards;
  int corrections_downwards;
  bool preadjustment_done;
  bool adjustment_done;
  std::vector<int> nob;
  std::string adjust_parameter;

  // checks
  void print_spins();
  void check_spins();

  void start();
  
  //- model related functions --------------------------------

  state_type min_state() { return std::numeric_limits<state_type>::min(); }
  state_type max_state() { return std::numeric_limits<state_type>::max(); }

  // state_type initial_state()                 { return state_type(0); }
  state_type create(state_type s, int=0)     { return s+1; }
  state_type annihilate(state_type s, int=0) { return s-1; }

  inline double creation_matrix_element(const state_type& state, const bool& c, site_descriptor site) { 
    int this_site_type = site_type(site);
    return (c ? matrix_element_raise_[this_site_type][state] * matrix_element_raise_[this_site_type][state] 
            : matrix_element_lower_[this_site_type][state] * matrix_element_lower_[this_site_type][state]);
  }

  inline double onsite_energy(const state_type& state, site_descriptor site) { 
    return site_matrix[inhomogeneous_site_type(site)][state];
  }

  inline double neighbor_energy(const state_type& state1, const state_type& state2, const bond_descriptor& b) {
    return diagonal_matrix_element[bond_type[b]][state1][state2];
  }   

  inline double hopping_matrix_element(const bond_descriptor& b) {
    return hopping_matrix[bond_type[b]];
  }

  // initialization
  void initialize_hamiltonian();
  boost::multi_array<double,4> bond_hamiltonian(const bond_descriptor&);
  std::vector<double> site_hamiltonian(const site_descriptor&);
  void print_hamiltonian();
  void create_observables();

  //- data ---------------------------------------------------

  std::vector<wormhead_type> worm_head;
  std::valarray<double> stat;
  int32_t steps; // steps already done
  int min_number;
  int max_number;
  double eta;
  int32_t thermal_sweeps;
  int skip_measurements;
  int measurements_done;
  bool have_worm;
  bool chain_kappa;
  int num_chains;

  int num_kinks;
  int worms_per_kink;
  double worms_per_update;
  double log_numeric_limits_double;
  
  std::vector<subinterval_info> subinterval;
  bool subinterval_valid;
  int current_head_num;

  double Sign;
  bool nonlocal;
  
  bool use_1D_stiffness ; //@#$br

  unsigned last_id_;

  //- kinks and related data structures ----------------------

  // iterator to first kink (in time)
  cyclic_iterator first_kink(site_descriptor i) { return cyclic_iterator(kinks[i],kinks[i].begin());} 
  state_type initial_state(site_descriptor i) 
  { return kinks[i].empty() ?  initial_state_[i] : kinks[i].rbegin()->state(); }
  kinkvector_type kinks;
  std::vector<state_type> initial_state_;
  std::valarray<double> green;
  std::vector<int> chain_number;
    
  alps::property_map<alps::bond_type_t,graph_type,int>::type bond_type;
  alps::property_map<alps::boundary_crossing_t,graph_type,alps::boundary_crossing>::type boundary_crossing;

  inline void erase_kink(int site, iterator w) 
  {
    kinks[site].erase(w);
  }

  inline iterator insert_kink(int site, iterator w, kink_type kink)
  {
    kinklist_type& l(kinks[site]);
#ifdef USE_SET
      w=l.insert(w,kink);
#else
    if (l.empty() || kink.time() < l.begin()->time())
      w=l.insert(l.begin(),kink);
    else if (w==l.begin())
      w=l.insert(l.end(),kink);
    else
      w=l.insert(w,kink);
#endif
    return w;
  }

  inline iterator move_kink(int site, iterator w, time_struct newtime)
  {
    kinklist_type& l(kinks[site]);
#ifdef USE_SET
    kink_type kink(*w);
    l.erase(w);
    kink.set_time(newtime);
    w=l.insert(kink).first;
#else
    iterator k=w;
    ++k;
    if (w==l.begin() && newtime > l.rbegin()->time()) {
      kink_type kink(*w);
      l.erase(w);
      w=l.insert(l.end(),kink);
    }
    else if (k==l.end() && newtime < l.begin()->time()) {
      kink_type kink(*w);
      l.erase(w);
      w=l.insert(l.begin(),kink);
    }
    w->set_time(newtime);
#endif
    return w;
  }
};

//- Unperturbed energy ----------------------------------------------------------

inline double WRun::H0(const site_descriptor& s1) {
  //
  //  determines the energy of the state on site s1 at time 0 taking into account 
  //  the states of all neighbors at that particular time.
  //  This function does not doublecount the bond energy terms.
  //

  // determine onsite energy
  state_type state1 = initial_state(s1);
  double Result = onsite_energy(state1,s1);

  if (nonlocal) {
    // determine neighbor energy
    state_type state2;
    neighbor_bond_iterator nbi, nbi_end;
    int nb=0;
    for(boost::tie(nbi, nbi_end) = neighbor_bonds(s1); nbi != nbi_end; ++nbi) {
      state2 = initial_state(neighbor(s1, nb));
      Result += 0.5*neighbor_energy(state1, state2,*nbi);
      nb++;
    }  
  }

  return Result;
}   // WRun::H0


#endif
