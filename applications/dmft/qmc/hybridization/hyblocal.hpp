/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                       Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
 *
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
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
#ifndef LOCAL_CONFIG_HPP
#define LOCAL_CONFIG_HPP

#include<alps/ngs/params.hpp>
#include"hybint.hpp"
#include"hybretintfun.hpp"
#include"hybsegment.hpp"
#include<set>
#include<vector>
//this is the class that handles everything in connection with the local
//impurity operators. It knows about segments, chemical potentials,
//interactions, and so on.

typedef std::set<segment> segment_container_t;
typedef std::map<double,int> state_map;

typedef class local_configuration{
public:
  local_configuration(const alps::params &p, int crank);
  double local_energy(const segment &seg, int orb,bool d_mu_only=false) const;
  double local_weight_change(const segment &seg, int orb, bool antisegment) const;
  int order(int orbital) const{ return segments_[orbital].size(); }
  int n_orbitals() const{ return n_orbitals_; }
  bool zero_order_orbital_occupied(int orbital) const{return zero_order_orbital_occupied_[orbital]; }
  void set_zero_order_orbital_occupied(int orbital, bool value){zero_order_orbital_occupied_[orbital]=value; }
  double find_next_segment_end_distance(double time, int orbital);
  double find_next_segment_start_distance(double time, int orbital);
  void insert_segment(const segment &new_segment, int orbital);
  void insert_antisegment(const segment &new_segment, int orbital);
  void remove_segment(const segment &new_segment, int orbital);
  void remove_antisegment(const segment &new_segment, int orbital);
  segment get_segment(int k, int orbital) const;
  bool exists(double t) const{ return times_set_.find(t)==times_set_.end()?false:true;}
  bool has_overlap(const segment &seg,const int orb);
  void get_segment_densities(std::vector<std::vector<std::vector<double> > > &n_tauprime)const;
  void get_F_prefactor(std::vector<std::map<double,double> > &F_prefactor)const;
  void measure_density(std::vector<double> &densities, double sign) const;
  double segment_density(int i) const;
  double measure_nn(int i, int j) const;
  void measure_nnw(int i, std::vector<double> &nnw, double sign) const;
  void get_density_vectors(std::vector<std::vector<double> > &n_vector) const;
  double density(int i, double tau) const;
  double mu(int orbital) {return mu_[orbital];}
  double interaction_density_integral(std::set<segment>::const_iterator &it) const;
  void state_map_segment_insert(state_map &states, const segment &s, int state) const;
  void measure_sector_statistics(std::vector<double> &sector_statistics, double sign) const;
  friend std::ostream &operator<<(std::ostream &os, const local_configuration &local_conf);
  
  //debug functions
  void check_consistency() const;
  double full_weight() const;

private:
  //private member functions
  double segment_overlap(const segment &seg1, const segment &seg2) const;

  //private variables
  int crank_;
  interaction_matrix U_;
  chemical_potential mu_;
  ret_int_fun K_;
  
  double beta_;
  int n_orbitals_;
  bool use_retarded_interaction_;
  std::vector<std::set<segment> > segments_;
  std::vector<bool > zero_order_orbital_occupied_; //special case for perturbation order zero, where the orbital can either be occupied or empty. True means it is occupied, false is empty.
  std::set<double> times_set_; //this is a map making sure we don't have any times double, which would otherwise confuse the commutators.
} local_configuration;

std::ostream &operator<<(std::ostream &os, const local_configuration &local_conf);
#ifndef COLORS
#define COLORS
#define cblack "\033[22;30m"
#define cred "\033[22;31m"
#define cgreen "\033[22;32m"
#define cbrown "\033[22;33m"
#define cblue "\033[22;34m"
#define cmagenta "\033[22;35m"
#define ccyan "\033[22;36m"
#define cgray "\033[22;37m"
#define cdgray "\033[01;30m"
#define clred "\033[01;31m"
#define clgreen "\033[01;32m"
#define clyellow "\033[01;33m"
#define clblue "\033[01;34m"
#define clmagenta "\033[01;35m"
#define clcyan "\033[01;36m"
#define cwhite "\033[01;37m"
#endif

#endif
