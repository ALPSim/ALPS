/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Si_mulations
*
* ALPS Applications
*
* Copyright (C) 2006-2010 by Lode Pollet <lpollet@physics.harvard.edu>,
*                            Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>
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

/* $Id: Lowa.h 3520 2010-03-03 16:55:00Z tamama $ */


/*
 *
 * 1) Code modification      -- partly done
 * 2) Replacing raw pointers -- not done yet  -- use std::vector whenever you allocate something
 *
 */



#ifndef LOWA_H
#define LOWA_H


#include "include/cdefine.h"
#include "include/typedef.h"
#include "include/bipartite_lattice.h"     
#include "include/kink.h"
#include "band_structure_calculations/band_structure_calculations_BHM.h"

#include <alps/scheduler/montecarlo.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <sstream>
#include <fstream>
#include <alps/hdf5.hpp>
#include <vector>
#include <valarray>
#include <list>
#include <algorithm>
#include <numeric>




namespace alps {
namespace applications {


class Lowa : public alps::scheduler::MCRun
{
public:
  // copyright
  static void print_copyright(std::ostream&);


  // typedefs
  typedef alps::applications::lowa::site_type            site_type;                 
  typedef alps::applications::lowa::index_type           index_type; 
  typedef alps::applications::lowa::fock_basis_type      fock_basis_type;
  typedef alps::applications::lowa::dim_type             dim_type;
  typedef alps::applications::lowa::time_type            time_type;
  typedef alps::applications::lowa::real_coordinate_type real_coordinate_type;
  typedef alps::applications::lowa::parm_type            parm_type;
  typedef alps::applications::lowa::obs_type             obs_type;
  typedef alps::applications::lowa::vertex_type          vertex_type;
  typedef alps::applications::lowa::particle_number_type particle_number_type;

  typedef alps::applications::lowa::bipartite_lattice<site_type,real_coordinate_type> lattice_type;     
  typedef alps::applications::band_structure_calculations_BHM<site_type,parm_type> band_structure_type;

  typedef alps::applications::lowa::kink<site_type,time_type> kink_type;
  typedef std::list<kink_type> kinks_type;
  typedef kinks_type::iterator kinks_iterator_type;


  // get_functions
  inline site_type            dist  (site_type const s1, site_type const s2)        {  return _lattice.dist(s1,s2);  }
  inline real_coordinate_type distsq(index_type index, dim_type which_dimension)    {  return _lattice.distsq(index,which_dimension);  }
  inline real_coordinate_type mid(dim_type which_dimension)                         {  return _lattice.mid(which_dimension); }
  inline bool                 is_border_site(index_type index)                      {  return _lattice.is_border_site(index); }
  inline site_type            nb(index_type index, dim_type which_direction)        {  return _lattice.nb(index,which_direction); }

  inline parm_type tx_plus(site_type index)  const  {  return _band_structure.tx_plus(index); }
  inline parm_type ty_plus(site_type index)  const  {  return _band_structure.ty_plus(index); }
  inline parm_type tz_plus(site_type index)  const  {  return _band_structure.tz_plus(index); }
  inline parm_type tx_minus(site_type index) const  {  return _band_structure.tx_minus(index); }
  inline parm_type ty_minus(site_type index) const  {  return _band_structure.ty_minus(index); }
  inline parm_type tz_minus(site_type index) const  {  return _band_structure.tz_minus(index); }
  inline parm_type U(site_type index)        const  {  return _band_structure.U(index); }
  inline parm_type V(site_type index)        const  {  return _band_structure.V(index); }

  particle_number_type get_Npart()                  { if (new_measurement) update_N_and_E();  return (_nrpart); }
  particle_number_type get_Nmeaspart()              { if (new_measurement) update_N_and_E();  return (_meas_nrpart); }
  fock_basis_type      get_density(site_type index)     { if (new_measurement) update_N_and_E();  return (_state[index]);}
  fock_basis_type      get_measdensity(site_type index) 
  { 
    if (new_measurement) update_N_and_E();  
    if (is_doublon_treated_as_hole)  {  return (_state[index] % 2);  } 
    else                             {  return (_state[index]);  }
  }

  obs_type get_kinetic_energy()   {if (new_measurement) update_N_and_E(); return (_Ekin);}
  obs_type get_potential_energy() {if (new_measurement) update_N_and_E(); return (_Epot);}
  obs_type get_energy()           {if (new_measurement) update_N_and_E(); return (_Epot + _Ekin);} 
 
  inline double worm_insertion_acceptance_ratio() { return (no_of_accepted_worm_insertions / no_of_proposed_worm_insertions); }


private:

// ### PHYSICAL CONSTANTS

  double pi;
  double tol;
  double hbar;
  double amu;


// ### INTRINSIC SIMULATION VARIABLES

  // member functions
  void   dostep();     // ### included in <lowa.dostep.h>

  double work_done() const       { return (is_thermalized() ? (sweeps-thermalization_sweeps)/double(total_sweeps) : 0.); }
  bool   is_thermalized() const  { return (sweeps >= thermalization_sweeps); }
  bool   change_parameter(const std::string& name,const alps::StringValue& value);

  // member objects
  uint32_t label;
  uint64_t counter_MEASURE, counter_MEASURE_GREEN, counter_TEST;
  uint64_t Nmeasure, Nmeasure_green, Ntest, Nsave;
  uint64_t sweeps, sweeps_green, thermalization_sweeps, total_sweeps;
  double   MCstep_total, MCold;
  std::time_t times1, times2;


// ### PRINT and I/O (ASCII/hdf5)   // *** Tama : to include them in <lowa.io.h>

  // member functions
  std::string obtain_filename(std::string original_filename)     {  return (original_filename + boost::lexical_cast<std::string>(label));  }
  std::string obtain_filename_h5(std::string original_filename)  {  return (original_filename + boost::lexical_cast<std::string>(label) + ".h5");  }
  void export_lowa_simulation(std::ostream&);
  void import_lowa_simulation(std::istream&);
  void save(alps::ODump& odump)  {}
  void load(alps::IDump& idump)  {}
  void save_lowa();
  void load_lowa();

  std::ostream& print_params(std::ostream&) const;
  std::ostream& print_configuration(std::ostream&);

  std::ifstream inFile;
  std::ofstream outFile;

  std::string filename1, filename1_trial;
  std::string filename_N, filename_mN, filename_dns0, filename_mdns0, filename_proj_cymdns, filename_cs_cymdns;
  std::string filename_mdns;                                                                           // for measure_time_series_density == true
  std::string filename_dns, filename_dns_trial;                                                        // for measure_time_series_density == false
  std::string filename_mdnsmat, filename_mdnsmat_inf;                                                  // for measure_time_series_density_matrix == true
  std::string filename_dnsmat, filename_dnsmat_inf, filename_dnsmat_trial, filename_dnsmat_inf_trial;  // for measure_time_series_density_matrix == false

  bool measure_time_series_density;
  bool measure_time_series_density_matrix;


// WORLDLINE DESCRIPTION (ie. LATTICE (space), BAND PARAMETERS (space), LIST OF INTERACTIONS (space-time), DENSITY at a "plane")   

  dim_type dim;
  dim_type zcoord;

  site_type *Ls;
  site_type _Nxy, _N, _site0;

  parm_type *v0, *vc;
  parm_type *lambda, *waist;
  parm_type beta, _mu, _as, mass, _tof;
  parm_type *phase;
  site_type Ncan;

  fock_basis_type _statemax;   // maximum occupation number
  parm_type       _Eoffset;       // intrinsic variable

  lattice_type        _lattice;             
  band_structure_type _band_structure;
  kinks_type*         _kinks_on_lattice;
  parm_type*          _mu_eff;


// ### MEASUREMENT OBSERVABLES

  // member functions
  void calc_N_and_E();
  inline void update_N_and_E()   {  calc_N_and_E();  new_measurement = false;  }
  inline void reset_av_dns()     {  _Z_dns = 0.;  for (site_type index = 0; index < _N; ++index)  {  av_dns[index] = 0.;  }  }
  inline void reset_av_dnsmat()  {  _Z_dnsmat = 0.;  for (site_type index = 0; index < _N; ++index)  {  av_dnsmat[index] = 0.;  av_dnsmat_inf[index] = 0.;  }  }
  void statebinning();
  inline void update_obs()  {  if (new_measurement)  {  update_N_and_E(); }  statebinning();  }
  void update_off_diag_obs();
  void take_diagonal_measurements();
  void take_offdiagonal_measurements();

  // member objects
  bool new_measurement;
  bool is_doublon_treated_as_hole;
  vertex_type          _nrvertex;
  particle_number_type _nrpart, _meas_nrpart;
  fock_basis_type*     _state;             // a particular snap-shot of occupation no on lattice
  obs_type             _Ekin, _Epot;
  obs_type*            av_dns;
  obs_type             _Z_dns;
  obs_type*            av_dnsmat;
  obs_type*            av_dnsmat_inf;
  obs_type             _Z_dnsmat;

  parm_type  _binsize;
  site_type  _nrbin;
  site_type* _proj_binnr;
  site_type* _proj_binfreq;
  std::valarray<time_type> _proj_binstate;
  std::valarray<time_type> _cs_binstate;


// ### ABOUT WORM DESCRIPTION

  // member objects
  kinks_iterator_type   worm_it;  // iterator pointing at the next element the worm will encounter
  kinks_iterator_type*  site_it;  // iterator for site
  kinks_iterator_type*  dummy_kink_it; // iterator pointing at dummy elements useful for measuring diag prop

  bool is_worm_diagonal;
  fock_basis_type _state_enclosed_within_wormpair_at_creation;
  kink_type worm_head;                     // in Russian called ira
  kink_type worm_tail;                     // in Russian called masha
  bool is_worm_moving_forward;             // direction of the worm, to higher of lower imag times
  bool is_worm_rising;                     // convention : density is increased or decreased between tail and head


// ### ABOUT WORM INSERTION, PROPAGATION, DELETION PROCESSING

  // member functions
  bool worm_from_insertion_to_deletion();                                    //  ###  included in <lowa.insertion_deletion.h>
  void worm_propagation(time_type&);                                         //  ###  included in <lowa.propagation.h>

  dim_type select_index_via_heatbatch_algorithm();                                                            // ###  included in <lowa.propagation.h>
  void setup_local_transition_weight_upon_deletion(site_type const);                                          // ###  included in <lowa.propagation.h>
  void setup_local_transition_weight_upon_insertion(obs_type const , site_type const, time_type const);       // ###  included in <lowa.propagation.h>
  bool is_kink_encountered_forward(time_type&, obs_type const, time_type&);                                   // ###  included in <lowa.propagation.h>
  bool is_kink_encountered_backward(time_type&, obs_type const, time_type&);                                  // ###  included in <lowa.propagation.h>

  bool is_time0_cyclically_between_time1_and_time2(time_type const, time_type const, time_type const) const;  // ###  included in <lowa.insertion_deletion.h>
  inline void move_kinks_iterator_cyclically_forward(kinks_iterator_type&, site_type const);                  // ###  included in <lowa.insertion_deletion.h>
  inline void move_kinks_iterator_cyclically_backward(kinks_iterator_type&, site_type const);                 // ###  included in <lowa.insertion_deletion.h>
  void reset_assoc_upon_insert(site_type, kinks_iterator_type);                                               // ###  included in <lowa.insertion_deletion.h>
  void reset_assoc_upon_delete(site_type, kinks_iterator_type);                                               // ###  included in <lowa.insertion_deletion.h>

  // member objects
  fock_basis_type* nbs;
  obs_type* local_transition_weight;
  obs_type* cummulative_local_transition_weight;
  obs_type  total_local_transition_weight;

  double no_of_accepted_worm_insertions, no_of_proposed_worm_insertions;
  double MCstep;


// ### TEST 
  bool test_conf();

public:
  // destructors
  ~Lowa()
  {
    if (dim == 1) {
      delete Ls;
      delete vc;
      delete v0;
      delete phase;
      delete waist;
      delete lambda;
    }
    else {
      delete [] Ls;
      delete [] vc;
      delete [] v0;
      delete [] phase;
      delete [] waist;
      delete [] lambda;
    }
    delete [] _kinks_on_lattice;
    delete [] dummy_kink_it;
    delete [] site_it;
    delete [] _proj_binnr;
    delete [] nbs, local_transition_weight, cummulative_local_transition_weight;
    delete [] _mu_eff;
    delete [] _state, av_dns, av_dnsmat, av_dnsmat_inf;
  }
    
  // constructors
  Lowa(const alps::ProcessList&,const alps::Parameters&,int);

  // init functions
  void init();

};


} // ending namespace applications
} // ending namespace alps


using alps::applications::Lowa;
typedef alps::scheduler::SimpleMCFactory<Lowa> Factory;

#endif
