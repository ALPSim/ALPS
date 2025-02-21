/****************************************************************************
*
* ALPS Project Applications: Directed Worm Algorithm  
*
* Copyright (C) 2013 - 2016 
*               by  Matthias Troyer  <troyer@phys.ethz.ch> ,
*                   Lode Pollet      <pollet@phys.ethz.ch> ,
*                   Ping Nang Ma     <tamama@yotcopi.com> 
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef DWA_HPP
#define DWA_HPP


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

#include <alps/expression.h>
#include <alps/numeric/vector_functions.hpp>
#include <alps/numeric/vector_valarray_conversion.hpp>
#include <alps/ngs.hpp>

#include "worldlines.hpp"

#include "../qmc.h"


class directed_worm_algorithm : public QMCRun<>     
{
public:
    typedef boost::uint64_t                count_type;
    typedef worldlines::line               line;
    typedef worldlines::location_type      location_type;

    //directed_worm_algorithm(alps::hdf5::archive & ar);   // Comment by Tama: I couldn't proceed on without Luka's NGS scheduler being ready.
    directed_worm_algorithm(const alps::ProcessList& p, const alps::Parameters& parms_, int n);

    static void print_copyright  (std::ostream & out);
    void        print_simulation (std::ostream & out);

    // i/o
    void save(alps::hdf5::archive & ar) const;
    void load(alps::hdf5::archive & ar);

private:
    /// private member functions 

    // regarding simulation backbone (ESSENTIAL)
    void  dostep();
    void  measure() {}  // I don't measure every sweep, and I do a checkpoint rather then only measure.      
    void  docheckpoint();

    // regarding simulation interprocess (ESSENTIAL)
    bool   wormhead_propagates_till_collision_with_wormtail (unsigned short wormpair_state, std::pair<neighbor_iterator, neighbor_iterator> const & neighbors_);
    void   insert_jump_or_bounce        (double diagonal_energy_relative_, std::pair<neighbor_iterator, neighbor_iterator> const & neighbors_);
    void   delete_relink_jump_or_bounce (std::pair<neighbor_iterator, neighbor_iterator> const & neighbors_);

    // regarding simulation performance 
    bool   is_thermalized()              const  {  return (_sweep_counter > _thermalization_sweeps); }
    double work_done()                   const  {  return is_thermalized() ? static_cast<double>(_sweep_counter-_thermalization_sweeps)/double(_total_sweeps) : 0.;  }
    double probability_worm_insertion()  const  {  return (1. - static_cast<double>(_sweep_failure_counter)/_sweep_counter);         }
    double probability_bounce()          const  {  return (static_cast<double>(_propagation_failure_counter)/_propagation_counter);  }

    // regarding model
    void initialize_onsite_hamiltonian();   // called from initialize_hamiltonian within
    void initialize_onbond_hamiltonian();     // called from initialize_hamiltonian within
    void initialize_hamiltonian();
    void reset_diagonal_cache(unsigned site_, double time_);
    void print_hamiltonian(std::ostream & out);

//    std::vector<double> onsite_hamiltonian(unsigned int site_type);
    inline boost::multi_array<double,4>  
                      bond_hamiltonian(const bond_descriptor & bond);
    boost::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double> > 
                      bond_ladder_hamiltonians(const bond_descriptor & bond);

    inline double onsite_energy (unsigned int site_, unsigned short state_) const  { return onsite_matrix[inhomogeneous_site_type(site_)][state_]; }     
    inline double onbond_energy (bond_descriptor const & bond_, unsigned short state0_, unsigned short state1_) const { return onbond_matrix[bond_type(bond_)][state0_][state1_]; }
    double diagonal_energy_relative (unsigned int site_, unsigned short state_, bool forward_, bool creation_) const;
    double diagonal_energy_relative (unsigned int site_, unsigned short state_, std::vector<bond_descriptor> const & bonds_, std::vector<unsigned short> const & targetstates_, bool forward_, bool creation_) const;


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

    // CHANGING HERE...
    int componentIndex2Index(std::vector<int> componentIndex) const  
    {  
        if (dimension() == 1)
            return componentIndex[0];
        else
            if (dimension() == 2)
                return (componentIndex[0] * component_coordinate_x100_lookup[1].size() + componentIndex[1]);
            else
                return (componentIndex[0] * component_coordinate_x100_lookup[1].size() * component_coordinate_x100_lookup[2].size() + componentIndex[1] * component_coordinate_x100_lookup[2].size() + componentIndex[2]);
    }

    // regarding lattice
    std::vector<double>   position(unsigned int site)                   const  { return position_lookup[site]; }
    double                position(unsigned int site, unsigned int dim) const  { return position_lookup[site][dim]; } 

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

    count_type  _thermalization_sweeps;
    count_type  _total_sweeps;
    count_type  _skip;

    std::pair<std::time_t, std::time_t>  _simulation_timer;  

    // regarding model
    bool is_diagonal_onsite_;

    std::vector<unsigned short> state_minimum;         // arg: site site type
    std::vector<unsigned short> state_maximum;         // arg: site site type

    std::vector<std::vector<double> > onsite_matrix;               // onsite matrix      - arg0: sitetype, arg1: state  
    std::vector<std::vector<std::vector<double> > > onbond_matrix; // onbond matrix      - arg0: bondtype, arg1: state_i , arg2: state_j 
    std::vector<std::vector<double> > site_oneup_matrix;           // site 1-up matrix   - arg0: sitetype, arg1: state
    std::vector<std::vector<double> > site_onedown_matrix;         // site 1-down matrix - arg0: sitetype, arg1: state

    std::vector<double> bond_strength_matrix;                      // bond strength matrix - arg: bondtype

    // regarding lookups
    std::vector<std::vector<double> >  position_lookup;
    std::vector<std::vector<int> >     component_coordinate_x100_lookup; 
    std::vector<int>                   site_lookup; 
    std::vector<double>                phase_lookup;

    // regarding lattice
    bool is_periodic_;
    std::vector<std::vector<double> > lattice_vector_;

    // regarding experiment
    bool finite_tof;
    bool finite_waist;

    // regarding worldline
    worldlines wl;
    wormpair   worm;

    // regarding restrictions
    bool restrict_measurements_total_particle_number_;

    int measure_total_particle_number_at_;

    // regarding measurements
    bool measure_;
    bool measure_simulation_speed_;
    bool measure_number2_;
    bool measure_energy2_;
    bool measure_density2_;
    bool measure_energy_density2_;
    bool measure_winding_number_;
    bool measure_local_num_kinks_;
    bool measure_local_density_;
    bool measure_local_density2_;
    bool measure_green_function_;

    std::vector<double>  _num_kinks_cache;
    std::vector<double>  _states_cache;
    std::vector<double>  _states2_cache;

    // regarding on-fly measurements
    double  green0;
    double  green1;
    double  nk0;
    double  nk0_tof;

    std::vector<double>  green;
    std::vector<double>  green_tof;

    // regarding caches (for optimization purpose)
#ifdef HEATBATH_ALGORITHM
    std::vector<location_type>  _neighborlocations_cache;
    std::vector<double>         _cummulative_weights_cache; 
#endif
    std::vector<bond_descriptor>  _neighborbonds_cache;
    std::vector<unsigned>         _neighborsites_cache;
    std::vector<unsigned short>   _neighborstates_cache;
    std::vector<double>           _neighbortimes_cache;
    std::vector<double>           _neighbortaus_cache;
};
#endif
