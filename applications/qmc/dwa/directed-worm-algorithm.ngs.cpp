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

#include "directed-worm-algorithm.ngs.hpp"


// ==================================================
// directed_worm_algorithm member functions
// ==================================================

void
  directed_worm_algorithm
    ::print_copyright (std::ostream & out)
    {
      out << "\n\n"
          << "\n/*****************************************************************************"
          << "\n*"
          << "\n* ALPS Project Applications: Directed Worm Algorithm (using NGS scheduler)"
          << "\n*"
          << "\n* Copyright (C) 2012 by Lode Pollet      <pollet@phys.ethz.ch>,"
          << "\n*                       Ping Nang Ma     <pingnang@phys.ethz.ch>,"
          << "\n*                       Matthias Troyer  <troyer@phys.ethz.ch>"
          << "\n*"
          << "\n* This software is part of the ALPS libraries, published under the ALPS"
          << "\n* Library License; you can use, redistribute it and/or modify it under"
          << "\n* the terms of the license, either version 1 or (at your option) any later"
          << "\n* version."
          << "\n*"
          << "\n* You should have received a copy of the ALPS Library License along with"
          << "\n* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also"
          << "\n* available from http://alps.comp-phys.org/."
          << "\n*"
          << "\n* THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR"
          << "\n* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,"
          << "\n* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT"
          << "\n* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE"
          << "\n* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,"
          << "\n* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER"
          << "\n* DEALINGS IN THE SOFTWARE."
          << "\n*"
          << "\n*****************************************************************************/"
          << "\n\n";
    }

void
  directed_worm_algorithm
    ::print_simulation (std::ostream & out)
    {
      using alps::numeric::operator<<;
      out << "Parameters\n"
          << "==========\n\n"
          << parameters
          << "\n\n";
      out << "Simulation\n"
          << "==========\n\n"
          << "\t Total sweeps                      : " << _total_sweeps                  << "\n"
          << "\t Sweep per measurement             : " << _sweep_per_measurement         << "\n"
          << "\t Sweep counter                     : " << _sweep_counter                 << "\n"
          << "\t Sweep failure counter             : " << _sweep_failure_counter         << "\n"
          << "\t Propagation counter               : " << _propagation_counter           << "\n"
          << "\t Propagation failure counter       : " << _propagation_failure_counter   << "\n"
          << "\n\n";
      out << "Lattice\n"
          << "=======\n\n"
          << "Name              : " << parameters["LATTICE"]                    << "\n"
          << "Dimension         : " << dimension()                              << "\n"
          << "Total Sites/Bonds : " << num_sites() << "\t/\t" << num_bonds()    << "\n"
          << "Periodic          : " << is_periodic_ << "\n"
#ifdef DEBUGMODE
          << "Graph             : " << "\n"
          << graph()
#endif
          << "\n\n";
      out << "Model\n"
          << "=====\n\n"
          << "Name         : " << model().name()   << "\n"
          << "Charge model : " << is_charge_model_ << "\n"
          << "Spin model   : " << is_spin_model_   << "\n"
          << "Site Basis   : " << model().basis().site_basis().name() << "\n"
          << "states       : " << diagonal_matrix_element.begin()->second[0] << "\n"
          << "\n\n";
      out << "Measurements\n"
          << "============\n"
          << "Measure                     : " << !measure_only_simulation_speed_ << "\n"
          << "Measure Winding Number      : " << measure_winding_number_         << "\n"
          << "Measure Local Energy        : " << measure_local_energy_           << "\n"
          << "Measure Local Density       : " << measure_local_density_          << "\n" 
          << "Measure Local Density^2     : " << measure_local_density2_         << "\n"
          << "Measure Green's Function    : " << measure_green_function_         << "\n"
          << "Measure Momentum Density    : " << measure_momentum_density_       << "\n"
          << "Measure TOF Image           : " << measure_tof_image_              << "\n"
          << "\n\n"; 
    }


directed_worm_algorithm
  ::directed_worm_algorithm(alps::hdf5::archive & ar)
    : qmcbase<>(ar) 
    // regarding MC simulation
    , _sweep_counter               (0)
    , _sweep_failure_counter       (0)
    , _propagation_counter         (0)
    , _propagation_failure_counter (0)
    , _total_sweeps                (this->parameters["TOTAL_SWEEPS"] | 10000000)
    , _sweep_per_measurement       (this->parameters["SKIP"] | 1)
    // regarding lattice
    , is_periodic_ (std::find((this->lattice().boundary()).begin(), (this->lattice().boundary()).end(), "open") == (this->lattice().boundary()).end())
    , num_component_momenta_ (1 + static_cast<int>(this->parameters["MOMENTUM_EXTENT"] | 10))
    // regarding worldline
    , wl (num_sites())
    // regarding experiment
    , finite_tof (is_periodic_ ? false : (0. != static_cast<double>(this->parameters["time_of_flight"] | 0.)))
    // regarding measurements
    , measure_only_simulation_speed_ (!(this->parameters["MEASURE"] | true))
    , measure_winding_number_        (is_periodic_ ? static_cast<bool>(this->parameters["MEASURE[Winding Number]"] | false) : false)
    , measure_local_density_         (this->parameters["MEASURE[Local Density]"] | false)
    , measure_local_density2_        (this->parameters["MEASURE[Local Density^2]"] | false)
    , measure_local_energy_          (this->parameters["MEASURE[Local Energy]"] | false)
    , measure_green_function_        (is_periodic_ ? static_cast<bool>(this->parameters["MEASURE[Green Function]"] | false) : false)
    , measure_momentum_density_      (this->parameters["MEASURE[Momentum Density]"] | false)
    , measure_tof_image_             (finite_tof ? static_cast<bool>(this->parameters["MEASURE[TOF Image]"] | false) : false)
  {
    // lattice enhancement
    using alps::numeric::operator*;
    for (int i=0; i<dimension(); ++i)
      lattice_vector_.push_back(*(basis_vectors().first + i) * lattice().extent()[i]);

    // other initialization
    std::cout << "\nInitialization stage 1 (Site states)\t... starting\n"; 
    initialize_site_states(); 
    std::cout << "\t\t... finished.\n"; 
 
    std::cout << "\nInitialization stage 2 (Hamiltonian)\t... starting\n"; 
    initialize_hamiltonian(); 
    std::cout << "\t\t... finished.\n"; 
	 
    std::cout << "\nInitialization stage 3 (Lookups)\t... starting\n"; 
    initialize_lookups(); 
    std::cout << "\t\t... finished.\n"; 
	 
    std::cout << "\nInitialization stage 4 (Measurements)\t... starting\n"; 
    initialize_measurements(); 
    std::cout << "\t\t... finished.\n"; 

    // regarding caches
#ifdef HEATBATH_ALGORITHM
    unsigned int maximum_num_neighbors=0;
    for (site_iterator it=sites().first; it!=sites().second; ++it)
      if (num_neighbors(*it) > maximum_num_neighbors)
        maximum_num_neighbors = num_neighbors(*it);
    
    _neighborlocations_cache.reserve(maximum_num_neighbors);
    _cummulative_weights_cache.reserve(maximum_num_neighbors+1);
#endif

    // print to screen
    print_simulation (std::cout);
#ifdef DEBUGMODE
    print_hamiltonian(std::cout);
    print_lookups    (std::cout);
    print_worldline  (std::cout);
#endif

    // start timer if needed...
    if (measure_only_simulation_speed_)
      std::time(&_simulation_timer.first);
  }


void
  directed_worm_algorithm
    ::save(alps::hdf5::archive & ar) const
    {
      // Worldlines
      ar << alps::make_pvp("/simulation/worldlines/num_sites", num_sites()); 
      ar << alps::make_pvp("/simulation/worldlines", wl);

      // Parameters
      ar << alps::make_pvp("/parameters", parameters);

      // Simulation results
      ar << alps::make_pvp("/simulation/results", measurements);

      // Other useful quantities
      ar << alps::make_pvp("/simulation/positions", position_lookup);
      if (measure_green_function_)
        ar << alps::make_pvp("/simulation/green_coordinates", lattice().distance_labels());
      if (measure_momentum_density_ || measure_tof_image_)
        ar << alps::make_pvp("/simulation/momenta", component_momenta_lookup);
    }

void 
  directed_worm_algorithm
    ::load(alps::hdf5::archive & ar)  
    {  
      // Worldlines
      ar >> alps::make_pvp("/simulation/worldlines", wl); 
      std::cout << "\n\nWorldlines configuration loaded from input file.\n\n";

      // Measurements
      ar >> alps::make_pvp("/simulation/results", measurements);
      std::cout << "\n\nPrevious measurements loaded from input file.\n\n";
    }

void
  directed_worm_algorithm
    ::initialize_hamiltonian()
    {
      // setup state_minimum and state_maximum
      std::cout << "\t\t\ti. State minimum/maximum \t... starting\n";
      for (int i=0; i<maximum_sitetype+1; ++i)
      {
        state_minimum.push_back(boost::math::iround(diagonal_matrix_element.begin()->second[i].front()));
        state_maximum.push_back(boost::math::iround(diagonal_matrix_element.begin()->second[i].back()));
      }
      std::cout << "\t\t\t\t... done.\n"; 

      // iterate all sites and build onsite matrix 
      std::cout << "\t\t\ti. Onsite matrix \t... starting\n";
      for (site_iterator it = sites().first; it != sites().second; ++it) {
        unsigned int this_site_type = inhomogeneous_site_type(*it);
        if (this_site_type >= onsite_matrix.size())
        { 
          using boost::numeric::operators::operator*;
          onsite_matrix.resize(this_site_type+1);
          onsite_matrix[this_site_type] = onsite_hamiltonian(this_site_type) * beta;
        }
      }
      std::cout << "\t\t\t\t... done.\n";

      // iterate all bonds and build site 1-up matrix, site 1-down matrix and bond strength matrix
      std::cout << "\t\t\tiii. Bond strength matrix \t... starting\n";
      std::vector<double> bond_strength_matrix_temporary;
      for (bond_iterator it = bonds().first; it != bonds().second; ++it) {
        unsigned int this_bond_type       = inhomogeneous_bond_type(*it);
        unsigned int this_sourcesite_type = site_site_type[source(*it)];
        unsigned int this_targetsite_type = site_site_type[target(*it)];

        if (std::max(this_sourcesite_type, this_targetsite_type) >= site_oneup_matrix.size())
        {
          boost::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double> >
            bond_ladder_matrices = bond_ladder_hamiltonians(*it); 
      
          if (this_sourcesite_type >= site_oneup_matrix.size())
          {
            site_oneup_matrix.resize(this_sourcesite_type+1);
            site_oneup_matrix[this_sourcesite_type]   = bond_ladder_matrices.get<0>();
          }

          if (this_sourcesite_type >= site_onedown_matrix.size())
          {
            site_onedown_matrix.resize(this_sourcesite_type+1);
            site_onedown_matrix[this_sourcesite_type] = bond_ladder_matrices.get<1>();
          }

          if (this_targetsite_type >= site_oneup_matrix.size())
          {
            site_oneup_matrix.resize(this_targetsite_type+1);
            site_oneup_matrix[this_targetsite_type]   = bond_ladder_matrices.get<2>();
          }

          if (this_targetsite_type >= site_onedown_matrix.size())
          {
            site_onedown_matrix.resize(this_targetsite_type+1);
            site_onedown_matrix[this_targetsite_type] = bond_ladder_matrices.get<3>();
          }
        }

        if (this_bond_type >= bond_strength_matrix_temporary.size())
        {
          bond_strength_matrix_temporary.resize(this_bond_type+1);
          bond_strength_matrix_temporary[this_bond_type] 
            = std::abs(bond_hamiltonian(*it)[0][1][1][0]) / (site_oneup_matrix[this_sourcesite_type][0] * site_onedown_matrix[this_targetsite_type][1]) * beta;
        }

        bond_strength_matrix.push_back(bond_strength_matrix_temporary[this_bond_type]); 
      }

      std::cout << "\t\t\t\t... done.\n"; 
    }

void
  directed_worm_algorithm
    ::print_hamiltonian(std::ostream & out)
    {
      using alps::numeric::operator<<;

      out << "Hamiltonian:\n"
          << "============\n\n"
          << model() << "\n"
          << "============\n\n";

      for (site_iterator it=sites().first; it!=sites().second; ++it)
      {
        out << "Site/Site Type/Position : " << *it << "\t/\t" << site_site_type[*it] << "\t/\t" << position(*it) << "\n\n"
            << "unsigned shorts     : " << (diagonal_matrix_element.begin()->second)[site_site_type[*it]]  << "\n\n"
            << "Onsite H   : " << onsite_matrix[inhomogeneous_site_type(*it)] << "\n\n"
            << "+ operator : " << site_oneup_matrix  [site_site_type[*it]] << "\n\n"
            << "- operator : " << site_onedown_matrix[site_site_type[*it]] << "\n\n"
            << "---------------------------------------------------\n\n";
      }

      for (bond_iterator it=bonds().first; it!=bonds().second; ++it)
      {
        out << "Bond/Bond Index/Bond Type : " << *it << "\t/\t" << index(*it) << "\t/\t" << bond_type(*it) << "\n\n"
            << "Source Site/Site Type/Position : " << source(*it) << "\t/\t" << site_site_type[source(*it)] << "\t/\t" << position(source(*it)) << "\n\n"
            << "Target Site/Site Type/Position : " << target(*it) << "\t/\t" << site_site_type[target(*it)] << "\t/\t" << position(target(*it)) << "\n\n"
            << "Bond strength : " << bond_strength_matrix[index(*it)] << "\n\n"
            << "---------------------------------------------------\n\n";
      }

      out << "\n\n";
    }

void
  directed_worm_algorithm
    ::initialize_lookups()
    {
      // position lookup table
      position_lookup.reserve(num_sites());
      for (unsigned int site=0; site<num_sites(); ++site)
      {
        std::vector<double> _position = coordinate(site);
        for (int i=0; i<dimension(); ++i)
          _position[i] -= static_cast<double>(lattice().extent()[i]-1)/2.;
        position_lookup.push_back(_position);
      }

      // component momenta lookup table 
      component_momenta_lookup.reserve(num_component_momenta_);
      for (unsigned int idx=0; idx < num_component_momenta(); ++idx)
        component_momenta_lookup.push_back(idx*M_PI/lattice().extent()[0]/(num_component_momenta()-1));

      // phase and phase lookup table
      if (finite_tof)
      {
        std::vector<double> inverse_tof = std::vector<double>(dimension(), 1./static_cast<double>(parameters["time_of_flight"]));
        phase_lookup.reserve(num_sites());
        for (unsigned int site=0; site<num_sites(); ++site) 
          phase_lookup.push_back(std::inner_product(inverse_tof.begin(), inverse_tof.end(), (alps::numeric::sq(position(site))).begin(), 0.));
      }
    }

void
  directed_worm_algorithm
    ::print_lookups(std::ostream & out) const
    {
      using alps::numeric::operator<<;
   
      out << "Lookup Tables:\n"
          << "=============\n\n";

      // position lookup table
      out << "Position Lookup Table:\n"
          << "----------------------\n\n";
      for (unsigned int site=0; site<num_sites(); ++site)
        out << "Site : " << site << " , Site position : " << position(site) << "\n";
      out << "\n";

      // component lookup table
      out << "Component Momenta Lookup Table:\n"
          << "-------------------------------\n\n";
      for (unsigned int idx=0; idx<num_component_momenta(); ++idx)
        out << "Index : " << idx << " , Component Momentum : " << component_momentum(idx) << "\n";
      out << "\n";

      // phase lookup table
      out << "Phase Lookup Table:\n"
          << "----------------------\n\n";
      for (unsigned int site=0; site<num_sites(); ++site)
        out << "Site : " << site << " , Phase : " << phase(site) << "\n";
      out << "\n";

      out << "\n\n";

      std::cin.get();
    }

void
  directed_worm_algorithm
    ::print_worldline(std::ostream & out) const
    {
      out << "Worldline:\n"
          << "==========\n\n"
          << wl << "\n"
          << "==========\n\n";
    }

void
  directed_worm_algorithm
    ::initialize_measurements()
    {
      // regarding measurements
      measurements 
        << alps::ngs::RealObservable   ("Total Particle Number")
        << alps::ngs::RealObservable   ("Total Particle Number^2")
        << alps::ngs::RealObservable   ("Density")
        << alps::ngs::RealObservable   ("Density^2")
        << alps::ngs::RealObservable   ("Hopping Energy")
        << alps::ngs::RealObservable   ("Onsite Energy")
        << alps::ngs::RealObservable   ("Energy")
        << alps::ngs::RealObservable   ("Energy^2")
        << alps::ngs::RealObservable   ("Hopping Energy Density")
        << alps::ngs::RealObservable   ("Onsite Energy Density")
        << alps::ngs::RealObservable   ("Energy Density")
        << alps::ngs::RealObservable   ("Energy Density^2")
        << alps::ngs::RealObservable   ("Total Vertex Number")
        << alps::ngs::RealObservable   ("Total Vertex Number^2")
      ;

      if (measure_winding_number_)
        measurements << alps::ngs::RealVectorObservable ("Winding Number^2");

      if (measure_local_density_)
      {
        _states_cache.resize(wl.num_sites(),0.);
        measurements << alps::ngs::SimpleRealVectorObservable ("Local Density");
      }

      if (measure_local_density2_)
      {
        _states2_cache.resize(wl.num_sites(),0.);
        measurements << alps::ngs::SimpleRealVectorObservable ("Local Density^2");
      }

      // regarding on-fly measurements
      reinitialize_on_fly_measurements();

      measurements
        << alps::ngs::RealObservable  ("Green's Function Onsite")
        << alps::ngs::RealObservable  ("Green's Function Neighbors")
        << alps::ngs::RealObservable  ("Zero Momentum Density")
        ;

      if (measure_green_function_)
      {
        green.resize(lattice().num_distances(),0.);
        measurements << alps::ngs::SimpleRealVectorObservable ("Green's Function");
      }

      if (measure_momentum_density_)
      {
        momentum_density.resize(num_component_momenta(),0.);
        measurements << alps::ngs::SimpleRealVectorObservable ("Momentum Density");
      }

      if (finite_tof)
      {
        measurements << alps::ngs::RealObservable  ("Zero TOF Image");

        if (measure_green_function_)
        {
          green_tof.resize(lattice().num_distances(),0.);
          measurements << alps::ngs::SimpleRealVectorObservable ("TOF Green's Function");
        }

        if (measure_tof_image_)
        {
          momentum_density_tof.resize(num_component_momenta(),0.);
          measurements << alps::ngs::SimpleRealVectorObservable ("TOF Image");
        }
      }
    }

void
  directed_worm_algorithm
    ::reinitialize_on_fly_measurements()
    {
      green_onsite          = 0.;
      green_neighbors       = 0.;
      zero_momentum_density = 0.;

      if (measure_green_function_)
        std::fill(green.begin(),green.end(),0.);
      if (measure_momentum_density_)
        std::fill(momentum_density.begin(),momentum_density.end(),0.);

      if (finite_tof)
      {
        zero_momentum_density_tof = 0.;

        if (measure_green_function_)
          std::fill(green_tof.begin(),green_tof.end(),0.);
        if (measure_tof_image_)
          std::fill(momentum_density_tof.begin(),momentum_density_tof.end(),0.);
      }
    }

std::vector<double> 
  directed_worm_algorithm
    ::onsite_hamiltonian( unsigned int this_site_type)
    {
      std::vector<alps::SiteTermDescriptor> all_site_terms = model().site_terms();
      std::map<std::string, alps::SiteTermDescriptor> this_site_terms;
      for (std::vector<alps::SiteTermDescriptor>::iterator it=all_site_terms.begin(); it!=all_site_terms.end(); ++it)
        if (it->match_type(site_site_type[this_site_type]))
          this_site_terms.insert(std::make_pair(it->name(), *it));

      std::vector<double> this_onsite_hamiltonian;
      for (std::map<std::string, alps::SiteTermDescriptor>::iterator it=this_site_terms.begin(); it!=this_site_terms.end(); ++it)
      {
        alps::Parameters this_parms;
        if(inhomogeneous_sites()) {
          alps::throw_if_xyz_defined(this_copy_of_params_is_reserved_for_the_old_scheduler_library_which_is_to_be_depreciated,this_site_type);      // check whether x, y, or z is set
          this_parms << coordinate_as_parameter(this_site_type); // set x, y and z
        }
        boost::multi_array<double,2> this_onsite_matrix =
          alps::get_matrix( double(), it->second
                          , model().basis().site_basis(this_site_type), this_parms);
        std::vector<double> this_onsite_matrix_diagonal;
        for (unsigned short i=0; i < this_onsite_matrix.shape()[0]; ++i)
          this_onsite_matrix_diagonal.push_back(this_onsite_matrix[i][i]);
        using boost::numeric::operators::operator+;
        if (this_onsite_hamiltonian.empty())
          this_onsite_hamiltonian = this_onsite_matrix_diagonal;
        else
          this_onsite_hamiltonian = this_onsite_hamiltonian + this_onsite_matrix_diagonal;
      }
      return this_onsite_hamiltonian;
    }   

inline boost::multi_array<double,4>
  directed_worm_algorithm
    ::bond_hamiltonian(const bond_descriptor & bond)
    {
      alps::Parameters this_parms;
      if(inhomogeneous_bonds()) {
        alps::throw_if_xyz_defined(this_copy_of_params_is_reserved_for_the_old_scheduler_library_which_is_to_be_depreciated,bond);      // check whether x, y, or z is set
        this_parms << coordinate_as_parameter(bond); // set x, y and z
      }

      return alps::get_matrix( double(), model().bond_term(bond_type(bond))
                             , model().basis().site_basis(site_type(source(bond)))
                             , model().basis().site_basis(site_type(target(bond)))
                             , this_parms);
    }

boost::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double> >
  directed_worm_algorithm
    ::bond_ladder_hamiltonians(const bond_descriptor & bond)
    {
      std::vector<double> this_sourcesite_oneup_hamiltonian;
      std::vector<double> this_sourcesite_onedown_hamiltonian;
      std::vector<double> this_targetsite_oneup_hamiltonian;
      std::vector<double> this_targetsite_onedown_hamiltonian;

      unsigned int source_num_states = model().basis().site_basis(source(bond)).num_states();
      unsigned int target_num_states = model().basis().site_basis(target(bond)).num_states();

      this_sourcesite_oneup_hamiltonian  .resize(source_num_states);
      this_sourcesite_onedown_hamiltonian.resize(source_num_states);
      this_targetsite_oneup_hamiltonian  .resize(target_num_states);
      this_targetsite_onedown_hamiltonian.resize(target_num_states);

      std::set<std::string> site_operators = model().bond_term(bond_type(bond)).operator_names(this_copy_of_params_is_reserved_for_the_old_scheduler_library_which_is_to_be_depreciated);

      for (std::set<std::string>::iterator site_operator_it = site_operators.begin(); site_operator_it != site_operators.end(); ++site_operator_it)
      {
        alps::Parameters this_parms;
        if(inhomogeneous_bonds()) {
          alps::throw_if_xyz_defined(this_copy_of_params_is_reserved_for_the_old_scheduler_library_which_is_to_be_depreciated,bond);      // check whether x, y, or z is set
          this_parms << coordinate_as_parameter(bond); // set x, y and z
        }
        boost::multi_array<double,2> this_sourcesite_operator_matrix 
          = alps::get_matrix(double(), alps::SiteOperator(*site_operator_it), model().basis().site_basis(source(bond)), this_parms);
        boost::multi_array<double,2> this_targetsite_operator_matrix
          = alps::get_matrix(double(), alps::SiteOperator(*site_operator_it), model().basis().site_basis(target(bond)), this_parms);

        for (unsigned int i=0; i < source_num_states-1; ++i)
        {
          this_sourcesite_oneup_hamiltonian[i]     += this_sourcesite_operator_matrix[i][i+1];
          this_sourcesite_onedown_hamiltonian[i+1] += this_sourcesite_operator_matrix[i+1][i]; 
        }
        for (unsigned int i=0; i < target_num_states-1; ++i)
        {
          this_targetsite_oneup_hamiltonian[i]     += this_targetsite_operator_matrix[i][i+1];
          this_targetsite_onedown_hamiltonian[i+1] += this_targetsite_operator_matrix[i+1][i];
        }
      }

      return boost::tuples::make_tuple( this_sourcesite_oneup_hamiltonian, this_sourcesite_onedown_hamiltonian
                                      , this_targetsite_oneup_hamiltonian, this_targetsite_onedown_hamiltonian
                                      );
    }

double
  directed_worm_algorithm
    ::onsite_energy_relative (unsigned int site_, unsigned short state_, bool forward_, bool creation_) const 
    {
      const double energy_offset = 0.1;  

      const double _onsite_energy        = onsite_energy(site_, state_);
      const double _onsite_energy_before = onsite_energy(site_, creation_ ? state_-1 : state_+1);

      return (forward_ ? _onsite_energy_before  - std::min(_onsite_energy_before, _onsite_energy) + energy_offset
                       : _onsite_energy         - std::min(_onsite_energy_before, _onsite_energy) + energy_offset); 
    }

double
  directed_worm_algorithm
    ::hopping_energy (bond_descriptor const & bond_, unsigned short targetstate_, bool increasing_) const 
    {  
      return (increasing_ ? (bond_strength_matrix[index(bond_)] * site_oneup_matrix[site_site_type[target(bond_)]][targetstate_] * site_onedown_matrix[site_site_type[target(bond_)]][targetstate_+1]) : (bond_strength_matrix[index(bond_)] * site_onedown_matrix[site_site_type[target(bond_)]][targetstate_] * site_oneup_matrix[site_site_type[target(bond_)]][targetstate_-1]));
    }

void
  directed_worm_algorithm
    ::update()
    {
      ++_sweep_counter;
#ifdef DEBUGMODE
//      if (_sweep_counter >= DEBUGMODE_START_COUNTER && _sweep_counter <= DEBUGMODE_END_COUNTER) 
      {
        std::cout << "\n\nSweep : " << _sweep_counter << "\t(Propagation sweep: " << _propagation_counter << ")"
                  << "\n============"
                  << "\n\n";
        //std::cin.get();
      }
#endif

      // Step 1A: Generate a random wormpair
      const unsigned int _site     = random()*num_sites();
      const double       _time     = random();
      const bool         _forward  = (random()<0.5);
      const bool         _creation = (random()<0.5);

      location_type      _location = wl.location(_site,_time);
#ifdef DEBUGMODE
//      if (_sweep_counter >= DEBUGMODE_START_COUNTER && _sweep_counter <= DEBUGMODE_END_COUNTER) 
      {
        std::cout << "\n\nA random wormpair is generated:"
                  << "\n----------------------------------"
                  << "\nSite       : " << _site
                  << "\nTime       : " << _time
                  << "\nForward    : " << _forward
                  << "\nCreation   : " << _creation
                  << "\nIncreasing : " << wormpair::increasing(_forward,_creation)
                  << "\nLocation is kink-unoccupied : " << wl.location_is_kink_unoccupied(_location,_time)
                  << "\n\n";
        //std::cin.get();
      }
#endif
      if (!wl.location_is_kink_unoccupied(_location,_time))
      { // redo this sweep if location is occupied by kink
        --_sweep_counter;
        return;  
      } 

      const bool checkpoint_sweep = (_sweep_counter % _sweep_per_measurement == 0);
      const bool _measure = (!measure_only_simulation_speed_);

      // Step 1B: Check state to see if we can insert wormpair 
      const unsigned short _state = wl.state_before(_location);

      if (  (_state == state_minimum[site_site_type[_site]] && !wormpair::increasing(_forward, _creation)) 
         || (_state == state_maximum[site_site_type[_site]] && wormpair::increasing(_forward, _creation))  
         )
      {
        ++_sweep_failure_counter;      // worm insertion is unsuccessful         
        if (_measure && _state != 0) 
        {
          green_onsite          += _state;
          zero_momentum_density += _state;

          if (measure_green_function_)
            green[0] += _state;
          if (measure_momentum_density_)
            for (unsigned int idx=0; idx<num_component_momenta();++idx)
              momentum_density[idx] += _state;

          if (finite_tof)
          {
            zero_momentum_density_tof += _state;

            if (measure_green_function_)
              green_tof[0] += _state;
            if (measure_tof_image_)
              for (unsigned int idx=0; idx<num_component_momenta();++idx)
                momentum_density_tof[idx] += _state;
          }
        }
#ifdef DEBUGMODE
//        if (_sweep_counter >= DEBUGMODE_START_COUNTER && _sweep_counter <= DEBUGMODE_END_COUNTER) 
        {
          std::cout << "\n\n"
                    << "\nWormpair creation fails!" 
                    << "\n========================"
                    << wl
                    << "\nEnded (Wormpair creation fails!)"
                    << "\n--------------------------------"
                    << "\n\n";
          //std::cin.get();
        }
#endif
        if (checkpoint_sweep) 
          docheckpoint(_measure);
        return;
      }

      // Step 1C: Insert wormpair
      worm = wormpair(_location, kink(_site,_time,_state), _forward, _creation);
#ifdef DEBUGMODE
//      if (_sweep_counter >= DEBUGMODE_START_COUNTER && _sweep_counter <= DEBUGMODE_END_COUNTER) 
      {
        std::cout << "\n\n"
                  << "\nWormpair creation!" 
                  << "\n=================="
                  << wl
                  << "\n\n"
                  << worm
                  << "\nEnded (Wormpair creation!)"
                  << "\n--------------------------------"
                  << "\n\n";
        //std::cin.get();
      }
#endif

      // Step 2: Wormhead propagates till collision with wormtail
      while(wormhead_propagates_till_collision_with_wormtail(worm.wormpair_state(), neighbors(worm.wormtail_site()), _measure));      
      const unsigned short termination_state = (worm.forward() ? worm.state_before() : worm.state());

      if (_measure && termination_state != 0)
      {
        green_onsite          += termination_state;
        zero_momentum_density += termination_state;

        if (measure_green_function_)
          green[0] += termination_state;
        if (measure_momentum_density_)
          for (unsigned int idx=0; idx<num_component_momenta();++idx)
            momentum_density[idx] += termination_state;
            
        if (finite_tof)
        {
          zero_momentum_density_tof += termination_state;

          if (measure_green_function_)
            green_tof[0] += termination_state;
          if (measure_tof_image_)
            for (unsigned int idx=0; idx<num_component_momenta();++idx)
              momentum_density_tof[idx] += termination_state;
        }
      }

      worm.wormhead_annihilates_wormtail();
#ifdef DEBUGMODE
//      if (_sweep_counter >= DEBUGMODE_START_COUNTER && _sweep_counter <= DEBUGMODE_END_COUNTER) 
      {
      std::cout << "\n\n"
                << "\nWormpair annihilation!"
                << "\n======================"
                << wl
                << "\nEnded (Wormpair annihilation!)"
                << "\n------------------------------------------------------"
                << "\n\n";
      //std::cin.get();
      }
#endif

      if (checkpoint_sweep) 
        docheckpoint(_measure);
      return;
    }

void
  directed_worm_algorithm
    ::docheckpoint(bool const & measure_)
    {
      std::cout << "Checkpoint Sweep " << _sweep_counter
                << " ... Probability : " << probability_worm_insertion() << " / " << probability_bounce() 
                ;

      // measurements
      if (measure_)
      {
        perform_diagonal_measurements();
        reinitialize_on_fly_measurements();
      }

      // simulation timer
      if (measure_only_simulation_speed_)
      {
        std::time(&_simulation_timer.second);
        std::cout << " ... speed = " << std::difftime(_simulation_timer.second, _simulation_timer.first)/_sweep_per_measurement; 
        std::time(&_simulation_timer.first);
      }      

      std::cout << "\n";
    }

void
  directed_worm_algorithm
    ::perform_diagonal_measurements()
    {
      // regarding measurements
      using boost::numeric::operators::operator*;
      using boost::numeric::operators::operator/;

      int total_density  = 0;
      int total_density2 = 0;
      double total_onsite_energy  = 0.;
      double total_hopping_energy = 0.;
      double total_energy  = 0.;
      double total_energy2 = 0.;
      double total_vertex_number  = 0.;
      double total_vertex_number2 = 0.;

      for (unsigned int site=0; site < num_sites(); ++site)
      {
        total_density  += wl.site_state(site);
        total_density2 += wl.site_state(site) * wl.site_state(site);

        if (measure_local_density_)
          _states_cache[site] = wl.site_state(site);
        if (measure_local_density2_)
          _states2_cache[site] = wl.site_state(site) * wl.site_state(site);

        const double this_onsite_energy  = onsite_energy(site, wl.site_state(site));
        const double this_hopping_energy = -static_cast<double>(wl.num_kinks(site)-1)/2;
        const double this_energy = this_onsite_energy+this_hopping_energy;
        total_onsite_energy  += this_onsite_energy;
        total_hopping_energy += this_hopping_energy;
        total_energy  += this_energy;
        total_energy2 += this_energy*this_energy;
        total_vertex_number  -= this_hopping_energy;  // due to beta=1 (intrinsically normalized...)
        total_vertex_number2 += this_hopping_energy * this_hopping_energy; 
      }
      total_onsite_energy  /= beta;
      total_hopping_energy /= beta;
      total_energy  /= beta;
      total_energy2 /= (beta*beta); 

      std::cout << " ... Measuring " 
                << " ... N = " << total_density 
                << " ... E = " << total_energy
                ;

      measurements["Total Particle Number"]     << static_cast<double>(total_density);
      measurements["Total Particle Number^2"]   << static_cast<double>(total_density * total_density);
      measurements["Density"]                   << static_cast<double>(total_density)/num_sites();
      measurements["Density^2"]                 << static_cast<double>(total_density2)/num_sites();
      measurements["Hopping Energy"]            << total_hopping_energy;
      measurements["Onsite Energy"]             << total_onsite_energy;
      measurements["Energy"]                    << total_energy;
      measurements["Energy^2"]                  << total_energy * total_energy;
      measurements["Hopping Energy Density"]    << total_hopping_energy/num_sites();
      measurements["Onsite Energy Density"]     << total_onsite_energy/num_sites();
      measurements["Energy Density"]            << total_energy/num_sites();
      measurements["Energy Density^2"]          << total_energy2/num_sites();
      measurements["Total Vertex Number"]       << total_vertex_number;
      measurements["Total Vertex Number^2"]     << total_vertex_number2;

      if (measure_winding_number_)
      {
        std::vector<double> winding_number(dimension(), 0.);
        for (bond_iterator it=bonds().first; it!=bonds().second; ++it)
        {
          const int net_number = wl.net_number_of_directed_hops(source(*it),target(*it));
          const std::vector<double> vec = bond_vector_relative(*it);
          for (int i=0; i<dimension(); ++i)
            winding_number[i] += net_number * vec[i];
        }
        measurements["Winding Number^2"]         << alps::numeric::sq(std::vector<double>(winding_number.begin(), winding_number.end()));
      }      

      if (measure_local_density_) 
        measurements["Local Density"]            << _states_cache;

      if (measure_local_density2_)
        measurements["Local Density^2"]          << _states2_cache;

      // regarding on-fly measurements
      measurements["Green's Function Onsite"]    << green_onsite/num_sites();
      measurements["Green's Function Neighbors"] << green_neighbors/num_sites();
      measurements["Zero Momentum Density"]      << zero_momentum_density/num_sites();

      if (measure_green_function_)
        measurements["Green's Function"]         << green/num_sites();
      if (measure_momentum_density_)
        measurements["Momentum Density"]         << momentum_density/num_sites();

      if (finite_tof)
      {
        measurements["Zero TOF Image"]           << zero_momentum_density_tof/num_sites();

        if (measure_green_function_)
          measurements["TOF Green's Function"]   << green_tof/num_sites();
        if (measure_tof_image_)
          measurements["TOF Image"]              << momentum_density_tof/num_sites();
      }
    }

bool
  directed_worm_algorithm
    ::wormhead_propagates_till_collision_with_wormtail(unsigned short wormpair_state_, std::pair<neighbor_iterator,neighbor_iterator> const & neighbors_, bool const & measure_)
    {
      ++_propagation_counter;

      // Step 1: Wormhead movement
      const double _time2next     = worm.time2next();
      const double _time2wormtail = worm.time2wormtail();

      const double _onsite_energy_relative = onsite_energy_relative(worm.site(), worm.state(), worm.forward(), worm.creation());

      double _deltatime = -std::log(1.- random())/_onsite_energy_relative;

      const bool _halted = (_deltatime >= _time2next);
      if (_halted)
        _deltatime = _time2next;

      double _newtime = worm.forward() ? worm.time() + _deltatime : worm.time() - _deltatime;
      bool   _winding_over_time = (_newtime < 0. || _newtime >= 1.);

      if      (!_halted && _winding_over_time)
        _newtime = wormpair::modone(_newtime);
      else if (_halted)
        _newtime = worm.forward() ? worm.next_time()-std::numeric_limits<double>::epsilon() : worm.next_time()+std::numeric_limits<double>::epsilon();

      worm.wormhead_moves_to_new_time(_newtime,_winding_over_time);
      
#ifdef DEBUGMODE
//      if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER 
//         && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
//         ) 
      {
        std::cout << "\n\n"
                  << "\nPropagation step : " << _propagation_counter
                  << "\n------------------------------------"
                  << "\nAttempting to move wormhead by a time of " << _deltatime << " ( " << _time2next << " / " << _time2wormtail << " ) "
                  << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                  << "\n\n";
        //std::cin.get();
        std::cout << "\nWormhead moves to time " << worm.time() << " ( " << (_halted ? "HALTED" : "NOT HALTED") << " ) "
                  << "\n============================================================================"
                  << "\n" << wl
                  << "\n\n"
                  << "\nWormpair"
                  << "\n--------"
                  << "\n"
                  << worm
                  << "\n"
                  << "\nEnded (Wormhead movement)!"
                  << "\n--------------------------"
                  << "\n\n";
        //std::cin.get();
      }
#endif

      if (measure_ && (_deltatime > _time2wormtail) && (worm.site() != worm.wormtail_site()))   
      {
        using alps::numeric::operator+;
        using alps::numeric::operator+=;
        using alps::numeric::operator*;

        if (worm.neighbor2wormtail())
          green_neighbors += wormpair_state_;
        zero_momentum_density += wormpair_state_;

        int green_distance;
        if (measure_green_function_)
          green_distance = lattice().distance(worm.wormtail_site(), worm.site());

        double component_displacement;
        if (measure_momentum_density_ || measure_tof_image_)   // Achtung!!! Bitte "green_distance" nicht hier benutzen! Das ist SEHR langsam! (Tamama)
        {
          component_displacement = position(worm.site(),0) - position(worm.wormtail_site(),0);
          if (is_periodic_)
            if (component_displacement < 0)
              component_displacement += lattice_vector_[0][0];
        }

        if (measure_green_function_)
          green[green_distance] += wormpair_state_;
        if (measure_momentum_density_)
          for (unsigned int idx=0; idx<num_component_momenta(); ++idx)
            momentum_density[idx] += wormpair_state_*std::cos(component_momentum(idx)*component_displacement);

        if (finite_tof)
        {
          const double tof_phase = phase(worm.site())-phase(worm.wormtail_site()); 

          zero_momentum_density_tof += wormpair_state_ * std::cos(tof_phase);

          if (measure_green_function_)
            green_tof[green_distance] += wormpair_state_ * std::cos(tof_phase);
          if (measure_tof_image_)
            for (unsigned int idx=0; idx<num_component_momenta(); ++idx)
              momentum_density_tof[idx] += wormpair_state_ * std::cos(component_momentum(idx)*component_displacement + tof_phase);
        }
      }

      if (!_halted)
      {
        // Step 2A: Wormhead either inserts vertex or rebounces
        insert_jump_or_bounce(_onsite_energy_relative,neighbors_); 
        return true;
      }

      else
      {
        // Step 2B: Wormhead touches wormtail (offdiagonal process ends)
        if (worm.wormhead_reaches_wormtail())
           return false;

        // Step 2C: Wormhead crosses kink
        if (worm.wormhead_is_same_type_as_next())
        {
          worm.wormhead_crosses_vertex();
#ifdef DEBUGMODE
//          if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
//             && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
//             )
          {
            std::cout << "\n\n"
                      << "\nWormhead crosses kink"
                      << "\n============================================================================"
                      << "\n" << wl
                      << "\n\n"
                      << worm
                      << "\n\n"
                      << "\nEnded (Wormhead crosses kink)!"
                      << "\n--------------------------"
                      << "\n\n";
            //std::cin.get();
          }
#endif
          return true;
        }
        // Step 2D: Wormhead either deletes/relinks vertex or rebounces
        else
        {
          delete_relink_jump_or_bounce(neighbors_);
          return true;
        }
      }
    }

void
  directed_worm_algorithm
    ::insert_jump_or_bounce (double const onsite_energy_relative_, std::pair<neighbor_iterator,neighbor_iterator> const & neighbors_)
    {
#ifndef HEATBATH_ALGORITHM
      neighbor_bond_iterator it=neighbor_bonds(worm.site()).first + random()*num_neighbors(worm.site());
      location_type _neighborlocation  = wl.location(target(*it), worm.time());
      const unsigned short _targetstate = (_neighborlocation.second-1)->state();
      const bool _increasing = worm.increasing();
      if (  ( worm.increasing() && _targetstate == state_maximum[site_site_type[target(*it)]])
         || (!worm.increasing() && _targetstate == state_minimum[site_site_type[target(*it)]])
         )  // new configuration is not possible... only bouncing allowed...
      {
        worm.wormhead_turns_around();
        ++_propagation_failure_counter;
        return;
      }
      else
      {
        const double _weight_new = hopping_energy(*it, _targetstate, worm.increasing());
        if (_weight_new >= onsite_energy_relative_ || random() < (_weight_new/onsite_energy_relative_))  // insert vertex
        {
          worm.wormhead_inserts_vertex_and_jumps_to_new_site(_neighborlocation);
          worm.set_neighbor2wormtail(std::find(neighbors_.first, neighbors_.second, worm.site()) != neighbors_.second);
          return;
        }
        else  // or bounce
        {
          worm.wormhead_turns_around();
          ++_propagation_failure_counter;
          return;
        }
      }

#else
      _neighborlocations_cache.clear();
      _cummulative_weights_cache.clear();

      double _total_weight = 0.;
      for (neighbor_bond_iterator _neighborbondit=neighbor_bonds(worm.site()).first; _neighborbondit!=neighbor_bonds(worm.site()).second; ++_neighborbondit)
      {
        location_type  _neighborlocation  = wl.location(target(*_neighborbondit), worm.time());
        _neighborlocations_cache.push_back(_neighborlocation);                                              

        const unsigned short _targetstate = (_neighborlocation.second-1)->state();
        const bool _increasing = worm.increasing();
        if (  ( worm.increasing() && _targetstate == state_maximum[site_site_type[target(*_neighborbondit)]])
           || (!worm.increasing() && _targetstate == state_minimum[site_site_type[target(*_neighborbondit)]])
           )
           _cummulative_weights_cache.push_back(_total_weight);
        else
        {
           _total_weight += hopping_energy(*_neighborbondit, _targetstate, worm.increasing());
           _cummulative_weights_cache.push_back(_total_weight);
        }
      }
   
      _total_weight += onsite_energy_relative_;
      _cummulative_weights_cache.push_back(_total_weight);

      const int which_neighbor = std::distance(_cummulative_weights_cache.begin(), std::lower_bound(_cummulative_weights_cache.begin(),_cummulative_weights_cache.end(), random()*_total_weight));

      if (which_neighbor == _neighborlocations_cache.size())
      {
        worm.wormhead_turns_around();
        ++_propagation_failure_counter;
        return;    
      }

      worm.wormhead_inserts_vertex_and_jumps_to_new_site(_neighborlocations_cache[which_neighbor]);
      worm.set_neighbor2wormtail(std::find(neighbors_.first, neighbors_.second, worm.site()) != neighbors_.second);
      return;
#endif
    }

void
  directed_worm_algorithm
    ::delete_relink_jump_or_bounce(std::pair<neighbor_iterator,neighbor_iterator> const & neighbors_)
    {
      location_type _sourcelocation = wl.location(worm.next_partnersite(), worm.next_time());

#ifndef HEATBATH_ALGORITHM
      const std::pair<neighbor_bond_iterator, neighbor_bond_iterator> _neighbor_bonds = neighbor_bonds(worm.next_partnersite());

      double _weight;
      for (neighbor_bond_iterator it=_neighbor_bonds.first; it!=_neighbor_bonds.second; ++it)
        if (target(*it) == worm.site())
        {
          _weight = hopping_energy(*it, (wl.location(target(*it), worm.next_time()).second-1)->state(), !worm.increasing());
          break;
        }
     
      neighbor_bond_iterator it=_neighbor_bonds.first + random()*num_neighbors(worm.next_partnersite());

      if (target(*it) == worm.site())    // delete vertex or bounce
      {
        const double _weight_new = onsite_energy_relative(worm.next_partnersite(),_sourcelocation.second->state(),!worm.forward(),worm.creation());

        if (_weight_new >= _weight || random() < (_weight_new/_weight))
        {
          worm.wormhead_deletes_vertex_and_jumps_to_new_site(_sourcelocation);
          worm.set_neighbor2wormtail(std::find(neighbors_.first, neighbors_.second, worm.site()) != neighbors_.second);
          return;
        }
        else
        {
          worm.wormhead_turns_around();
          ++_propagation_failure_counter;
          return;
        }
      }
      else     // relink vertex
      {
        location_type _neighborlocation  = wl.location(target(*it), worm.next_time());
        const unsigned short _targetstate = (_neighborlocation.second-1)->state();
        const bool _increasing = !worm.increasing();
        if (  ( _increasing && _targetstate == state_maximum[site_site_type[target(*it)]])
           || (!_increasing && _targetstate == state_minimum[site_site_type[target(*it)]])
           )
        {
          worm.wormhead_turns_around();
          ++_propagation_failure_counter;
          return;
        }
        else
        {
          const double _weight_new = hopping_energy(*it, _targetstate, _increasing);

          if (_weight_new >= _weight || random() < (_weight_new/_weight))
          {
            worm.wormhead_relinks_vertex_and_jumps_to_new_site(_sourcelocation, _neighborlocation);
            worm.set_neighbor2wormtail(std::find(neighbors_.first, neighbors_.second, worm.site()) != neighbors_.second);
            return;
          }
          else
          {
            worm.wormhead_turns_around();
            ++_propagation_failure_counter;
            return;
          }
        }
      }

#else
      _neighborlocations_cache.clear();
      _cummulative_weights_cache.clear();

      double _total_weight = 0.;
      for (neighbor_bond_iterator _neighborbondit=neighbor_bonds(worm.next_partnersite()).first; _neighborbondit!=neighbor_bonds(worm.next_partnersite()).second; ++_neighborbondit)
      {
        location_type _neighborlocation  = wl.location(target(*_neighborbondit), worm.next_time());

        _neighborlocations_cache.push_back(_neighborlocation);

        const unsigned short _targetstate = (_neighborlocation.second-1)->state();
        const bool _increasing = !worm.increasing();
        if (  ( _increasing && _targetstate == state_maximum[site_site_type[target(*_neighborbondit)]])
           || (!_increasing && _targetstate == state_minimum[site_site_type[target(*_neighborbondit)]])
           )
           _cummulative_weights_cache.push_back(_total_weight);
        else
        {
           _total_weight += hopping_energy(*_neighborbondit, _targetstate, _increasing);
           _cummulative_weights_cache.push_back(_total_weight);
        }
      }

      _total_weight += onsite_energy_relative(worm.next_partnersite(),_sourcelocation.second->state(),!worm.forward(),worm.creation());
      _cummulative_weights_cache.push_back(_total_weight);

      const int which_neighbor = std::distance(_cummulative_weights_cache.begin(), std::lower_bound(_cummulative_weights_cache.begin(),_cummulative_weights_cache.end(), random()*_total_weight));

      if (which_neighbor == _neighborlocations_cache.size())
      {
        worm.wormhead_deletes_vertex_and_jumps_to_new_site(_sourcelocation);
        worm.set_neighbor2wormtail(std::find(neighbors_.first, neighbors_.second, worm.site()) != neighbors_.second);
        return;
      }

      if (_neighborlocations_cache[which_neighbor].first->begin()->siteindicator() == worm.site())
      {
        worm.wormhead_turns_around();
        ++_propagation_failure_counter;
        return;
      }

      worm.wormhead_relinks_vertex_and_jumps_to_new_site(_sourcelocation, _neighborlocations_cache[which_neighbor]);  
      worm.set_neighbor2wormtail(std::find(neighbors_.first, neighbors_.second, worm.site()) != neighbors_.second);
      return;
#endif
    }

// ==================================================
// main 
// ==================================================

int main(int argc, char** argv)
{
  // read options from command line
  alps::mcoptions options(argc, argv);

  // initialize simulation and load the parameters from input file
  alps::hdf5::archive simulation_input = alps::hdf5::archive(options.input_file);
  directed_worm_algorithm simulation(simulation_input);

  // try loading worldlines and measurement from the input file
  try {
    simulation.load(simulation_input);
  }
  catch(...) {}
 
  // run the simulation
  simulation.run(alps::stop_callback(options.time_limit));

  // save worldlines and measurement to output file
  alps::hdf5::archive simulation_output = alps::hdf5::archive(options.output_file, 'a');
  simulation.save(simulation_output);

  return 0;
}
