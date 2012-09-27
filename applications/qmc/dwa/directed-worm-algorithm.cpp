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

#include "directed-worm-algorithm.hpp"


namespace alps {
namespace applications {

void
  directed_worm_algorithm
    ::print_copyright (std::ostream & out)
    {
      out << "\n\n"
          << "\n/*****************************************************************************"
          << "\n*"
          << "\n* ALPS Project Applications: Directed Worm Algorithm"
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

      out << "Parameters\n"
          << "==========\n\n"
          << parms
          << "\n\n";

      out << "Simulation\n"
          << "==========\n\n"
          << "\t Total thermalization sweeps       : " << _total_thermalization_sweeps   << "\n"
          << "\t Total sweeps                      : " << _total_sweeps                  << "\n"
          << "\t Sweep per measurement             : " << _sweep_per_measurement         << "\n"
          << "\t Sweep counter                     : " << _sweep_counter                 << "\n"
          << "\t Sweep failure counter             : " << _sweep_failure_counter         << "\n"
          << "\t Propagation counter               : " << _propagation_counter           << "\n"
          << "\t Propagation failure counter       : " << _propagation_failure_counter   << "\n"
          << "\n\n";

      out << "Lattice\n"
          << "=======\n\n"
          << "Name              : " << static_cast<std::string>(parms["LATTICE"]) << "\n"
          << "Dimension         : " << dimension()                                << "\n"
          << "Total Sites/Bonds : " << num_sites() << "\t/\t" << num_bonds()      << "\n"
          << "Periodic          : " << !inhomogeneous() << "\n"
#ifdef DEBUGMODE
          << "Graph             : " << "\n"
          << graph()
#endif
          << "\n\n";
    }

directed_worm_algorithm
  ::directed_worm_algorithm
    ( const alps::ProcessList &  processes
    , const alps::Parameters  &  parameters
    , int                        processnode
    )
    : QMCRun<>(processes, parameters, processnode)

    // regarding MC simulation
    , _sweep_counter               (0)
    , _sweep_failure_counter       (0)
    , _propagation_counter         (0)
    , _propagation_failure_counter (0)
    , _total_thermalization_sweeps (parameters.value_or_default("THERMALIZATION",0))
    , _total_sweeps                (parameters.value_or_default("TOTAL_SWEEPS",10000000))   
    , _sweep_per_measurement       (parameters.value_or_default("SKIP",1))

    // regarding worldline
    , worldline (num_sites(), parameters.value_or_default("LINE_CAPACITY",1))

    // regarding experiment
    , is_experiment (static_cast<bool>(parameters.value_or_default("experiment",false)))

    // regarding measurements
    , _measure_now (false)

    // regarding on-fly measurements
    , _on_fly_measurements_reinitialized (false)
  {
    initialize_site_states();
    initialize_hamiltonian();
    initialize_measurements();
  }

void
  directed_worm_algorithm
    ::start()
    {
      print_simulation (std::cout);
#ifdef DEBUGMODE
      print_hamiltonian(std::cout);
      print_worldline  (std::cout);
#endif
    }

void
  directed_worm_algorithm
    ::initialize_hamiltonian()
    {
      // setup state_minimum and state_maximum
      for (int i=0; i<maximum_sitetype+1; ++i)
      {
        state_minimum.push_back(alps::applications::numeric::iround(diagonal_matrix_element.begin()->second[i].front()));
        state_maximum.push_back(alps::applications::numeric::iround(diagonal_matrix_element.begin()->second[i].back()));
      }

      // this part takes care of experimental details
      if (is_experiment && model().name() == "trapped boson Hubbard")
      {
        phase = std::vector<double>(dimension(), 
          1e-8 * static_cast<double>(parms["mass"])*amu*alps::numeric::sq(static_cast<double>(parms["lambda"])) 
          / (8*static_cast<double>(parms["time_of_flight"])*hbar));
      }

      // this part takes care of disorders
      if (model().name() == "trapped boson Hubbard")
      {
        percentage_error_hopping_t.resize(num_bonds(), 1. + (2*random()-1)*static_cast<double>(parms.value_or_default("t_disorder",0)));
        percentage_error_onsite_U .resize(num_sites(), 1. + (2*random()-1)*static_cast<double>(parms.value_or_default("U_disorder",0)));
        percentage_error_onsite_mu.resize(num_sites(), 1. + (2*random()-1)*static_cast<double>(parms.value_or_default("mu_disorder",0)));
      }

      // iterate all sites and build onsite matrix 
      for (site_iterator it = sites().first; it != sites().second; ++it) {
        unsigned int this_site_type = inhomogeneous_site_type(*it);
        if (this_site_type >= onsite_matrix.size())
        { 
          using boost::numeric::operators::operator*;
          onsite_matrix.resize(this_site_type+1);
          onsite_matrix[this_site_type] = onsite_hamiltonian(*it) * beta;
        }
      }

      // iterate all bonds and build site 1-up matrix, site 1-down matrix and bond strength matrix
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

      if (model().name() == "trapped boson Hubbard")
      {
        using boost::numeric::operators::operator*;
        bond_strength_matrix = bond_strength_matrix * percentage_error_hopping_t; 
      }
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
            << "States     : " << (diagonal_matrix_element.begin()->second)[site_site_type[*it]]  << "\n\n"
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
    ::print_worldline(std::ostream & out) const
    {
      out << "Worldline:\n"
          << "==========\n\n"
          << worldline << "\n"
          << "==========\n\n";
    }

void
  directed_worm_algorithm
    ::initialize_measurements()
    {
      // regarding measurements
      measurements 
        << alps::IntTimeSeriesObservable    ("Total Particle Number")      // The entire time series is stored for thermalization indication
        << alps::IntObservable              ("Total Particle Number^2")
        << alps::RealObservable             ("Density")
        << alps::RealObservable             ("Density^2")
        << alps::RealObservable             ("Hopping Energy")
        << alps::RealObservable             ("Onsite Energy")
        << alps::RealObservable             ("Energy")
        << alps::RealObservable             ("Energy^2")
        << alps::RealObservable             ("Energy Density")
        << alps::RealObservable             ("Energy Density^2")
        << alps::SimpleIntVectorObservable  ("Local Density")
        << alps::SimpleIntVectorObservable  ("Local Density^2")
        << alps::SimpleRealVectorObservable ("Local Hopping Energy")
        << alps::SimpleRealVectorObservable ("Local Onsite Energy")
        << alps::SimpleRealVectorObservable ("Local Energy")
        << alps::SimpleRealVectorObservable ("Local Energy^2")
      ;

      // regarding on-fly measurements
      measurements 
        << alps::SimpleRealVectorObservable ("Green's Function")
        << alps::RealObservable             ("Onsite Green's Function")
        << alps::RealVectorObservable       ("N.N. Green's Function")
        << alps::RealObservable             ("Condensate Fraction")
      ;

      if (is_experiment && model().name() == "trapped boson Hubbard")
      {
      measurements
        << alps::SimpleRealVectorObservable ("TOF Green's Function")
        << alps::RealObservable             ("TOF Onsite Green's Function")
        << alps::RealVectorObservable       ("TOF N.N. Green's Function")
        << alps::RealObservable             ("TOF Condensate Fraction")
      ;
      }

      if (!inhomogeneous())
      measurements << alps::IntVectorObservable ("Winding Number^2");
    }

void
  directed_worm_algorithm
    ::reinitialize_on_fly_measurements()
    {
      // density matrix
      green.clear();
      green.resize(num_sites(), 0.);

      if (is_experiment && model().name() == "trapped boson Hubbard")
      {
      green_tof.clear();
      green_tof.resize(num_sites(), 0.);
      }

      // winding number
      winding_number.clear();
      winding_number.resize(dimension(), 0);
    }

std::vector<double> 
  directed_worm_algorithm
    ::onsite_hamiltonian(const site_descriptor& site) 
    {
      unsigned int this_site_type = site_type(site);

      std::vector<alps::SiteTermDescriptor> all_site_terms = model().site_terms();
      std::map<std::string, alps::SiteTermDescriptor> this_site_terms;
      for (std::vector<alps::SiteTermDescriptor>::iterator it=all_site_terms.begin(); it!=all_site_terms.end(); ++it)
        if (it->match_type(site_site_type[this_site_type]))
          this_site_terms.insert(std::make_pair(it->name(), *it));

      std::vector<double> this_onsite_hamiltonian;

      // For trapped boson Hubbard model
      if (model().name() == "trapped boson Hubbard")
      {
        { // U term
          boost::multi_array<double,2> this_onsite_matrix =
            alps::get_matrix( double(), this_site_terms.find("U term")->second
                            , model().basis().site_basis(this_site_type), parms);
          std::vector<double> this_onsite_matrix_diagonal;
          for (state_type i=0; i < this_onsite_matrix.shape()[0]; ++i)
            this_onsite_matrix_diagonal.push_back(this_onsite_matrix[i][i]);
          using boost::numeric::operators::operator+;
          using alps::numeric::operator*;
          this_onsite_matrix_diagonal = this_onsite_matrix_diagonal * percentage_error_onsite_U[site];
          if (this_onsite_hamiltonian.empty())
            this_onsite_hamiltonian = this_onsite_matrix_diagonal;
          else
            this_onsite_hamiltonian = this_onsite_hamiltonian + this_onsite_matrix_diagonal;
        }
        { // mu term
          boost::multi_array<double,2> this_onsite_matrix =
            alps::get_matrix( double(), this_site_terms.find("mu term")->second
                            , model().basis().site_basis(this_site_type), parms);
          std::vector<double> this_onsite_matrix_diagonal;
          for (state_type i=0; i < this_onsite_matrix.shape()[0]; ++i)
            this_onsite_matrix_diagonal.push_back(this_onsite_matrix[i][i]);
          using boost::numeric::operators::operator+;
          using alps::numeric::operator*;
          this_onsite_matrix_diagonal = this_onsite_matrix_diagonal * percentage_error_onsite_mu[site];
          if (this_onsite_hamiltonian.empty())
            this_onsite_hamiltonian = this_onsite_matrix_diagonal;
          else
            this_onsite_hamiltonian = this_onsite_hamiltonian + this_onsite_matrix_diagonal;
        }
        { // trapping term
          using boost::numeric::operators::operator+;
          using alps::numeric::operator*;
          using alps::applications::numeric::norm;
          this_onsite_hamiltonian = this_onsite_hamiltonian 
                                  + (static_cast<double>(parms["K"]) * norm(position(site))) * diagonal_matrix_element["n"][site_site_type[site]];
        }
        return this_onsite_hamiltonian;
      }

      // For any other models
      {
        for (std::map<std::string, alps::SiteTermDescriptor>::iterator it=this_site_terms.begin(); it!=this_site_terms.end(); ++it)
        {
          boost::multi_array<double,2> this_onsite_matrix =
            alps::get_matrix( double(), it->second
                            , model().basis().site_basis(this_site_type), parms);
          std::vector<double> this_onsite_matrix_diagonal;
          for (state_type i=0; i < this_onsite_matrix.shape()[0]; ++i)
            this_onsite_matrix_diagonal.push_back(this_onsite_matrix[i][i]);
          using boost::numeric::operators::operator+;
          if (this_onsite_hamiltonian.empty())
            this_onsite_hamiltonian = this_onsite_matrix_diagonal;
          else
            this_onsite_hamiltonian = this_onsite_hamiltonian + this_onsite_matrix_diagonal;
        }
        return this_onsite_hamiltonian;
      }
    }   

inline boost::multi_array<double,4>
  directed_worm_algorithm
    ::bond_hamiltonian(const bond_descriptor & bond)
    {
      return alps::get_matrix( double(), model().bond_term(bond_type(bond))
                             , model().basis().site_basis(site_type(source(bond)))
                             , model().basis().site_basis(site_type(target(bond)))
                             , parms);
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

      std::set<std::string> site_operators = model().bond_term(bond_type(bond)).operator_names(parms);

      for (std::set<std::string>::iterator site_operator_it = site_operators.begin(); site_operator_it != site_operators.end(); ++site_operator_it)
      {
        boost::multi_array<double,2> this_sourcesite_operator_matrix 
          = alps::get_matrix(double(), SiteOperator(*site_operator_it), model().basis().site_basis(source(bond)), parms);
        boost::multi_array<double,2> this_targetsite_operator_matrix
          = alps::get_matrix(double(), SiteOperator(*site_operator_it), model().basis().site_basis(target(bond)), parms);

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

inline double
  directed_worm_algorithm
    ::onsite_energy_relative (Kink const & kink_, bool const forward_) const 
    {
      double energy_offset = 0.1;   // Note: Ask Lode for a better automation method

      double _onsite_energy       = onsite_energy(kink_.site(), kink_.state());
      double _onsite_energy_after = onsite_energy(kink_.site(), kink_.state_after());

      return (forward_ ? _onsite_energy       - std::min(_onsite_energy, _onsite_energy_after) + energy_offset
                       : _onsite_energy_after - std::min(_onsite_energy, _onsite_energy_after) + energy_offset); 
    }

inline double
  directed_worm_algorithm
    ::hopping_energy (bond_descriptor const & bond_, State const targetstate_, bool const increasing_) const 
    {  
      return (increasing_ ? (bond_strength_matrix[index(bond_)] * site_oneup_matrix[site_site_type[target(bond_)]][targetstate_] * site_onedown_matrix[site_site_type[target(bond_)]][targetstate_+1]) : (bond_strength_matrix[index(bond_)] * site_onedown_matrix[site_site_type[target(bond_)]][targetstate_] * site_oneup_matrix[site_site_type[target(bond_)]][targetstate_-1]));
    }

void
  directed_worm_algorithm
    ::dostep()
    {
      //
      // 1) The dostep() function attempts to do an update on the (closed) worldline configuration through a series of worm updates.
      // 2) It starts off with wormpair creation. If unsuccessful, dostep() will be return-ed.
      // 3) Then, by moving the wormhead "detailed-balance"-ly, vertices are created, deleted, or relinked according to the parameters characterized by the boson hubbard model.
      // 4) It ends off by annihilating the wormpair, and a new worldline configuration arises most probably.
      // 5) dostep() will then return-ed.
      //

      ++_sweep_counter;

      // are we ready to perform measurements?
      if (_sweep_counter % _sweep_per_measurement == 0)
        _measure_now = true;

      // shall we reinitialize on fly measurements?
      if (!_on_fly_measurements_reinitialized)
      {
        reinitialize_on_fly_measurements();
        _on_fly_measurements_reinitialized = true;
      }

#ifdef DEBUGMODE
      if (_sweep_counter >= DEBUGMODE_START_COUNTER && _sweep_counter <= DEBUGMODE_END_COUNTER) {
        std::cout << "\n\nSweep : " << _sweep_counter << "\t(Propagation sweep: " << _propagation_counter << ")"
                  << "\n============"
                  << "\n\n";
      }
#endif

      // Step 1: Attempt to create a worm-pair
      bool _wormpair_has_found_suitable_insertion_location = false;
      while (!_wormpair_has_found_suitable_insertion_location)
      {
        Site _currentsite       = randomsite();
        Time _currenttime       = randomtime();
        bool _currentcreation   = randombool();
        bool _currentforward    = randombool();

        _wormpair_has_found_suitable_insertion_location = worldline.unoccupied(_currentsite, _currenttime);
        if (_wormpair_has_found_suitable_insertion_location) 
        {
          State _state = worldline.state(_currentsite, _currenttime);

          // Check if state is at minimum or maximum ?
          if (  (_state == state_minimum[site_site_type[_currentsite]] && !alps::applications::increasing(_currentforward, _currentcreation)) 
             || (_state == state_maximum[site_site_type[_currentsite]] && alps::applications::increasing(_currentforward, _currentcreation))  
             )
          {
            ++_sweep_failure_counter;      // worm insertion is unsuccessful          
            green[0] += _state;
            if (is_experiment && model().name() == "trapped boson Hubbard")
            green_tof[0] += _state;
#ifdef DEBUGMODE
            if (_sweep_counter >= DEBUGMODE_START_COUNTER && _sweep_counter <= DEBUGMODE_END_COUNTER) {
            std::cout << "\n\n"
                      << "\nWormpair creation fails!" 
                      << "\n========================"
                      << worldline
                      << "\nEnded (Wormpair creation fails!)"
                      << "\n--------------------------------"
                      << "\n\n";
            }
#endif
            if (_measure_now) {
              perform_diagonal_measurements();
              _measure_now = false;
              _on_fly_measurements_reinitialized = false;
            }
            return;
          }

          State _newstate = alps::applications::increasing(_currentforward, _currentcreation) ? _state + 1 : _state - 1;

          wormpair = std::make_pair<Kink, Kink>
               ( Kink(true, _currentcreation,  _currentforward,  _newstate, _currenttime, _currentsite, _currentsite, false)     // wormhead
               , Kink(true, !_currentcreation, !_currentforward, _state,    _currenttime, _currentsite, _currentsite, true)      // wormtail
               );
          if (!_currentforward)
            swap_state(wormpair.first, wormpair.second);

          worldline.wormpair_insertion(wormpair);
        }
      }
#ifdef DEBUGMODE
      if (_sweep_counter >= DEBUGMODE_START_COUNTER && _sweep_counter <= DEBUGMODE_END_COUNTER) {
        std::cout << "\n\n"
                  << "\nWormpair creation!" 
                  << "\n=================="
                  << worldline
                  << "\n"
                  << "\n\n"
                  << "\nWormhead"
                  << "\n--------"
                  << wormpair.first
                  << "\n"
                  << "\nWormtail"
                  << "\n--------"
                  << wormpair.second
                  << "\n"
                  << "\nEnded (Wormpair creation!)"
                  << "\n--------------------------"
                  << "\n\n";
      }
#endif

      // Step 2: Propagates the wormhead until it collides with wormtail
      //         Note -- density matrix is being measured everytime wormhead passes through equal time plane
      State  _wormpair_state = (wormpair.second.creation() ? wormpair.first.state() : wormpair.second.state() );
      while(wormhead_propagation_until_halted_by_wormtail(_wormpair_state));      

#ifdef DEBUGMODE
       if (_sweep_counter >= DEBUGMODE_START_COUNTER && _sweep_counter <= DEBUGMODE_END_COUNTER) {
       std::cout << "\n\n"
                 << "\nMeasuring Green's function... (annilation)"
                 << "\n\n";
       }
#endif

      wormpair.first = worldline.wormhead(); 

      green[0] += (wormpair.first.forward() ? wormpair.first.state() : wormpair.first.state_after());
      if (is_experiment && model().name() == "trapped boson Hubbard")
      green_tof[0] += (wormpair.first.forward() ? wormpair.first.state() : wormpair.first.state_after());

      worldline.wormpair_removal();
#ifdef DEBUGMODE
      if (_sweep_counter >= DEBUGMODE_START_COUNTER && _sweep_counter <= DEBUGMODE_END_COUNTER) {
      std::cout << "\n\n"
                << "\nWormpair annihilation!"
                << "\n======================"
                << worldline
                << "\nEnded (Wormpair annihilation!)"
                << "\n------------------------------------------------------"
                << "\n\n";
      }
#endif
      if (_measure_now) {
        perform_diagonal_measurements();
        _measure_now = false;
        _on_fly_measurements_reinitialized = false;
      }
      return;
    }

void
  directed_worm_algorithm
    ::perform_diagonal_measurements()
    {
      // regarding measurements
      using boost::numeric::operators::operator+;
      using boost::numeric::operators::operator*;
      using boost::numeric::operators::operator/;

      std::vector<int>    _states           = alps::applications::numeric::vector_cast<int>(worldline.states());
      std::vector<int>    _states2          = _states * _states;
      std::vector<double> _hopping_energies = alps::applications::numeric::vector_cast<double>(worldline.vertices())/(-2.*beta);
      std::vector<double> _onsite_energies  = onsite_energies(_states);
      std::vector<double> _energies         = _hopping_energies + _onsite_energies;
      std::vector<double> _energies2        = _energies * _energies;

      int total_particle_number = std::accumulate(_states.begin(), _states.end(), 0);

      double total_hopping_energy  = std::accumulate(_hopping_energies.begin(), _hopping_energies.end(), 0.); 
      double total_onsite_energy   = std::accumulate(_onsite_energies.begin(), _onsite_energies.end(), 0.);     
      double total_energy          = total_hopping_energy + total_onsite_energy;

      std::cout << "Measuring " 
//                << _sweep_counter  
                << " ... N = " << total_particle_number 
//                << " ... Eo= " << total_hopping_energy
//                << " ... Ed= " << total_onsite_energy
//                << " ... E = " << total_energy
                << "\n";

      measurements["Total Particle Number"]     << total_particle_number;
      measurements["Total Particle Number^2"]   << total_particle_number * total_particle_number;
      measurements["Density"]                   << static_cast<double>(total_particle_number)/num_sites();
      measurements["Density^2"]                 << static_cast<double>(std::accumulate(_states2.begin(), _states2.end(), 0))/num_sites();
      measurements["Hopping Energy"]            << total_hopping_energy;
      measurements["Onsite Energy"]             << total_onsite_energy;
      measurements["Energy"]                    << total_energy;
      measurements["Energy^2"]                  << total_energy * total_energy;
      measurements["Energy Density"]            << total_energy/num_sites();
      measurements["Energy Density^2"]          << std::accumulate(_energies2.begin(), _energies2.end(), 0.)/num_sites();
      measurements["Local Density"]             << alps::numeric::vector2valarray(_states);
      measurements["Local Density^2"]           << alps::numeric::vector2valarray(_states2);
      measurements["Local Hopping Energy"]      << alps::numeric::vector2valarray(_hopping_energies);
      measurements["Local Onsite Energy"]       << alps::numeric::vector2valarray(_onsite_energies);
      measurements["Local Energy"]              << alps::numeric::vector2valarray(_energies);
      measurements["Local Energy^2"]            << alps::numeric::vector2valarray(_energies2);

      // regarding on-fly measurements
      measurements["Green's Function"]          << alps::numeric::vector2valarray(green);

      std::vector<double> _nn_green;
      _nn_green.reserve(dimension());
      for (int i=0; i<dimension(); ++i)
        _nn_green.push_back(green[nearest_neighbor_idx(0)[i]]);

      measurements["Onsite Green's Function"]   << green[0];
      measurements["N.N. Green's Function"]     << alps::numeric::vector2valarray(_nn_green);

      if (!inhomogeneous())
      measurements["Condensate Fraction"]       << std::accumulate(green.begin(), green.end(), 0.);
      else
      measurements["Condensate Fraction"]       << (2. * std::accumulate(green.begin(), green.end(), 0.)) - green[0];

      if (is_experiment && model().name() == "trapped boson Hubbard")
      {
      measurements["TOF Green's Function"]          << alps::numeric::vector2valarray(green_tof);

      std::vector<double> _nn_green_tof;
      _nn_green_tof.reserve(dimension());
      for (int i=0; i<dimension(); ++i)
        _nn_green_tof.push_back(green_tof[nearest_neighbor_idx(0)[i]]);

      measurements["TOF Onsite Green's Function"]   << green_tof[0];
      measurements["TOF N.N. Green's Function"]     << alps::numeric::vector2valarray(_nn_green_tof);

      if (!inhomogeneous())
      measurements["TOF Condensate Fraction"]       << std::accumulate(green_tof.begin(), green_tof.end(), 0.);
      else
      measurements["TOF Condensate Fraction"]       << (2. * std::accumulate(green_tof.begin(), green_tof.end(), 0.)) - green_tof[0];
      }

      if (!inhomogeneous())
      measurements["Winding Number^2"]              << alps::numeric::vector2valarray(winding_number*winding_number);
    }

bool
  directed_worm_algorithm
    ::wormhead_propagation_until_halted_by_wormtail(State const & wormpair_state_)
    {
      ++_propagation_counter;

      double _onsite_energy_relative = onsite_energy_relative(worldline.wormhead(), worldline.forward());

      // Step 1: Wormhead movement
      bool _halted;
      Time _deltatime       = -std::log(1.- randomtime())/_onsite_energy_relative;
      Time _deltatime_pause;
      Time _deltatime_halt;

      if (worldline.forward())
      {
        _deltatime_pause = alps::applications::numeric::mod(wormpair.second.time() - worldline.time());  
        if (worldline.wormhead_touches_end())
        {
          _deltatime_halt = 1.- worldline.time();
#ifdef DEBUGMODE
          if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER 
             && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
             ) 
          {
          std::cout << "\n\n"
                    << "\nPropagation step : " << _propagation_counter
                    << "\n------------------------------------"
                    << "\nAttempting to move the wormhead over a time interval of " << _deltatime << " ( " << _deltatime_pause << " / " << _deltatime_halt << " ) "
                    << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                    << "\n\n";

          }
#endif
          if (_deltatime < _deltatime_halt)
          {
            worldline.wormhead_moves_to_new_time_without_winding(worldline.time()+_deltatime);
            _deltatime_pause -= _deltatime;
            _halted = false;
          }
          else
          {
            worldline.wormhead_moves_to_new_time_with_winding(0.);
            _deltatime_pause -= _deltatime_halt;
            _deltatime       -= _deltatime_halt;

            _deltatime_halt = worldline.forwardlineit()->time() - worldline.time();

            if (_deltatime < _deltatime_halt)
            {
              worldline.wormhead_moves_to_new_time_without_winding(worldline.time()+_deltatime);
              _deltatime_pause -= _deltatime;
              _halted = false;
            }
            else
            {
              worldline.wormhead_moves_to_new_time_without_winding(worldline.time()+_deltatime_halt);
              _deltatime_pause -= _deltatime_halt;
              _halted = true;
            }
          }
        }
        else
        {
          _deltatime_halt = worldline.forwardlineit()->time() - worldline.time();
#ifdef DEBUGMODE
          if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
             && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
             ) 
          {
          std::cout << "\n\n"
                    << "\nPropagation step : " << _propagation_counter
                    << "\n------------------------------------"
                    << "\nAttempting to move the wormhead over a time interval of " << _deltatime << " ( " << _deltatime_pause << " / " << _deltatime_halt << " ) "
                    << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                    << "\n\n";
          }
#endif
          if (_deltatime < _deltatime_halt)
          {
            worldline.wormhead_moves_to_new_time_without_winding(worldline.time()+_deltatime);
            _deltatime_pause -= _deltatime;
            _halted = false;
          }
          else
          {
            worldline.wormhead_moves_to_new_time_without_winding(worldline.time()+_deltatime_halt);
            _deltatime_pause -= _deltatime_halt;
            _halted = true;
          }
        }
      }
      else  // backward
      {
        _deltatime_pause = alps::applications::numeric::mod(worldline.time() - wormpair.second.time());
        if (worldline.wormhead_touches_begin())
        {
          _deltatime_halt = worldline.time();
#ifdef DEBUGMODE
          if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
             && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
             ) 
          {
          std::cout << "\n\n"
                    << "\nPropagation step : " << _propagation_counter
                    << "\n------------------------------------"
                    << "\nAttempting to move the wormhead over a time interval of " << _deltatime << " ( " << _deltatime_pause << " / " << _deltatime_halt << " ) "
                    << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                    << "\n\n";
          }
#endif
          if (_deltatime < _deltatime_halt)
          {
            worldline.wormhead_moves_to_new_time_without_winding(worldline.time()-_deltatime);
            _deltatime_pause -= _deltatime;
            _halted = false;
          }
          else
          {
            worldline.wormhead_moves_to_new_time_with_winding(1.);
            _deltatime_pause -= _deltatime_halt;
            _deltatime       -= _deltatime_halt;

            _deltatime_halt = worldline.time() - worldline.backwardlineit()->time();

            if (_deltatime < _deltatime_halt)
            {
              worldline.wormhead_moves_to_new_time_without_winding(worldline.time()-_deltatime);
              _deltatime_pause -= _deltatime;
              _halted = false;
            }
            else
            {
              worldline.wormhead_moves_to_new_time_without_winding(worldline.time()-_deltatime_halt);
              _deltatime_pause -= _deltatime_halt;
              _halted = true;
            }
          }
        }
        else
        {
          _deltatime_halt = worldline.time() - worldline.backwardlineit()->time();
#ifdef DEBUGMODE
          if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
             && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
             )  
          {
          std::cout << "\n\n"
                    << "\nPropagation step : " << _propagation_counter
                    << "\n------------------------------------"
                    << "\nAttempting to move the wormhead over a time interval of " << _deltatime << " ( " << _deltatime_pause << " / " << _deltatime_halt << " ) "
                    << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                    << "\n\n";
          }
#endif
          if (_deltatime < _deltatime_halt)
          {
            worldline.wormhead_moves_to_new_time_without_winding(worldline.time()-_deltatime);
            _deltatime_pause -= _deltatime;
            _halted = false;
          }
          else
          {
            worldline.wormhead_moves_to_new_time_without_winding(worldline.time()-_deltatime_halt);
            _deltatime_pause -= _deltatime_halt;
            _halted = true;
          }
        }
      }

#ifdef DEBUGMODE
      if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
         && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
         ) 
      {
      std::cout << "\n\n"
                << "\nWormhead moves to time " << worldline.time() << " ( " << (_halted ? "HALTED" : "NOT HALTED") << " ) "
                << "\n============================================================================"
                << "\n" << worldline
                << "\n\n"
                << "\nWormhead"
                << "\n--------"
                << worldline.wormhead()
                << "\n"
                << "\nWormtail"
                << "\n--------"
                << wormpair.second
                << "\n"
                << "\nEnded (Wormhead movement)!"
                << "\n--------------------------"
                << "\n\n";
      }
#endif

      if (_deltatime_pause < 0.)   // Measure Green's function
      {
        using boost::numeric::operators::operator-;
        std::vector<int> wormpair_displacement = alps::applications::numeric::iround(coordinate(worldline.site()) - coordinate(wormpair.second.site()));

        if (inhomogeneous())
          green[idx(alps::numeric::abs(wormpair_displacement))]  += (0.5 * wormpair_state_);
        else
          green[idx(periodicize(wormpair_displacement))]         += wormpair_state_;

        if (is_experiment && model().name() == "trapped boson Hubbard")
        {
          using boost::numeric::operators::operator*;
          using alps::numeric::sq;
          using alps::applications::numeric::inner_product;

          if (inhomogeneous())
            green_tof[idx(alps::numeric::abs(wormpair_displacement))] += 
              (0.5 * wormpair_state_ * std::cos(inner_product(phase,sq(position(worldline.site()))-sq(position(wormpair.second.site())))));
          else
            green_tof[idx(periodicize(wormpair_displacement))] +=
              (wormpair_state_ * std::cos(inner_product(phase,sq(position(worldline.site()))-sq(position(wormpair.second.site())))));
        }

#ifdef DEBUGMODE
        if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
           && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
           )
        {
        using alps::numeric::operator<<;
        std::cout << "\n\n"
                  << "\nMeasuring Green's function within offdiagonal process... "  
                  << "\n"
                  << "\nInformation: "
                  << "\n"
                  << "\nWormhead site:         " << worldline.site() << "\t( " << coordinate(worldline.site()) << " ) "
                  << "\nWormtail site:         " << wormpair.second.site() << "\t( " << coordinate(wormpair.second.site()) << " )"
                  << "\nWormpair displacement: " << wormpair_displacement 
                  << "\nDensity Matrix Index:  " << (inhomogeneous() ? idx(alps::numeric::abs(wormpair_displacement)) : idx(wormpair_displacement % lattice().extent()))
                  << "\n\n";
        }
#endif
      }

      if (!_halted)
      {
        // Step 2A: Wormhead either inserts vertex or rebounces
        Site _site = worldline.site();
        insert_jump_or_bounce(_onsite_energy_relative); 
        if (!inhomogeneous())
        {
          using alps::numeric::operator+=;
          using boost::numeric::operators::operator-;
          using alps::applications::numeric::iround;
          winding_number += winding(iround(coordinate(worldline.site()) - coordinate(_site)));
        }
        return true;
      }

      else
      {
        LineIterator _vertexlineit = (worldline.forward() ? worldline.forwardlineit() : worldline.backwardlineit());

        // Step 2B: Wormhead touches wormtail (offdiagonal process ends)
        if (_vertexlineit->wormtail())
           return false;

        // Step 2C: Wormhead crosses kink
        if (  ( worldline.creation() &&  _vertexlineit->creation())
           || (!worldline.creation() && !_vertexlineit->creation())
           )
        {
          worldline.wormhead_crosses_vertex();
#ifdef DEBUGMODE
          if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
             && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
             )
          {
          std::cout << "\n\n"
                    << "\nWormhead crosses kink"
                    << "\n============================================================================"
                    << "\n" << worldline
                    << "\n\n"
                    << "\nWormhead"
                    << "\n--------"
                    << worldline.wormhead()
                    << "\n"
                    << "\nWormtail"
                    << "\n--------"
                    << wormpair.second
                    << "\n"
                    << "\nEnded (Wormhead crosses kink)!"
                    << "\n--------------------------"
                    << "\n\n";
          }
#endif
          return true;
        }
        // Step 2D: Wormhead either deletes/relinks vertex or rebounces
        else
        {
          Site _site = worldline.site();
          delete_relink_jump_or_bounce();
          if (!inhomogeneous())
          {
            using alps::numeric::operator+=;
            using boost::numeric::operators::operator-;
            using alps::applications::numeric::iround;
            winding_number += winding(iround(coordinate(worldline.site()) - coordinate(_site)));
          }
          return true;
        }
      }
    }

void
  directed_worm_algorithm
    ::insert_jump_or_bounce (double const onsite_energy_relative_)
    {
      std::vector<LineIterator> _neighborlineitsprev;  
      _neighborlineitsprev.reserve(num_neighbors(worldline.site()));

      std::vector<double>  _weights;
      _weights.reserve(num_neighbors(worldline.site()));

#ifdef DEBUGMODE
      if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
         && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
         )
      std::cout << "\n\nAttempting to insert into 1 of these neighbors : ";
#endif

      for (neighbor_bond_iterator _neighborbondit = neighbor_bonds(worldline.site()).first; _neighborbondit != neighbor_bonds(worldline.site()).second; ++_neighborbondit)
      {
#ifdef DEBUGMODE
        if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
           && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
           )
        std::cout << "\t" << target(*_neighborbondit);
#endif
        LineIterator  _neighborlineitprev  = --(worldline.lineit(target(*_neighborbondit), worldline.time())); 
        _neighborlineitsprev.push_back(_neighborlineitprev);                                              

        State _targetstate = _neighborlineitprev->state_after();
        bool  _increasing = alps::applications::increasing(worldline.forward(), worldline.creation());
        if (  ( _increasing && _targetstate == state_maximum[site_site_type[target(*_neighborbondit)]])
           || (!_increasing && _targetstate == state_minimum[site_site_type[target(*_neighborbondit)]])
           )
           _weights.push_back(0);
        else
           _weights.push_back(hopping_energy(*_neighborbondit, _targetstate, _increasing));
      }
#ifdef DEBUGMODE
      if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
         && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
         )
      std::cout << "\n\n";
#endif

      double _total_weight = std::accumulate(_weights.begin(), _weights.end(), onsite_energy_relative_);
      double _weight       = random()*_total_weight;
#ifdef DEBUGMODE
      if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
         && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
         )
      {
      std::cout << "\n\nTransition weights:\t" << onsite_energy_relative_ << "\t";
      std::copy(_weights.begin(), _weights.end(), std::ostream_iterator<double>(std::cout,"\t"));
      std::cout << "Total (selected) =\t" << _total_weight << " ( " << _weight << " ) ";
      std::cout << "\n------------------------------------------------" << "\n\n";
      }
#endif

      if (_weight < onsite_energy_relative_) 
      {
        worldline.wormhead_changes_its_direction();
        ++_propagation_failure_counter;
#ifdef DEBUGMODE
        if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
           && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
           )
        {
        std::cout << "\n\n"
                  << "\nWormhead bounces\n"
                  << "\n================"
                  << "\n" << worldline
                  << "\n\n"
                  << "\nWormhead"
                  << "\n--------"
                  << worldline.wormhead()
                  << "\n"
                  << "\nWormtail"
                  << "\n--------"
                  << wormpair.second
                  << "\n"
                  << "\nEnded (Wormhead bounces)!)"
                  << "\n--------------------------"
                  << "\n\n";
        }
#endif           
        return;    
      }
      _weight -= onsite_energy_relative_;

      for (neighbors_size_type _which_neighbor=0; _which_neighbor < num_neighbors(worldline.site()); ++_which_neighbor) 
      {
        if (_weight < _weights.back())
        {
          worldline.wormhead_jumps_to_new_site_and_inserts_vertex(_neighborlineitsprev.back());
#ifdef DEBUGMODE
          if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
             && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
             )
          {
          std::cout << "\n\n"
                    << "\nInserting an interaction to site " << _neighborlineitsprev.back()->site()
                    << "\n====================================================================="
                    << "\n" << worldline
                    << "\nWormhead"
                    << "\n--------"
                    << worldline.wormhead()
                    << "\n"
                    << "\nWormtail"
                    << "\n--------"
                    << wormpair.second
                    << "\n"
                    << "\nEnded (Inserting an interaction)!)"
                    << "\n----------------------------------"
                    << "\n\n";
          }
#endif
          return ;
        } 
        _weight -= _weights.back();
        _neighborlineitsprev .pop_back();
        _weights             .pop_back();
      }
    }

void
  directed_worm_algorithm
    ::delete_relink_jump_or_bounce()
    {
      LineIterator  _vertexlineit = worldline.forward() ? worldline.forwardlineit() : worldline.backwardlineit();
      LineIterator  _sourcelineit = worldline.lineit(_vertexlineit->partnersite(), _vertexlineit->time());
      bool          _forward      = !worldline.forward();    // Note: Here, we consider the direction of the reverse move.

      std::vector<LineIterator> _neighborlineitsprev;  
      _neighborlineitsprev.reserve(num_neighbors(_sourcelineit->site()));

      std::vector<double>  _weights;              
      _weights.reserve(num_neighbors(_sourcelineit->site()));
#ifdef DEBUGMODE
      if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
         && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
         )
      std::cout << "\n\nAttempting to delete/relink to 1 of these neighbors : ";
#endif
      for (neighbor_bond_iterator _neighborbondit = neighbor_bonds(_sourcelineit->site()).first; _neighborbondit != neighbor_bonds(_sourcelineit->site()).second; ++_neighborbondit)
      {
#ifdef DEBUGMODE
        if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
           && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
           )
        std::cout << "\t" << target(*_neighborbondit);
#endif
        LineIterator _neighborlineitprev  = --worldline.lineit(target(*_neighborbondit), _sourcelineit->time());
        _neighborlineitsprev.push_back(_neighborlineitprev);

        State  _targetstate = _neighborlineitprev->state_after();
        if (target(*_neighborbondit) == worldline.site())
          _targetstate = (_forward ? worldline.state_after() : worldline.state());
#ifdef DEBUGMODE
        if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
           && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
           )
        std::cout << "(" << _targetstate << ")";
#endif
        bool  _increasing = alps::applications::increasing(_forward, _sourcelineit->creation());
        if (  ( _increasing && _targetstate == state_maximum[site_site_type[target(*_neighborbondit)]])
           || (!_increasing && _targetstate == state_minimum[site_site_type[target(*_neighborbondit)]])
           )
           _weights.push_back(0);
        else
           _weights.push_back(hopping_energy(*_neighborbondit, _targetstate, _increasing));
      }
#ifdef DEBUGMODE
      if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
         && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
         )
      std::cout << "\n\n";
#endif

      double _onsite_energy_relative = onsite_energy_relative(*_sourcelineit, _forward);
      double _total_weight = std::accumulate(_weights.begin(), _weights.end(), _onsite_energy_relative);
      double _weight       = random() * _total_weight;
#ifdef DEBUGMODE
      if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
         && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
         )
      {
      std::cout << "\n\nTransition weights:\t" << _onsite_energy_relative << "\t";
      std::copy(_weights.begin(), _weights.end(), std::ostream_iterator<double>(std::cout,"\t"));
      std::cout << "Total (selected) =\t" << _total_weight << " ( " << _weight << " ) ";
      std::cout << "\n------------------------------------------------" << "\n\n";
      }
#endif

      if (_weight < _onsite_energy_relative)
      {
        worldline.wormhead_deletes_vertex_and_jumps_to_new_site(_sourcelineit);
#ifdef DEBUGMODE
        if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
           && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
           )
        {
        std::cout << "\n\n"
                  << "\nWormhead deletes vertex"
                  << "\n======================="
                  << "\n" << worldline
                  << "\nWormhead"
                  << "\n--------"
                  << worldline.wormhead()
                  << "\n"
                  << "\nWormtail"
                  << "\n--------"
                  << wormpair.second
                  << "\n"
                  << "\nEnded (Deleting a vertex)!)"
                  << "\n----------------------------------"
                  << "\n\n";
        }
#endif
        return;
      }
      _weight -= _onsite_energy_relative;

      for (neighbors_size_type _which_neighbor=0; _which_neighbor < num_neighbors(_sourcelineit->site()); ++_which_neighbor)
      {
        if (_weight < _weights.back())
        {
          LineIterator  _targetlineitprev  = _neighborlineitsprev.back();
          if (_targetlineitprev->site() == worldline.site())
          {
            worldline.wormhead_changes_its_direction();
            ++_propagation_failure_counter;
#ifdef DEBUGMODE
            if (  _sweep_counter       >= DEBUGMODE_START_COUNTER              && _sweep_counter       <= DEBUGMODE_END_COUNTER
               && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
               )  
            {
            std::cout << "\n\n"
                      << "\nWormhead rebounces" 
                      << "\n=================="
                      << "\n" << worldline
                      << "\nWormhead"
                      << "\n--------"
                      << worldline.wormhead()
                      << "\n"
                      << "\nWormtail"
                      << "\n--------"
                      << wormpair.second
                      << "\n"                    
                      << "\nEnded (Wormhead rebounces)!)"
                      << "\n----------------------------"                    
                      << "\n\n";
            }
#endif
            return;
          }
          worldline.wormhead_relinks_vertex_and_jumps_to_new_site(_sourcelineit, _targetlineitprev);  
#ifdef DEBUGMODE
          if (  _sweep_counter       >= DEBUGMODE_START_COUNTER   && _sweep_counter       <= DEBUGMODE_END_COUNTER
             && _propagation_counter >= DEBUGMODE_START_PROPAGATION_COUNTER  && _propagation_counter <= DEBUGMODE_END_PROPAGATION_COUNTER
             )
          {          
          std::cout << "\n\n"
                    << "\nWormhead relinks vertex to site " << _targetlineitprev->site() << " (which_neighbor ? " << _which_neighbor << " )"              
                    << "\n=================================================================================================="
                    << "\n" << worldline
                    << "\nWormhead"
                    << "\n--------"
                    << worldline.wormhead()
                    << "\n"
                    << "\nWormtail"
                    << "\n--------"
                    << wormpair.second
                    << "\n"
                    << "\nEnded (Wormhead relink)!)"
                    << "\n-------------------------"
                    << "\n\n";
          }
#endif
          return;
        }
        _weight -= _weights.back();
        _neighborlineitsprev .pop_back();
        _weights             .pop_back();
      }
    }

} // ending namespace applications
} // ending namespace alps


int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
  try 
  {
#endif
   return alps::scheduler::start(argc, argv, alps::scheduler::SimpleMCFactory<alps::applications::directed_worm_algorithm>());
#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (std::exception& exc) 
  {
    std::cerr << exc.what() << "\n";
    alps::comm_exit(true);
    return -1;
  }
  catch (...) 
  {
    std::cerr << "Fatal Error: Unknown Exception!\n";
    return -2;
  }
#endif
}
