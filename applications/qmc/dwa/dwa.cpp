/****************************************************************************
*
* ALPS Project Applications: Directed Worm Algorithm  
*
* ALPS Project Applications: Directed Worm Algorithm  
*
* Copyright (C) 2013 by Matthias Troyer  <troyer@phys.ethz.ch> ,
*                       Lode Pollet      <pollet@phys.ethz.ch> ,
*                       Ping Nang Ma     <pingnang@phys.ethz.ch> 
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

#include "dwa.hpp"
#include <alps/scheduler.h>
#include <alps/osiris/comm.h>


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
          << "\n* ALPS Project Applications: Directed Worm Algorithm"
          << "\n*"
          << "\n* Copyright (C) 2013 by Matthias Troyer  <troyer@phys.ethz.ch> ," 
          << "\n*                       Lode Pollet      <pollet@phys.ethz.ch> ,"
          << "\n*                       Ping Nang Ma     <pingnang@phys.ethz.ch>"
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
          << parms
          << "\n\n";
      out << "Simulation\n"
          << "==========\n\n"
          << "\t Sweeps                              : " << _total_sweeps                  << "\n"
          << "\t Thermalization sweeps               : " << _thermalization_sweeps         << "\n"
          << "\t Skip (Sweeps per measurement)       : " << _skip                          << "\n"
          << "\t Sweep counter                       : " << _sweep_counter                 << "\n"
          << "\t Sweep failure counter               : " << _sweep_failure_counter         << "\n"
          << "\t Propagation counter                 : " << _propagation_counter           << "\n"
          << "\t Propagation failure counter         : " << _propagation_failure_counter   << "\n"
          << "\n"
          << "\t Interested in time-of-flight images : " << finite_tof                     << "\n"
          << "\t Waist correction                    : " << finite_waist                   << "\n"
          << "\n\n";
      out << "Lattice\n"
          << "=======\n\n"
          << "Name              : " << parms["LATTICE"]                         << "\n"
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
          << "MEASURE                          : " << measure_                        << "\n"
          << "MEASURE[Simulation Speed]        : " << measure_simulation_speed_       << "\n"
          << "MEASURE[Total Particle Number^2] : " << measure_number2_                << "\n"
          << "MEASURE[Energy^2]                : " << measure_energy2_                << "\n"
          << "MEASURE[Density^2]               : " << measure_density2_               << "\n"
          << "MEASURE[Energy Density^2]        : " << measure_energy_density2_        << "\n"
          << "MEASURE[Winding Number]          : " << measure_winding_number_         << "\n"
          << "MEASURE[Local Kink:Number]       : " << measure_local_num_kinks_        << "\n"
          << "MEASURE[Local Density]           : " << measure_local_density_          << "\n" 
          << "MEASURE[Local Density^2]         : " << measure_local_density2_         << "\n"
          << "MEASURE[Green Function]          : " << measure_green_function_         << "\n"
          << "\n\n"; 
    }


directed_worm_algorithm

//  ::directed_worm_algorithm(alps::hdf5::archive & ar)    // Luka's NGS scheduler is not ready.
//    : qmcbase<>(ar)                                      // Luka's NGS scheduler is not ready.

  ::directed_worm_algorithm(const alps::ProcessList& plist_, const alps::Parameters& parms_, int n)
    : QMCRun<>(plist_, parms_, n)

    // regarding MC simulation
    , _sweep_counter               (0)
    , _sweep_failure_counter       (0)
    , _propagation_counter         (0)
    , _propagation_failure_counter (0)
    , _thermalization_sweeps       (parms_.value_or_default("THERMALIZATION",0))
    , _total_sweeps                (parms_.value_or_default("SWEEPS",10000000))
    , _skip                        (parms_.value_or_default("SKIP",1))
    // regarding lattice
    , is_periodic_ (std::find((this->lattice().boundary()).begin(), (this->lattice().boundary()).end(), "open") == (this->lattice().boundary()).end())
    // regarding experiment
    , finite_tof   (is_periodic_ ? false : (parms_.defined("tof_phase")))
    , finite_waist (parms_.defined("waist"))
    // regarding worldline
    , wl (num_sites())
    // regarding measurements
    , measure_                       (parms_.value_or_default("MEASURE",true))
    , measure_simulation_speed_      (parms_.value_or_default("MEASURE[Simulation Speed",true))
    , measure_number2_               (parms_.value_or_default("MEASURE[Total Particle Number^2]",false))
    , measure_energy2_               (parms_.value_or_default("MEASURE[Energy^2]",false))
    , measure_density2_              (inhomogeneous() ? false : static_cast<bool>(parms_.value_or_default("MEASURE[Density^2]",false)))
    , measure_energy_density2_       (inhomogeneous() ? false : static_cast<bool>(parms_.value_or_default("MEASURE[Energy Density^2]",false)))
    , measure_winding_number_        (!is_periodic_  ? false :  static_cast<bool>(parms_.value_or_default("MEASURE[Winding Number]",false)))
    , measure_local_num_kinks_       (parms_.value_or_default("MEASURE[Local Kink: Number]",false))
    , measure_local_density_         (parms_.value_or_default("MEASURE[Local Density]",false))
    , measure_local_density2_        (parms_.value_or_default("MEASURE[Local Density^2]",false))
    , measure_green_function_        (parms_.value_or_default("MEASURE[Green Function]",false))
    // regarding development
    , alps_dwa_development_model_parsing_mode  (parms_.value_or_default("DEVELOPMENT[MODEL PARSING]", -1))
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
    if (measure_simulation_speed_)
      std::time(&_simulation_timer.first);
  }


void
  directed_worm_algorithm
    ::save(alps::hdf5::archive & ar) const
    {
      std::cout << "\nSaving simulation ... starting ...\n\n";

      // Worldlines
      std::cout << "\t\ti. worldlines \t\t ... starting ... \n";
      //ar << alps::make_pvp("/simulation/worldlines/num_sites", num_sites()); 
      //ar << alps::make_pvp("/simulation/worldlines", wl);
      wl.save(ar);
      std::cout << "\t\t\t\t ... done.\n\n";

      // Parameters
      std::cout << "\t\tii. parameters \t\t ... starting ... \n";
      ar << alps::make_pvp("/parameters", parms);
      std::cout << "\t\t\t\t ... done.\n\n";

      // Counters
      ar << alps::make_pvp("/simulation/sweep_counter", _sweep_counter);

      // Simulation results
      std::cout << "\t\tiii. measurements \t\t ... starting ... \n";
      ar << alps::make_pvp("/simulation/results", measurements);
      std::cout << "\t\t\t\t ... done.\n\n";

      // Other useful quantities
      std::cout << "\t\tiv. other useful quantities \t\t ... starting ... \n";
      //ar << alps::make_pvp("/simulation/positions", position_lookup);
      //if (measure_green_function_)
      //  ar << alps::make_pvp("/simulation/green_coordinates", lattice().distance_labels());
      std::cout << "\t\t\t\t ... done.\n\n";

      std::cout << "\t\t ... done.\n\n";
    }


void 
  directed_worm_algorithm
    ::load(alps::hdf5::archive & ar)  
    {  
      // Worldlines
      if (ar.is_group("/simulation/worldlines")) {
        wl.load(ar); 
        std::cout << "\nWorldlines configuration loaded from input file.\n";
        std::cout << "\t\t- valid ? " << wl.is_valid(state_maximum[0]) << "\n";
      }
      else
        std::cout << "\nNo worldlines configuration found.\n";

      // Counters
      if (ar.is_group("/simulation/sweep_counter")) {
        ar >> alps::make_pvp("/simulation/sweep_counter", _sweep_counter);
      }

      // Measurements
      if (ar.is_group("/simulation/results")) {
        ar >> alps::make_pvp("/simulation/results", measurements);
        std::cout << "\nPrevious measurements loaded from input file.\n";
      }
      else
        std::cout << "\nNo measurements found.\n";
    }


// Tama Ma is here...
void
directed_worm_algorithm
::initialize_onsite_hamiltonian()
{
  for (site_iterator it = sites().first; it != sites().second; ++it) {
    unsigned this_site_type = site_type(*it);
    unsigned this_inhomogeneous_site_type = inhomogeneous_site_type(*it);
    if (this_inhomogeneous_site_type >= onsite_matrix.size()) 
    {
      alps::Parameters parms2(parms);
      if (inhomogeneous_sites()) {
        throw_if_xyz_defined(parms,*it);   // check whether x, y, or z is set
        parms2 << coordinate_as_parameter(*it); // set x, y and z
      }

      boost::multi_array<double,2> this_site_matrix =
        alps::get_matrix(double(), model().site_term(this_site_type),
             model().basis().site_basis(this_site_type), parms2); 

      std::vector<double> this_site_matrix_diagonal(this_site_matrix.shape()[0]);
      for (unsigned idx=0; idx < this_site_matrix_diagonal.size(); ++idx)
         this_site_matrix_diagonal[idx] = this_site_matrix[idx][idx] * beta;

      onsite_matrix.push_back(this_site_matrix_diagonal);  // Note: It HAS BEEN multiplied by beta here... 
    }
  }
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

      if (alps_dwa_development_model_parsing_mode == 0)
      {
        initialize_onsite_hamiltonian();
      }
      else if (alps_dwa_development_model_parsing_mode == -1)
      {
        //// tamama's fast evaluate scheme (else evaluating x, y, z for say 100^3 lattice will take forever)
        //// Step 1: partial evaluate everything up to x,y,z terms
        std::cout << "\t\t\t\t\ta. Step 1: partial evaluate everything up to x,y,z terms \t... starting\n";
        std::vector<alps::SiteTermDescriptor> _siteterms = model().site_terms();
  
        std::vector<std::vector<double> > _bare_values;
        std::vector<std::vector<alps::expression::Expression<double> > > _bare_x_expressions;
        std::vector<std::vector<alps::expression::Expression<double> > > _bare_y_expressions;
        std::vector<std::vector<alps::expression::Expression<double> > > _bare_z_expressions;
        std::vector<std::vector<alps::expression::Expression<double> > > _bare_expressions;  // leftover
        bool _can_evaluate_individual_bare_expressions = true;
        bool _can_evaluate_bare_expressions = true;
        for (unsigned int i=0; i < maximum_sitetype+1; ++i) {
          alps::expression::Expression<double> _bare_expression(_siteterms[i].term());
          _bare_expression.flatten();
          _bare_expression.simplify();
  
          //std::cout << "\nsite_type = " << i << "\n";
          //std::cout << "\n_bare_expression.output(std::cout):\n";
          //_bare_expression.output(std::cout);
          //std::cout << "\n";
          //std::cin.get(); 
  
          _bare_values.push_back(std::vector<double>());
          _bare_x_expressions.push_back(std::vector<alps::expression::Expression<double> >());
          _bare_y_expressions.push_back(std::vector<alps::expression::Expression<double> >());
          _bare_z_expressions.push_back(std::vector<alps::expression::Expression<double> >());
          _bare_expressions.push_back(std::vector<alps::expression::Expression<double> >());
  
          alps::site_basis<short> _states(model().basis().site_basis());
          _bare_expressions.reserve(_states.size());
          _bare_values.reserve(_states.size());
          for (unsigned int m=0; m<_states.size(); ++m) {
            alps::expression::Expression<double> this_expression(_bare_expression);
            alps::expression::Expression<double> this_x_expression;
            alps::expression::Expression<double> this_y_expression;
            alps::expression::Expression<double> this_z_expression;
            alps::SiteOperatorEvaluator<short,double> _site_basis_evaluator(_states[m], model().basis().site_basis(), parms, _siteterms[i].site());
            this_expression.partial_evaluate(_site_basis_evaluator);   // evaluate all operator terms...
            this_expression.simplify();
            alps::expression::Term<double> this_zeroth_term = this_expression.zeroth_term();
            if (this_zeroth_term.can_evaluate())
            {
              _bare_values[i].push_back(this_zeroth_term.value());
              this_expression -= this_zeroth_term;
              this_expression.simplify();   // contains x,y,z terms only
            }
            else
              _bare_values[i].push_back(0.);
  
            if (!this_expression.can_evaluate())
              _can_evaluate_individual_bare_expressions = false;
  
            this_x_expression = this_expression.expression_dependent_only_on("x");
            _bare_x_expressions[i].push_back(this_x_expression);
            this_expression -= this_x_expression;
            this_expression.simplify();
  
            this_y_expression = this_expression.expression_dependent_only_on("y");
            _bare_y_expressions[i].push_back(this_y_expression);
            this_expression -= this_y_expression;
            this_expression.simplify();
  
            this_z_expression = this_expression.expression_dependent_only_on("z");
            _bare_z_expressions[i].push_back(this_z_expression);
            this_expression -= this_z_expression;
            this_expression.simplify();
  
            if (!this_expression.can_evaluate())
              _can_evaluate_bare_expressions = false;
  
            _bare_expressions[i].push_back(this_expression);
  
            //std::cout << "\nsite_type = " << i << "\tstate = " << m << "\n";
            //std::cout << "\nthis_expression.output(std::cout):\n";
            //this_expression.output(std::cout);
            //std::cout << "\n";
            //std::cin.get(); 
          }        
        }
  
        //std::cout << "\n";
        //for (unsigned int i=0; i<_bare_values.size(); ++i) {
        //  using alps::numeric::operator<<;
        //  std::cout << "_bare_values[" << i << "] : " << _bare_values[i] << "\n";
        //}
        //std::cin.get();
  
        //std::cout << "\n";
        //for (unsigned int i=0; i<_bare_expressions.size(); ++i) {
        //for (unsigned int m=0; m<_bare_expressions[i].size(); ++m) {
        //  std::cout << "_bare_expressions[" << i << "][" << m << "] : ";
        //  std::cout << "[x] : ";
        //  _bare_x_expressions[i][m].output(std::cout);
        //  std::cout << ", [y] : ";
        //  _bare_y_expressions[i][m].output(std::cout);
        //  std::cout << ", [z] : ";
        //  _bare_z_expressions[i][m].output(std::cout);
        //  std::cout << ", [others]: ";
        //  _bare_expressions[i][m].output(std::cout);
        //  std::cout << "\n";
        //}
        //}
        //std::cin.get();
        std::cout << "\t\t\t\t\t\t\t ... done.\n";
  
        //// Step 2: Initialize onsite matrix 
        std::cout << "\t\t\t\t\tb. Step 2: Initialize onsite matrix \t ... starting:\n";
  
        ////// _bare_values -> onsite_matrix
        for (site_iterator it = sites().first; it != sites().second; ++it) {
          unsigned int this_site_type = inhomogeneous_site_type(*it);
          if (this_site_type >= onsite_matrix.size()) 
            onsite_matrix.push_back(_bare_values[site_type(*it)]);  // Note: It has not been yet multiplied by beta here... 
        }
        //std::cout << "\n";
        //for (unsigned int i=0; i<onsite_matrix.size(); ++i) {
        //  using alps::numeric::operator<<;
        //  std::cout << "onsite_matrix[" << i << "] : " << onsite_matrix[i] << "\n";
        //}
        //std::cin.get();
  
        ////// _bare_???_expressions, ... -> onsite_matrix
        if (!_can_evaluate_individual_bare_expressions && inhomogeneous_sites())  // it only makes sense that x,y,z are for inhomogeneous lattice...
        {
          if (dimension() >= 1)   // _bare_x_expression -> onsite_matrix
          {
            std::map<double, unsigned int> coordinate_type;
            for (site_iterator it = sites().first; it != sites().second; ++it) 
              coordinate_type.insert(std::make_pair(coordinate(*it)[0],coordinate_type.size()));
            
            std::vector<std::vector<std::vector<double> > > _bare_directional_values;     // dim1: coordinate_type, dim2: site_type, dim3: state
            _bare_directional_values.resize(coordinate_type.size());
            for (std::map<double, unsigned int>::iterator mit = coordinate_type.begin(); mit != coordinate_type.end(); ++mit) {
            for (unsigned int i=0; i<_bare_x_expressions.size(); ++i) {
            _bare_directional_values[mit->second].push_back(std::vector<double>());
            _bare_directional_values[mit->second][i].reserve(_bare_x_expressions[i].size());
            for (unsigned int m=0; m<_bare_x_expressions[i].size(); ++m) {
               alps::expression::Expression<double> this_expression = _bare_x_expressions[i][m];
               alps::Parameters this_parms;
               this_parms["x"] = boost::lexical_cast<std::string>(mit->first);
               alps::expression::ParameterEvaluator<double> this_evaluator(this_parms);
               this_expression.partial_evaluate(this_evaluator);
               this_expression.simplify();
               _bare_directional_values[mit->second][i].push_back(this_expression.value());
            }
            }
            }
            //for (std::map<double, unsigned int>::iterator mit = coordinate_type.begin(); mit != coordinate_type.end(); ++mit) {
            //for (unsigned int i=0; i<_bare_x_expressions.size(); ++i) {
            //for (unsigned int m=0; m<_bare_x_expressions[i].size(); ++m) {
            //  std::cout << "x : " << mit->first << " , _bare_x_expressions[" << i << "][" << m << "] : ";
            //  _bare_x_expressions[i][m].output(std::cout);
            //  std::cout << " , _bare_directional_value : " << _bare_directional_values[mit->second][i][m] << "\n";
            //}
            //}
            //std::cin.get();
            //}
            //std::cin.get();
  
            std::vector<std::vector<double> > onsite_matrix_plus;
            for (site_iterator it = sites().first; it != sites().second; ++it)
            {
              unsigned int this_site_type = inhomogeneous_site_type(*it);
              if (this_site_type >= onsite_matrix_plus.size())  
                onsite_matrix_plus.push_back(_bare_directional_values[coordinate_type[coordinate(*it)[0]]][site_type(*it)]);
            }
            //for (unsigned int i=0; i<onsite_matrix.size(); ++i) {
            //  using alps::numeric::operator<<;
            //  std::cout << "onsite_matrix     [" << i << "] : " << onsite_matrix[i] << "\n";
            //  std::cout << "onsite_matrix_plus[" << i << "] : " << onsite_matrix_plus[i] << "\n";
            //}
            //std::cin.get();
  
            for (unsigned int i=0; i<onsite_matrix.size(); ++i) {
              using boost::numeric::operators::operator+;
              onsite_matrix[i] = onsite_matrix[i] + onsite_matrix_plus[i];
            }
            //for (unsigned int i=0; i<onsite_matrix.size(); ++i) {
            //  using alps::numeric::operator<<;
            //  std::cout << "onsite_matrix     [" << i << "] : " << onsite_matrix[i] << "\n";
            //}
            //std::cin.get();
          }
  
          if (dimension() >= 2)   // _bare_y_expression -> onsite_matrix
          {
            std::map<double, unsigned int> coordinate_type;
            for (site_iterator it = sites().first; it != sites().second; ++it) 
              coordinate_type.insert(std::make_pair(coordinate(*it)[1],coordinate_type.size()));
            
            std::vector<std::vector<std::vector<double> > > _bare_directional_values;     // dim1: coordinate_type, dim2: site_type, dim3: state
            _bare_directional_values.resize(coordinate_type.size());
            for (std::map<double, unsigned int>::iterator mit = coordinate_type.begin(); mit != coordinate_type.end(); ++mit) {
            for (unsigned int i=0; i<_bare_y_expressions.size(); ++i) {
            _bare_directional_values[mit->second].push_back(std::vector<double>());
            _bare_directional_values[mit->second][i].reserve(_bare_y_expressions[i].size());
            for (unsigned int m=0; m<_bare_y_expressions[i].size(); ++m) {
               alps::expression::Expression<double> this_expression = _bare_y_expressions[i][m];
               alps::Parameters this_parms;
               this_parms["y"] = boost::lexical_cast<std::string>(mit->first);
               alps::expression::ParameterEvaluator<double> this_evaluator(this_parms);
               this_expression.partial_evaluate(this_evaluator);
               this_expression.simplify();
               _bare_directional_values[mit->second][i].push_back(this_expression.value());
            }
            }
            }
            //for (std::map<double, unsigned int>::iterator mit = coordinate_type.begin(); mit != coordinate_type.end(); ++mit) {
            //for (unsigned int i=0; i<_bare_y_expressions.size(); ++i) {
            //for (unsigned int m=0; m<_bare_y_expressions[i].size(); ++m) {
            //  std::cout << "y : " << mit->first << " , _bare_y_expressions[" << i << "][" << m << "] : ";
            //  _bare_y_expressions[i][m].output(std::cout);
            //  std::cout << " , _bare_directional_value : " << _bare_directional_values[mit->second][i][m] << "\n";
            //}
            //}
            //std::cin.get();
            //}
            //std::cin.get();
  
            std::vector<std::vector<double> > onsite_matrix_plus;
            for (site_iterator it = sites().first; it != sites().second; ++it)
            {
              unsigned int this_site_type = inhomogeneous_site_type(*it);
              if (this_site_type >= onsite_matrix_plus.size())  
                onsite_matrix_plus.push_back(_bare_directional_values[coordinate_type[coordinate(*it)[1]]][site_type(*it)]);
            }
            //for (unsigned int i=0; i<onsite_matrix.size(); ++i) {
            //  using alps::numeric::operator<<;
            //  std::cout << "onsite_matrix     [" << i << "] : " << onsite_matrix[i] << "\n";
            //  std::cout << "onsite_matrix_plus[" << i << "] : " << onsite_matrix_plus[i] << "\n";
            //}
            //std::cin.get();
  
            for (unsigned int i=0; i<onsite_matrix.size(); ++i) {
              using boost::numeric::operators::operator+;
              onsite_matrix[i] = onsite_matrix[i] + onsite_matrix_plus[i];
            }
            //for (unsigned int i=0; i<onsite_matrix.size(); ++i) {
            //  using alps::numeric::operator<<;
            //  std::cout << "onsite_matrix     [" << i << "] : " << onsite_matrix[i] << "\n";
            //}
            //std::cin.get();
          }
  
          if (dimension() >= 3)   // _bare_z_expression -> onsite_matrix
          {
            std::map<double, unsigned int> coordinate_type;
            for (site_iterator it = sites().first; it != sites().second; ++it) 
              coordinate_type.insert(std::make_pair(coordinate(*it)[2],coordinate_type.size()));
            
            std::vector<std::vector<std::vector<double> > > _bare_directional_values;     // dim1: coordinate_type, dim2: site_type, dim3: state
            _bare_directional_values.resize(coordinate_type.size());
            for (std::map<double, unsigned int>::iterator mit = coordinate_type.begin(); mit != coordinate_type.end(); ++mit) {
            for (unsigned int i=0; i<_bare_z_expressions.size(); ++i) {
            _bare_directional_values[mit->second].push_back(std::vector<double>());
            _bare_directional_values[mit->second][i].reserve(_bare_z_expressions[i].size());
            for (unsigned int m=0; m<_bare_z_expressions[i].size(); ++m) {
               alps::expression::Expression<double> this_expression = _bare_z_expressions[i][m];
               alps::Parameters this_parms;
               this_parms["z"] = boost::lexical_cast<std::string>(mit->first);
               alps::expression::ParameterEvaluator<double> this_evaluator(this_parms);
               this_expression.partial_evaluate(this_evaluator);
               this_expression.simplify();
               _bare_directional_values[mit->second][i].push_back(this_expression.value());
            }
            }
            }
            //for (std::map<double, unsigned int>::iterator mit = coordinate_type.begin(); mit != coordinate_type.end(); ++mit) {
            //for (unsigned int i=0; i<_bare_z_expressions.size(); ++i) {
            //for (unsigned int m=0; m<_bare_z_expressions[i].size(); ++m) {
            //  std::cout << "z : " << mit->first << " , _bare_z_expressions[" << i << "][" << m << "] : ";
            //  _bare_z_expressions[i][m].output(std::cout);
            //  std::cout << " , _bare_directional_value : " << _bare_directional_values[mit->second][i][m] << "\n";
            //}
            //}
            //std::cin.get();
            //}
            //std::cin.get();
  
            std::vector<std::vector<double> > onsite_matrix_plus;
            for (site_iterator it = sites().first; it != sites().second; ++it)
            {
              unsigned int this_site_type = inhomogeneous_site_type(*it);
              if (this_site_type >= onsite_matrix_plus.size())  
                onsite_matrix_plus.push_back(_bare_directional_values[coordinate_type[coordinate(*it)[2]]][site_type(*it)]);
            }
            //for (unsigned int i=0; i<onsite_matrix.size(); ++i) {
            //  using alps::numeric::operator<<;
            //  std::cout << "onsite_matrix     [" << i << "] : " << onsite_matrix[i] << "\n";
            //  std::cout << "onsite_matrix_plus[" << i << "] : " << onsite_matrix_plus[i] << "\n";
            //}
            //std::cin.get();
  
            for (unsigned int i=0; i<onsite_matrix.size(); ++i) {
              using boost::numeric::operators::operator+;
              onsite_matrix[i] = onsite_matrix[i] + onsite_matrix_plus[i];
            }
            //for (unsigned int i=0; i<onsite_matrix.size(); ++i) {
            //  using alps::numeric::operator<<;
            //  std::cout << "onsite_matrix     [" << i << "] : " << onsite_matrix[i] << "\n";
            //}
            //std::cin.get();
          }
  
          if (!_can_evaluate_bare_expressions)  // other expressions that consists of mixed x,y,z terms... 
          {
            std::vector<std::vector<double> > onsite_matrix_plus;
            for (site_iterator it = sites().first; it != sites().second; ++it)
            {    
              unsigned int this_site_type = inhomogeneous_site_type(*it);
              if (this_site_type >= onsite_matrix_plus.size())  
              {
                alps::Parameters this_parms;              
                this_parms << coordinate_as_parameter(this_site_type); // set x, y and z
  
                std::vector<double> this_values_plus;
                for (unsigned int m=0; m<_bare_expressions[site_type(*it)].size(); ++m) {
                  alps::expression::Expression<double> this_expression = _bare_expressions[site_type(*it)][m];
                  alps::expression::ParameterEvaluator<double> this_evaluator(this_parms);
                  this_expression.partial_evaluate(this_evaluator);
                  this_expression.simplify();
                  this_values_plus.push_back(this_expression.value()); 
                }
  
                onsite_matrix_plus.push_back(this_values_plus);
              }
            }    
            //for (unsigned int i=0; i<onsite_matrix.size(); ++i) {
            //  using alps::numeric::operator<<;
            //  std::cout << "onsite_matrix     [" << i << "] : " << onsite_matrix[i] << "\n";
            //  std::cout << "onsite_matrix_plus[" << i << "] : " << onsite_matrix_plus[i] << "\n";
            //}
            //std::cin.get();
  
            for (unsigned int i=0; i<onsite_matrix.size(); ++i) {
              using boost::numeric::operators::operator+;
              onsite_matrix[i] = onsite_matrix[i] + onsite_matrix_plus[i];
            }
            //for (unsigned int i=0; i<onsite_matrix.size(); ++i) {
            //  using alps::numeric::operator<<;
            //  std::cout << "onsite_matrix     [" << i << "] : " << onsite_matrix[i] << "\n";
            //}
            //std::cin.get();
          }
        }
  
        for (unsigned int i=0; i<onsite_matrix.size(); ++i) {
          using boost::numeric::operators::operator*;
          onsite_matrix[i] = onsite_matrix[i] * beta;
        }

      }  // ending if (alps_dwa_development_model_parsing_mode == -1)

      std::cout << "\t\t\t\t\t\t\t ... done.\n";

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

      // modify bond strength matrix if tx/t, ty/t, tz/t exists
      if (parms.defined("tx_t") || parms.defined("ty_t") || parms.defined("tz_t"))
      {
        std::vector<double> t_relative(dimension(), 1.);
        if (parms.defined("tx_t"))  t_relative[0] = parms.value_or_default("tx_t", 1.); 
        if (parms.defined("ty_t"))  t_relative[1] = parms.value_or_default("ty_t", 1.);
        if (parms.defined("tz_t"))  t_relative[2] = parms.value_or_default("tz_t", 1.);

        //using alps::numeric::operator<<;
        //std::cout << "\nt (relative) : " << t_relative << "\n";
        //std::cin.get(); 
     
        for (bond_iterator it=bonds().first; it!=bonds().second; ++it)
        {
          for (int i=0; i<dimension(); ++i)
            if (bond_vector(*it)[i] != 0.)
              bond_strength_matrix[index(*it)] *= t_relative[i];
        }
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

      // waist correction
      if (finite_waist)
      {
        double dummy_V0 = static_cast<double>(parms["V0"]);
        double dummy_w0 = static_cast<double>(parms["waist"]);
        std::vector<double> delta_VT(num_sites(), dummy_V0);
        for (unsigned int i=0; i<num_sites(); ++i) {
          double dummy = 2. * std::inner_product(position_lookup[i].begin(), position_lookup[i].end(), position_lookup[i].begin(), 0.) / (dummy_w0*dummy_w0);
          delta_VT[i] -= dummy_V0 * (dummy + std::exp(-dummy)); 
          delta_VT[i] *= beta;
        }

        for (unsigned int i=0; i<num_sites(); ++i) {
          using boost::numeric::operators::operator+;
          using boost::numeric::operators::operator*;
          onsite_matrix[i] = onsite_matrix[i] + delta_VT[i] * (diagonal_matrix_element.begin()->second[site_type(i)]);
        }
      }

      // component coordinate lookup table 
      std::vector<std::set<int> > component_coordinate_x100_set(dimension());
      for (unsigned int site=0; site<num_sites(); ++site)
        for (int i=0; i<dimension(); ++i)
          component_coordinate_x100_set[i].insert(boost::math::iround(100*coordinate(site)[i]));  
      component_coordinate_x100_lookup.resize(dimension());
      for (int i=0; i<dimension(); ++i)
        component_coordinate_x100_lookup[i] = std::vector<int>(component_coordinate_x100_set[i].begin(), component_coordinate_x100_set[i].end());

      // site lookup table
      int site_lookup_size = 1;
      for (int i=0; i<dimension(); ++i)
        site_lookup_size *= component_coordinate_x100_lookup[i].size();
      site_lookup.resize(site_lookup_size, -1);

      for (int site=0; site<num_sites(); ++site)
      {
        std::vector<int> site_idx(dimension()); 
        for (unsigned int i=0; i<dimension(); ++i)
          site_idx[i] = std::distance(component_coordinate_x100_lookup[i].begin(), std::lower_bound(component_coordinate_x100_lookup[i].begin(), component_coordinate_x100_lookup[i].end(), boost::math::iround(100*coordinate(site)[i])));
        site_lookup[componentIndex2Index(site_idx)] = site;
      }

      // phase and phase lookup table
      if (finite_tof)
      {
        std::vector<double> inverse_tof;
        std::string tof_phase = boost::lexical_cast<std::string>(parms["tof_phase"]);
        tof_phase.erase(std::remove(tof_phase.begin(), tof_phase.end(), '['), tof_phase.end());
        tof_phase.erase(std::remove(tof_phase.begin(), tof_phase.end(), ']'), tof_phase.end());
        tof_phase.erase(std::remove(tof_phase.begin(), tof_phase.end(), ' '), tof_phase.end());

        if (tof_phase.find_first_of(',') == std::string::npos)
          inverse_tof = std::vector<double>(dimension(), boost::lexical_cast<double>(tof_phase));
        else {
          while(tof_phase.find_first_of(',') != std::string::npos) {
            inverse_tof.push_back(boost::lexical_cast<double>(std::string(tof_phase, 0, tof_phase.find_first_of(','))));
            tof_phase.erase(0, tof_phase.find_first_of(',')+1);
          }
          inverse_tof.push_back(boost::lexical_cast<double>(tof_phase));
          inverse_tof.resize(dimension());
        } 
          
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

      // phase lookup table
      out << "Position / Phase Lookup Table:\n"
          << "----------------------\n\n";
      for (unsigned int site=0; site<num_sites(); ++site)  
        out << "Site : " << site << " , Site position : " << position(site) << " , Phase : " << phase(site) << "\n";
      out << "\n";

      // component coordinate lookup table
      out << "Component coordinate (x100) Lookup Table:\n"
          << "----------------------\n\n";
      for (unsigned int i=0; i<dimension(); ++i)
      {
        out << "dim = " << i << " : ";
        for (unsigned int j=0; j<component_coordinate_x100_lookup[i].size(); ++j)
          out << component_coordinate_x100_lookup[i][j] << "  ";
        out << "\n";
      }
      out << "\n\n";
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
        << alps::RealObservable   ("Total Particle Number")
        << alps::RealObservable   ("Energy")
        << alps::RealObservable   ("Energy:Vertex")
        << alps::RealObservable   ("Energy:Onsite")
      ;

      if (!inhomogeneous()) 
      {
        measurements 
          << alps::RealObservable ("Density")
          << alps::RealObservable ("Energy Density")
          << alps::RealObservable ("Energy Density:Vertex")
          << alps::RealObservable ("Energy Density:Onsite")
        ;
      }

      if (measure_number2_)
        measurements << alps::RealObservable ("Total Particle Number^2");
      if (measure_energy2_)
        measurements << alps::RealObservable ("Energy^2");
      if (measure_density2_)
        measurements << alps::RealObservable ("Density^2");
      if (measure_energy_density2_)
        measurements << alps::RealObservable ("Energy Density^2");

      if (measure_winding_number_) 
      {
        measurements << alps::RealVectorObservable ("Winding Number^2");
        measurements << alps::RealObservable       ("Stiffness");
      }

      if (measure_local_num_kinks_) 
      {
        _num_kinks_cache.resize(wl.num_sites(),0.);
        measurements << alps::SimpleRealVectorObservable ("Local Kink:Number");
      }

      if (measure_local_density_)
      {
        _states_cache.resize(wl.num_sites(),0.);
        measurements << alps::SimpleRealVectorObservable ("Local Density");
      }

      if (measure_local_density2_)
      {
        _states2_cache.resize(wl.num_sites(),0.);
        measurements << alps::SimpleRealVectorObservable ("Local Density^2");
      }

      // regarding on-fly measurements
      reinitialize_on_fly_measurements();

      measurements
        << alps::RealObservable  ("Green Function:0")
        << alps::RealObservable  ("Green Function:1")
        << alps::RealObservable  ("Momentum Distribution:0")
        ;

      if (measure_green_function_)
      {
        green.resize(num_sites(),0.);
        measurements << alps::SimpleRealVectorObservable ("Green Function");
      }

      if (finite_tof)
      {
        measurements << alps::RealObservable  ("Momentum Distribution:TOF:0");

        if (measure_green_function_)
        {
          green_tof.resize(num_sites(),0.);
          measurements << alps::SimpleRealVectorObservable ("Green Function:TOF");
        }

      }
    }

void
  directed_worm_algorithm
    ::reinitialize_on_fly_measurements()
    {
      green0 = 0.;
      green1 = 0.;
      nk0    = 0.;
      if (measure_green_function_)
        std::fill(green.begin(),green.end(),0.);

      if (finite_tof)
      {
        nk0_tof = 0.;
        if (measure_green_function_)
          std::fill(green_tof.begin(),green_tof.end(),0.);
      }
    }

inline boost::multi_array<double,4>
  directed_worm_algorithm
    ::bond_hamiltonian(const bond_descriptor & bond)
    {
      alps::Parameters this_parms;
      if(inhomogeneous_bonds()) {
        alps::throw_if_xyz_defined(parms,bond);      // check whether x, y, or z is set
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

      std::set<std::string> site_operators = model().bond_term(bond_type(bond)).operator_names(parms);

      for (std::set<std::string>::iterator site_operator_it = site_operators.begin(); site_operator_it != site_operators.end(); ++site_operator_it)
      {
        alps::Parameters this_parms;
        if(inhomogeneous_bonds()) {
          alps::throw_if_xyz_defined(parms,bond);      // check whether x, y, or z is set
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
    ::dostep()
    {
      ++_sweep_counter;

#ifdef DEBUGMODE
      // Step 0: Check whether worldlines configuration is valid? (DEBUG)
      if (_sweep_counter > 11179463)
        if (_sweep_counter % 1 == 0)
          if (wl.is_valid(state_maximum[0]))
            std::cout << "Worldlines configuration is valid up to the beginning of sweep " << _sweep_counter << "\n";
          else {
            std::cout << "Worldlines configuration is NOT valid at sweep " << _sweep_counter << "\n";
            std::cin.get();
          }
#endif

#define DEBUG_SWEEP_COUNTER 11179463

#ifdef DEBUGMODE
      if (_sweep_counter == DEBUG_SWEEP_COUNTER) 
      {
        std::cout << "\n\nSweep : " << _sweep_counter << "\t(Propagation sweep: " << _propagation_counter << ")"
                  << "\n============"
                  << "\n\n";

        std::cin.get();
      }
#endif

      // Step 1A: Generate a random wormpair
      const unsigned int _site     = random()*num_sites();
      const double       _time     = random();
      const bool         _forward  = (random()<0.5);
      const bool         _creation = (random()<0.5);

      location_type      _location = wl.location(_site,_time);
#ifdef DEBUGMODE
      if (_sweep_counter == DEBUG_SWEEP_COUNTER)
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

      if (!is_thermalized() && _sweep_counter % _skip == 0)
        std::cout << "Thermalizing " << _sweep_counter << " / " << _thermalization_sweeps << "\n";
                
      const bool checkpoint_sweep = (is_thermalized() && _sweep_counter % _skip == 0);

      // Step 1B: Check state to see if we can insert wormpair 
      const unsigned short _state = wl.state_before(_location);

      if (  (_state == state_minimum[site_site_type[_site]] && !wormpair::increasing(_forward, _creation)) 
         || (_state == state_maximum[site_site_type[_site]] && wormpair::increasing(_forward, _creation))  
         )
      {
        ++_sweep_failure_counter;      // worm insertion is unsuccessful         
        if (is_thermalized() && measure_ && _state != 0) 
        {
          green0 += _state;
          nk0    += _state;
          if (measure_green_function_)
            green[0] += _state;

          if (finite_tof)
          {
            nk0_tof += _state;
            if (measure_green_function_)
              green_tof[0] += _state;
          }
        }
#ifdef DEBUGMODE
        if (_sweep_counter == DEBUG_SWEEP_COUNTER)
        {
          std::cout << "\n\n"
                    << "\nWormpair creation fails!" 
                    << "\n========================";
          wl.output(std::cout, _site);
          std::cout << "\nEnded (Wormpair creation fails!)"
                    << "\n--------------------------------"
                    << "\n\n";
          //std::cin.get();
        }
#endif
        if (checkpoint_sweep) 
          docheckpoint();
        return;
      }

      // Step 1C: Insert wormpair
      worm = wormpair(_location, kink(_site,_time,_state), _forward, _creation);
#ifdef DEBUGMODE
      if (_sweep_counter == DEBUG_SWEEP_COUNTER)
      {
        std::cout << "\n\n"
                  << "\nWormpair creation!" 
                  << "\n==================";
        wl.output(std::cout, _site);
        std::cout << "\n\n"
                  << worm
                  << "\nEnded (Wormpair creation!)"
                  << "\n--------------------------------"
                  << "\n\n";
        //std::cin.get();
      }
#endif

      // Step 2: Wormhead propagates till collision with wormtail
      while(wormhead_propagates_till_collision_with_wormtail(worm.wormpair_state(), neighbors(worm.wormtail_site())));      
      const unsigned short termination_state = (worm.forward() ? worm.state_before() : worm.state());

      if (is_thermalized() && measure_ && termination_state != 0)
      {
        green0 += termination_state;
        nk0    += termination_state;
        if (measure_green_function_)
          green[0] += termination_state;
            
        if (finite_tof)
        {
          nk0_tof += termination_state;
          if (measure_green_function_)
            green_tof[0] += termination_state;
        }
      }

      worm.wormhead_annihilates_wormtail();
#ifdef DEBUGMODE
      if (_sweep_counter == DEBUG_SWEEP_COUNTER)
      {
        std::cout << "\n\n"
                  << "\nWormpair annihilation!"
                  << "\n======================";
        wl.output(std::cout, _site);
        std::cout << "\nEnded (Wormpair annihilation!)"
                  << "\n------------------------------------------------------"
                  << "\n\n";
        //std::cin.get();
      }
#endif

      if (checkpoint_sweep) 
        docheckpoint();
      return;
    }

void
  directed_worm_algorithm
    ::docheckpoint()
    {
#ifdef DEBUGMODE
      // check worldlines
      if (wl.is_valid(state_maximum[0])) 
        std::cout << " Valid : ";
#endif 

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
      if (measure_simulation_speed_)
      {
        std::time(&_simulation_timer.second);
        std::cout << " ... speed = " << std::difftime(_simulation_timer.second, _simulation_timer.first); 
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

      int    total_particle_number = 0;
      double total_energy_vertex   = 0.;
      double total_energy_onsite   = 0.;
      int    total_density2        = 0;
      double total_energy2         = 0.;

      for (unsigned int site=0; site < num_sites(); ++site)
      {
        total_particle_number += wl.site_state(site);
        const double this_energy_vertex = -static_cast<double>(wl.num_kinks(site)-1)/2;
        const double this_energy_onsite = onsite_energy(site, wl.site_state(site));
        total_energy_vertex   += this_energy_vertex;
        total_energy_onsite   += this_energy_onsite;
        if (measure_density2_)
          total_density2 += wl.site_state(site) * wl.site_state(site);
        if (measure_energy_density2_)
          total_energy2  += (total_energy_vertex+total_energy_onsite)*(total_energy_vertex+total_energy_onsite);
        if (measure_local_num_kinks_)
          _num_kinks_cache[site] = wl.num_kinks(site);
        if (measure_local_density_)
          _states_cache[site]    = wl.site_state(site);
        if (measure_local_density2_)
          _states2_cache[site]   = wl.site_state(site) * wl.site_state(site);
      }
      total_energy_vertex /= beta;
      total_energy_onsite /= beta;
      if (measure_energy_density2_)
        total_energy2 /= (beta*beta); 

      std::cout << " ... Measuring " 
                << " ... N = " << total_particle_number 
                ;
     
      double total_energy = total_energy_vertex + total_energy_onsite;

      measurements["Total Particle Number"]     << static_cast<double>(total_particle_number);
      measurements["Energy"]                    << total_energy;
      measurements["Energy:Vertex"]             << total_energy_vertex;
      measurements["Energy:Onsite"]             << total_energy_onsite;

      if (!inhomogeneous())
      {
        measurements["Density"]                 << static_cast<double>(total_particle_number)/num_sites();
        measurements["Energy Density"]          << total_energy/num_sites();
        measurements["Energy Density:Vertex"]   << total_energy_vertex/num_sites();
        measurements["Energy Density:Onsite"]   << total_energy_onsite/num_sites();
      }

      if (measure_number2_)
        measurements["Total Particle Number^2"] << static_cast<double>(total_particle_number * total_particle_number);
      if (measure_energy2_)
        measurements["Energy^2"]                << total_energy * total_energy;
      if (measure_density2_)
        measurements["Density^2"]               << static_cast<double>(total_density2)/num_sites();
      if (measure_energy_density2_)
        measurements["Energy Density^2"]        << total_energy2/num_sites();

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
        measurements["Winding Number^2"] << alps::numeric::vector2valarray(alps::numeric::sq(std::vector<double>(winding_number.begin(), winding_number.end())));
        measurements["Stiffness"]        << std::inner_product(winding_number.begin(), winding_number.end(), winding_number.begin(), 0.)/(dimension()*beta);
      }      

      if (measure_local_num_kinks_)
        measurements["Local Kink:Number"] << alps::numeric::vector2valarray(_num_kinks_cache);
      if (measure_local_density_) 
        measurements["Local Density"]     << alps::numeric::vector2valarray(_states_cache);
      if (measure_local_density2_)
        measurements["Local Density^2"]   << alps::numeric::vector2valarray(_states2_cache);

      // regarding on-fly measurements
      measurements["Green Function:0"]        << green0/_skip;
      measurements["Green Function:1"]        << green1/_skip;
      measurements["Momentum Distribution:0"] << nk0/_skip;

      if (measure_green_function_)
        measurements["Green Function"]      << alps::numeric::vector2valarray(green/_skip);

      if (finite_tof)
      {
        measurements["Momentum Distribution:TOF:0"] << nk0_tof/_skip;
        if (measure_green_function_)
          measurements["Green Function:TOF"]        << alps::numeric::vector2valarray(green_tof/_skip);
      }
    }

bool
  directed_worm_algorithm
    ::wormhead_propagates_till_collision_with_wormtail(unsigned short wormpair_state_, std::pair<neighbor_iterator,neighbor_iterator> const & neighbors_)
    {
      ++_propagation_counter;

#ifdef DEBUGMODE
      // Step 0: Check extended worldlines configuration
      if (_sweep_counter == DEBUG_SWEEP_COUNTER)
      {
        if (!worm.is_valid(state_maximum[0])) {
          std::cout << "Extended worldlines configuration is NOT valid at " << _propagation_counter << "\n";
          std::cin.get();
        }
      }
#endif

#define DEBUG_PROPAGATION_COUNTER 12383395

      // Step 1: Wormhead movement
      const double _time2next     = worm.time2next();
      const double _time2wormtail = worm.time2wormtail();

      const double _onsite_energy_relative = onsite_energy_relative(worm.site(), worm.state(), worm.forward(), worm.creation());

      double _deltatime = -std::log(1.- random())/_onsite_energy_relative;

#ifdef DEBUGMODE
      if (_sweep_counter == DEBUG_SWEEP_COUNTER)
      if (_propagation_counter >= DEBUG_PROPAGATION_COUNTER)
      {
        std::cout << "\n\n"
                  << "\nPropagation step : " << _propagation_counter
                  << "\n------------------------------------"
                  << "\nWorm :\n"
                  << worm
                  << "\n\n"
                  << "\nAttempting to move wormhead by a time of " << _deltatime << " ( " << _time2next << " / " << _time2wormtail << " ) "
                  << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                  << "\n\n";
        std::cin.get();
      }
#endif


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
      if (_sweep_counter == DEBUG_SWEEP_COUNTER)
      if (_propagation_counter >= DEBUG_PROPAGATION_COUNTER)
      {
        std::cout << "\nWormhead moves to time " << worm.time() << " ( " << (_halted ? "HALTED" : "NOT HALTED") << " ) "
                  << "\n============================================================================";
        wl.output(std::cout, worm.site());
        std::cout << "\n\n"
                  << "\nWormpair"
                  << "\n--------"
                  << "\n"
                  << worm
                  << "\n"
                  << "\nEnded (Wormhead movement)!"
                  << "\n--------------------------"
                  << "\n\n";
        std::cin.get();
      }
#endif

      if (is_thermalized() && measure_ && (_deltatime > _time2wormtail) && (worm.site() != worm.wormtail_site()))   
      {
        using alps::numeric::operator+;
        using alps::numeric::operator+=;
        using alps::numeric::operator*;

        int green_distance;
        std::vector<int> green_idx(dimension());
        if (measure_green_function_) {
          if (is_periodic_)
            green_distance = lattice().distance(worm.wormtail_site(), worm.site());
          else { 
            using boost::numeric::operators::operator-;
            std::vector<double> green_coordinate = alps::numeric::abs(coordinate(worm.wormtail_site()) - coordinate(worm.site()));
            // determine green_distance
            for (unsigned int i=0; i<dimension(); ++i)
              green_idx[i] = std::distance(component_coordinate_x100_lookup[i].begin(), std::lower_bound(component_coordinate_x100_lookup[i].begin(), component_coordinate_x100_lookup[i].end(), boost::math::iround(100*green_coordinate[i])));
            green_distance = site_lookup[componentIndex2Index(green_idx)];
          }
        }

        if (worm.neighbor2wormtail())
          green1 += wormpair_state_;
        nk0 += wormpair_state_;
        if (measure_green_function_)
          green[green_distance] += wormpair_state_;

        if (finite_tof)
        {
          const double tof_phase = phase(worm.site())-phase(worm.wormtail_site()); 
          nk0_tof += wormpair_state_ * std::cos(tof_phase);
          if (measure_green_function_)
            green_tof[green_distance] += wormpair_state_ * std::cos(tof_phase);
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
          if (_sweep_counter == DEBUG_SWEEP_COUNTER)
          if (_propagation_counter >= DEBUG_PROPAGATION_COUNTER)
          {
            std::cout << "\n\n"
                      << "\nWormhead crosses kink"
                      << "\n============================================================================"
                      << "\n";
            wl.output(std::cout, worm.site());
            std::cout << "\n\n"
                      << worm
                      << "\n\n"
                      << "\nEnded (Wormhead crosses kink)!"
                      << "\n--------------------------"
                      << "\n\n";
            std::cin.get();
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
      if (  ( worm.increasing() && _targetstate == state_maximum[site_site_type[target(*it)]])
         || (!worm.increasing() && _targetstate == state_minimum[site_site_type[target(*it)]])
         || (!wl.location_is_kink_unoccupied(_neighborlocation, worm.time()))
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

#ifdef DEBUGMODE
      if (_sweep_counter == DEBUG_SWEEP_COUNTER)
      if (_propagation_counter >= DEBUG_PROPAGATION_COUNTER)
      {
        std::cout << "\n\nDelete/relink/bounce?"
                  << "\n======================="
                  << "\nworm :\n"
                  << worm << "\n"
                  << "\nwl(wormsite):"
                  << "\n-------------";
        wl.output(std::cout, worm.site());
        std::cout << "\n-------------";
        std::cout << "\nwl(sourcesite):"
                  << "\n-------------";
        wl.output(std::cout, _sourcelocation.first->begin()->siteindicator());
        std::cout << "\n-------------";
        std::cin.get();
      }
#endif


#ifndef HEATBATH_ALGORITHM
      const std::pair<neighbor_bond_iterator, neighbor_bond_iterator> _neighbor_bonds = neighbor_bonds(worm.next_partnersite());

      double _weight;
      for (neighbor_bond_iterator it=_neighbor_bonds.first; it!=_neighbor_bonds.second; ++it)
        if (target(*it) == worm.site())
        {
          _weight = hopping_energy(*it, (wl.location(target(*it), worm.next_time()).second-1)->state(), !worm.increasing());
          break;
        }

#ifdef DEBUGMODE
      if (_sweep_counter == DEBUG_SWEEP_COUNTER)
      if (_propagation_counter >= DEBUG_PROPAGATION_COUNTER)
      {
        std::cout << "\nWeight (old) = " << _weight << "\n";
        std::cin.get();
      }
#endif
     
      neighbor_bond_iterator it=_neighbor_bonds.first + random()*num_neighbors(worm.next_partnersite());

      if (target(*it) == worm.site())    // delete vertex or bounce
      {
        const double _weight_new = onsite_energy_relative(worm.next_partnersite(),_sourcelocation.second->state(),!worm.forward(),worm.creation());

#ifdef DEBUGMODE
        if (_sweep_counter == DEBUG_SWEEP_COUNTER)
        if (_propagation_counter >= DEBUG_PROPAGATION_COUNTER)
        {
          std::cout << "\nWeight (new) = " << _weight_new << "\n";
          std::cin.get();
        }
#endif

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
           || (!wl.location_is_kink_unoccupied(_neighborlocation, worm.next_time()))
           )
        {
          worm.wormhead_turns_around();
          ++_propagation_failure_counter;
          return;
        }
        else
        {
          const double _weight_new = hopping_energy(*it, _targetstate, _increasing);

#ifdef DEBUGMODE
          if (_sweep_counter == DEBUG_SWEEP_COUNTER)
          if (_propagation_counter >= DEBUG_PROPAGATION_COUNTER)
          {
            std::cout << "\nWeight (new) = " << _weight_new << "\n";
            std::cin.get();
          }
#endif

          if (_weight_new >= _weight || random() < (_weight_new/_weight))
          {
#ifdef DEBUGMODE
            if (_sweep_counter == DEBUG_SWEEP_COUNTER)
            if (_propagation_counter >= DEBUG_PROPAGATION_COUNTER)
            {
              std::cout << "\nwl(wormsite):"
                        << "\n-------------";
              wl.output(std::cout, worm.site());
              std::cout << "\n-------------";
              std::cout << "\nwl(sourcesite):"
                        << "\n-------------";
              wl.output(std::cout, _sourcelocation.first->begin()->siteindicator());
              std::cout << "\n-------------";
              std::cout << "\nwl(targetsite):"
                        << "\n-------------";
              wl.output(std::cout, _neighborlocation.first->begin()->siteindicator());
              std::cout << "\n-------------";
              std::cin.get();
            }
#endif 
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
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
   return alps::scheduler::start(argc,argv,alps::scheduler::SimpleMCFactory<directed_worm_algorithm>());
#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (std::exception& exc) {
    std::cerr << exc.what() << "\n";
      alps::comm_exit(true);
      return -1;
    }
  catch (...) {
    std::cerr << "Fatal Error: Unknown Exception!\n";
    return -2;
  }
#endif
}

/*
 * Luka's NGS scheduler is not ready yet, Tama is forced to revert back to the old scheduler.
 *
int main(int argc, char** argv)
{
  // read options from command line
  alps::mcoptions options(argc, argv);

  // initialize simulation and load the parameters from input file
  alps::hdf5::archive simulation_input = alps::hdf5::archive(options.input_file);
  directed_worm_algorithm simulation(simulation_input);

  // try loading whatever is there in the input file
  simulation.load(simulation_input);
 
  // run the simulation
  simulation.run(alps::stop_callback(options.time_limit));

  // save worldlines and measurement to output file
  alps::hdf5::archive simulation_output = alps::hdf5::archive(options.output_file, 'a');
  simulation.save(simulation_output);

  return 0;
}
*/
