/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2005 by Fabien Alet <alet@comp-phys.org>,
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

/* $Id: SSE.hpp,v 1.26 2006/04/29 01:32:28 troyer Exp $ */

#ifndef SSE_HPP
#define SSE_HPP


// Here do define
// #define AUTO_TIMES

#include "../qmc.h"
#include "SSE.Classes.hpp"
#include "SSE.Histogram.hpp"
#include <boost/cstdint.hpp>
#include <boost/multi_array.hpp>
#include <boost/vector_property_map.hpp>
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
        bool is_thermalized() const {return thermalized;}
        double work_done() const;
        bool change_parameter(const std::string&, const StringValue&);
        void store_histograms(string);
    
protected :
    // Common parameters
    histogram<double> final_weight_fract;
        histogram<double> histo;  //Histogram for current_number_of_non_identity
           histogram<double> histoup; //Histogram for upwalker
    vector<double> timeup;  // Time up (measured as proposed by Stefan Wessel)
    vector<double> timedown;// Time down (measured as proposed by Stefan Wessel)
    vector<double> realtimeup;
    vector<double> realtimedown;
    double thistime;
    double realtime;
    bool variable_length_representation;
        uint32_t upwalker;
    uint32_t number_of_worms_per_sweep; 
    int simulation_phase;  //1: WL-RUN, 2: OE-RUN (one step)

    // WL parameters
    double logf;
        unsigned int logf_steps_total;
        unsigned int logf_step;
        double minimum_histogram;
        double flatness_treshold;
        uint32_t use_zhou_bhatt;
    histogram<double> g; // Density of states for current_number_of_non_identity    
    
    // OE parameters
        unsigned int OE_nb_thermalization;     
       

    void finish_measurements();
    void store_histo_vars();
        void do_WL_step();
        void do_OE_step();
        void update_final_weights();
       // virtual bool almost_thermalized();
    void create_observables();

    void initialize_current_simulation();
        double ensemble_weight_fraction(int norder);
    void load_external_data();

        typedef QMCRun<> super_type;
        /************* Worm weights ******************/
        double worm_plus(uint8_t s,uint8_t t) 
        { 
                return s< matrix_element_raise_[t].size() ? 
                        matrix_element_raise_[t][s] : 0.;
        }
        
        double worm_minus(uint8_t s,uint8_t t) 
        { 
                return s< matrix_element_lower_[t].size() ? 
                        matrix_element_lower_[t][s] : 0.;
        }
        
        /******** Simulation parameters ***********************
         * It may be desirable todo more than ~4.3e9 sweeps. *
         * Hence use data type uint64_t for nb_steps and     *
         * steps_done_total                                   *
         *****************************************************/
        
        double epsilon;
        double total_shift;
        uint32_t cutoff_L; 
        uint32_t norder_min;
        uint32_t norder_max;
        
        uint64_t block_sweeps;
        uint64_t block_sweeps_total;
        uint64_t nb_measurement_steps;
        uint64_t measurements_done;
        uint64_t nb_steps;
        uint64_t steps_done_total;
        unsigned int all_done;
       
        bool perform_measurements;
        bool stop_thermalization;
        bool thermalized;
                       
        double worm_size;
        bool FORCE_STANDARD_DIRECTED_LOOPS;
        bool FORCE_HEATBATH;
        double nb_worms_thermalization;
        int count_worms_thermalization;

        int number_of_bond_types; 
        uint8_t maximum_number_of_states;
        vector<uint8_t> site_number_of_states; 
        vector<uint8_t> site_state; 

        vector<std::pair<uint8_t,uint8_t> > site_type_for_bond_type_;
        std::map<int,int> number_states_for_site_type_;
        vector<uint32_t> number_of_bonds_for_bond_type;

        std::vector<std::vector<double> > matrix_element_raise_;
        std::vector<std::vector<double> > matrix_element_lower_;
        std::vector<std::vector<double> > matrix_element_n_;
        std::vector<std::vector<double> > matrix_element_Sz_;

        std::map<int,boost::multi_array<double,4> > matrix_element;
        boost::multi_array<double, 3> proba_diagonal;
        boost::multi_array<double, 7> proba_worm;
 
        vector<double> energy_offset;
        vector<vertex_type> operator_string;
        vector<vertex_type> buffer_string;      //These two are used for variable length only
        vector<vertex_type> identity_string;

        uint32_t current_number_of_non_identity;
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
        void initialize_site_states();
        void initialize_hamiltonian();
        void print_arrays();
        void determine_bonds_offset();
        void initialize_diagonal_update_probabilities();
        void initialize_simulation();
      

        /******** SSE.Update.cpp ********/
        void print_operator_string();
        void diagonal_update();
        virtual void diagonal_update_iteration_end(int i);
        virtual void diagonal_update_sweep_end();

        void link_legs();
        std::pair<uint32_t,uint8_t> initial_vertex();
        void worm_update();
        inline uint8_t return_exit_leg(uint32_t,uint8_t);
        void do_update();

        

        /******** SSE.Directed.cpp ********/
        double the_sol[4][4];
        bool check_out(double[4]);
        bool find_sol(double[4],double[4],double[4],double[4]);
        double return_weight_after_flip(uint32_t, uint8_t[4],uint8_t,bool);
        double return_weight_before_flip(uint32_t, uint8_t[4],uint8_t,bool);
        double return_worm_weight_before(uint32_t,uint8_t[4],uint8_t,bool);
        double return_worm_weight_after(uint32_t,uint8_t[4],uint8_t,bool);
        bool calculate_scattering_matrix(uint32_t,uint8_t[4],bool);
        bool calculate_all_scattering_matrices(uint32_t,bool);            
        void put_in(uint32_t,uint8_t[4],uint8_t);
        void calculate_proba_worm();
        void calculate_heat_bath_matrix(uint32_t,uint8_t[4],bool);
        
        /******** SSE.Measurements.cpp ********/
       // void create_observables();
        void do_measurements();
  
        std::vector<unsigned int> original_bond_type;
        
};  

#endif
