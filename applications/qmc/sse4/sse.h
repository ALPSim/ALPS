/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2009-2010 by Sergei Isakov <isakov@itp.phys.ethz.ch>
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

#ifndef __SSE_RUN_H__
#define __SSE_RUN_H__

#include <alps/scheduler/montecarlo.h>
#include <alps/scheduler/measurement_operators.h>
#include <alps/osiris/dump.h>
#include <alps/lattice.h>
#include <alps/alea.h>

#include "../qmc.h"
#include "lattice.h"
#include "model.h"
#include "measurement.h"
#include "sse_alg.h"

class SSE_run : public QMCRun<> {    
public:
    typedef SSE_run self_type;
    typedef QMCRun<> super_type;
    typedef Lattice<self_type> lattice_type;
    typedef Model<self_type, Lattice<self_type> > model_type;
    typedef SSE_alg<lattice_type, model_type, self_type> algorithm_type;
    
    SSE_run(alps::ProcessList const& w, alps::Parameters const& params,
            int n, bool issymbolic = false) :
        super_type(w, params, n, issymbolic),
        lattice(params, *this),
        model(params, *this, lattice),
        algorithm(lattice, model, *this, params)
    {
        nsweeps = params.defined("SWEEPS") ?
            boost::uint64_t(params["SWEEPS"]) : (params.defined("MCS") ? 
                boost::uint64_t(params["MCS"]) : boost::uint64_t(params["Steps"]));
        
        nthermalization = params.defined("THERMALIZATION") ?
                boost::uint64_t(params["THERMALIZATION"]) : (params.defined("thermalization") ?
                boost::uint64_t(params["thermalization"]) : nsweeps / 10);
                
        measurements_skip = parms.value_or_default("SKIP", 1);

        if (is_signed())
            std::cout << "Warning: hamiltonian has a sign problem.\n";
            
        sweeps_done = 0;
        measurements_done = 0;
    }
    
    void save(alps::ODump& dump) const
    {
        dump << nthermalization << sweeps_done << measurements_done;
        algorithm.save(dump);
    }
    
    void load(alps::IDump& dump)
    {
        if (super_type::where.empty())
            super_type::measurements.compact();
        else {
            dump >> nthermalization >> sweeps_done >> measurements_done;
            algorithm.load(dump);
        }
    }
    
    void dostep()
    {
//        if (sweeps_done >= nthermalization + nsweeps)
//            return;
        
        algorithm.do_step();
        
        if (is_thermalized() && ++measurements_done == measurements_skip) {
            algorithm.do_measurement();
            measurements_done = 0;
        }

        ++sweeps_done;
    }
    
    bool is_thermalized() const
    {
        return sweeps_done >= nthermalization;
    }
    
    double work_done() const
    {
        if (is_thermalized())
            return (sweeps_done - nthermalization) / double(nsweeps);
        return 0.;
    }
    
    bool change_parameter(std::string const& name, alps::StringValue const& value)
    {
        boost::uint64_t new_nsweeps = 0;
        
        if (name == "SWEEPS")
            new_nsweeps = boost::uint64_t(value);
        if (name == "MCS")
            new_nsweeps = boost::uint64_t(value);
        if (name == "Steps")
            new_nsweeps = boost::uint64_t(value);

        if (new_nsweeps > 0) {
            nsweeps = new_nsweeps;
            return true;
        }
        
        // cannot change anyting else
        return false;
    }
    
    static void print_copyright(std::ostream& out)
    {
        out << "Quantum Monte Carlo simulations using the SSE algorithm v. 4.1.1\n"
            << "  available from http://alps.comp-phys.org/\n"
            << "  copyright (c) 2003-2010 by Sergei Isakov <isakov@itp.phys.ethz.ch>\n\n";
    }
    
    bool is_thermalization_done(double percentage) const
    {
        return double(sweeps_done) > percentage * nthermalization;
    }
    
    int mrandom_int(int n)
    {
        return super_type::random_int(n);
    }
    
    double mrandom_real()
    {
        return super_type::random_01();
    }
    
    alps::ObservableSet& measurements()
    {
        return super_type::measurements;
    }
    
    bool is_signed() const
    {
        return super_type::is_signed_;
    }
    
    void measure_green_function(bool flag)
    {
        super_type::measure_green_function_ = flag;
    }
    
    bool measure_green_function() const
    {
        return super_type::measure_green_function_;
    }
    
    bool do_measurement_origin() const
    {
        if (super_type::measurement_origin_)
            return true;
        else
            return false;
    }
    
    int measurement_origin() const
    {
        return super_type::measurement_origin_.get();
    }
    
    std::vector<unsigned> const& distance_mult() const
    {
        return super_type::distance_mult;
    }
    
    bool measure_site_compressibility() const
    {
        return super_type::measure_site_compressibility_;
    }
    
    bool do_common_measurements(double sign,
        std::vector<state_type> const& state, std::valarray<double> const& localint)
    {
        return super_type::do_common_measurements(sign, state, localint);
    }
    
    void create_common_observables()
    {
        super_type::create_common_observables();
    }
    
    void initialize_site_states()
    {
        super_type::initialize_site_states();
    }
    
    double beta() const
    {
        return super_type::beta;
    }
    
    std::vector<std::vector<double> > const& phys_states_n()
    {
        return super_type::diagonal_matrix_element["n"];
    }
private:
    lattice_type lattice;
    model_type model;
    algorithm_type algorithm;
    
    boost::uint64_t nsweeps;
    boost::uint64_t nthermalization;
    boost::uint64_t sweeps_done;
    unsigned measurements_skip;
    unsigned measurements_done;
};

#endif
