/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2010 by Sergei Isakov <isakov@itp.phys.ethz.ch>
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

#ifndef __SSE_ALG_H__
#define __SSE_ALG_H__

#include <cmath>
#include <vector>
#include <limits>

#include <alps/osiris/dump.h>

#include "lattice.h"
#include "model.h"
#include "sse_alg_def.h"
#include "sse_worm_prob.h"

template<typename L, typename M, typename W>
class SSE_alg {
public:
    typedef L lattice_type;
    typedef M model_type;
    typedef W worker_type;
    
    typedef typename lattice_type::lat_unit_type lat_unit_type;
    typedef typename lattice_type::lat_units_type lat_units_type;
    typedef typename lattice_type::lat_unit_sites_type lat_unit_sites_type;
    typedef typename model_type::vertex_type vertex_type;
    typedef typename model_type::lat_unit_state_type lat_unit_state_type;
    
    typedef typename worker_type::state_type state_type;
    
    typedef Worm_prob<lattice_type, model_type> worm_prob_type;
    typedef Measurement<lattice_type, model_type, worker_type> measurement_type;
    
    typedef typename worm_prob_type::wprob_type wprob_type;
    typedef typename worm_prob_type::wprob_table_type wprob_table_type;
    
    SSE_alg(lattice_type const& lattice, model_type const& model,
            worker_type& worker, alps::Parameters const& params) :
        lattice(lattice),
        model(model),
        worker(worker),
        lat_units(lattice.lat_units()),
        nbstates(model.nbstates()),
        measurement(lattice, model, worker, params, state, opstring),
        green(measurement.get_green())
    {
        nworms = params.value_or_default("NUMBER_OF_WORMS_PER_SWEEP", 0);
        nops = params.value_or_default("INITIAL_CUTOFF", 10);
        abort_factor = params.value_or_default("WORM_ABORT", 0.0);
        
        if (abort_factor != 0.0)
            worker.measure_green_function(false);
            
        measure_green_function = false;
        
        nnonzero = 0;
        
        cc = 0.0;
        for (unsigned i = 0; i < lat_units.size(); ++i)
            cc += model.c(lat_units[i].type);
        cc *= worker.beta();

        // worm probabilities
        worm_prob_type worm_prob(params, lattice, model);
        worm_prob.calc_probabilities(wprobs, -1, 1);
        
        // diagonal probabilities
        d10_probs.resize(model.nvertices());
        d01_probs.resize(model.nvertices());
        for (unsigned vi = 0; vi < model.nvertices(); ++vi) {
            vertex_type const& vertex = model.vertex(vi);
            if (!vertex.diagonal)
                continue;

            double w = model.c(vertex.unit_type) - vertex.me;
            if (w < 1e-15) {
                // zero weight
                d01_probs[vi] = 0.0;
                continue;
            }

            d01_probs[vi] = w * worker.beta() * lat_units.size();
            d10_probs[vi] = 1.0 / d01_probs[vi];
        }
        
        nsites = lattice.nsites();
                
        state.resize(nsites);
        first.resize(nsites);
        last.resize(nsites);
        opstring.resize(nops);
        op_indices.resize(nops);
        
        // startconf
        for (unsigned i = 0; i < nsites; ++i)
            state[i] = worker.mrandom_int(nbstates[lattice.sitei2alps_type(i)]);
        for (op_iterator op = opstring.begin() ; op != opstring.end(); ++op)
            op->vertex_index = IDENTITY;
            
        nworms_therm = 0;
        count_therm = 0;
    }
    
    void save(alps::ODump& dump) const
    {
        dump << nworms << state << opstring << nops << nnonzero
                 << nworms_therm << count_therm;
    }
    
    void load(alps::IDump& dump)
    {
        dump >> nworms >> state >> opstring >> nops >> nnonzero
                >> nworms_therm >> count_therm;
        op_indices.resize(nops);
    }
    
    void do_step()
    {
        diagonal_update();
        
        if (!worker.is_thermalized()) {
            if (nnonzero == 0)
                update_free_states();
            else {
                if (nworms > 0) {
                    for (unsigned i = 0; i < nworms; ++i)
                        worm_update();
                } else {
                    // adjust nworms
                    double len1 = 0.0;
                    boost::uint64_t nworms_therm2 = 0;
                    do {
                        ++nworms_therm2;
                        len1 += worm_update();
                    } while (len1 < 2 * nnonzero);

                    if (worker.is_thermalization_done(0.8)) {
                        nworms_therm += nworms_therm2;
                        ++count_therm;
                    }
                }

                opstring2state();
                adjust_opstring_len();
            }
        } else {
            measure_green_function = worker.measure_green_function();
                    
            if (nworms == 0) {
                if (count_therm > 0)
                    nworms = unsigned(nworms_therm / count_therm);
                else
                    nworms = 10;
                std::cout << "NUMBER_OF_WORMS_PER_SWEEP: " << nworms << "\n";
            }
            measurement.set_nworms(nworms);
            
            if (nnonzero == 0)
                update_free_states();
                
            for (unsigned i = 0; i < nworms; ++i)
                worm_update();

            opstring2state();
        }
    }
        
    void do_measurement()
    {
        measurement.do_measurement(nops, nnonzero, cc);
    }        
private:
    SSE_alg();
    
    lattice_type const& lattice;
    model_type const& model;
    worker_type& worker;
    lat_units_type const& lat_units;
    std::vector<unsigned> const& nbstates;
    
    measurement_type measurement;
    
    double cc;
    
    unsigned nworms;
    
    boost::uint64_t nworms_therm;
    boost::uint64_t count_therm;
    
    unsigned nsites;
    
    std::vector<state_type> state;
    std::vector<Operator> opstring;
    std::vector<Operator> opstring_copy;
    std::vector<unsigned> op_indices;

    std::vector<unsigned> first;
    std::vector<unsigned> last;
    
    wprob_table_type wprobs;
    std::vector<double> d01_probs;
    std::vector<double> d10_probs;
    
    unsigned nops;
    unsigned nnonzero;
    
    double abort_factor;
    
    bool measure_green_function;
    std::valarray<double>& green;
    
    void diagonal_update()
    {
        std::memset(&first[0], 0xff, nsites * sizeof(unsigned));
        std::memset(&last[0], 0xff, nsites * sizeof(unsigned));

        unsigned oi = 0;
        unsigned opcount = 0;
        for (op_iterator op = opstring.begin(); op != opstring.end(); ++op, ++oi) {
            if (op->vertex_index == IDENTITY) {
                // try 0->1 update
                
                unsigned ui = worker.mrandom_int(lat_units.size());
                lat_unit_type const& lat_unit = lat_units[ui];

                lat_unit_state_type unit_state;
                for (unsigned i = 0; i < UNIT_SIZE; ++i)
                    unit_state[i] = state[lat_unit.sites[i]];

                unsigned vertex_index = model.diag_vertex_index(unit_state,
                    lat_unit.type, lat_unit.sites);
                if (vertex_index == model_type::INVALID_VERTEX)
                    continue;
                                    
                double p = d01_probs[vertex_index];
                if (p >= (nops - nnonzero)
                        || p > (worker.mrandom_real() * (nops - nnonzero))) {
                    op->vertex_index = vertex_index;
                    op->unit_ref = ui;
                    ++nnonzero;
                } else
                    continue;
            } else {
                lat_unit_type const& lat_unit = lat_units[op->unit_ref];
                vertex_type const& vertex = model.vertex(op->vertex_index);

                if (vertex.diagonal) {
                    // try 1->0 update

                    double p = d10_probs[op->vertex_index] * (nops - nnonzero + 1);
                    if (p >= 1.0 || p >= worker.mrandom_real()) {
                        op->vertex_index = IDENTITY;
                        --nnonzero;
                        continue;
                    }
                } else
                    // propagate state
                    for (unsigned i = 0; i < UNIT_SIZE; ++i)
                        state[lat_unit.sites[i]] = vertex.state[UNIT_SIZE + i];
            }
                    
            // build vertex list
            op_indices[opcount++] = oi;
            unsigned cur_leg = 2 * UNIT_SIZE * oi;
            lat_unit_sites_type const& sites = lat_units[op->unit_ref].sites;
            for (unsigned j = 0; j < UNIT_SIZE; ++j) {
                unsigned s = sites[j];

                unsigned last_leg = last[s];

                if (last_leg == MAX_NUMBER)
                    first[s] = cur_leg;
                else {
                    const unsigned S2 = 2 * UNIT_SIZE;
                    opstring[cur_leg / S2].linked[cur_leg % S2] = last_leg;
                    opstring[last_leg / S2].linked[last_leg % S2] = cur_leg;
                }

                last[s] = cur_leg + UNIT_SIZE;
                ++cur_leg;
            }
        }
        
        // vertex list --- force periodic boundary conditions
        for (unsigned i = 0; i < nsites; ++i) {
            unsigned first_leg = first[i];
            if (first_leg == MAX_NUMBER)
                continue;
            
            unsigned last_leg = last[i];
            const unsigned S2 = 2 * UNIT_SIZE;
            opstring[first_leg / S2].linked[first_leg % S2] = last_leg;
            opstring[last_leg / S2].linked[last_leg % S2] = first_leg;
        }

        if (worker.is_thermalized() && nnonzero == nops)
            throw std::runtime_error("nnonzero reached the maximum value. "
                    "Please increase THERMALIZATION or set INITIAL_CUTOFF "
                    "to some large enough value.");
    }
    
    boost::uint64_t worm_update()
    {
        unsigned start_site, start_level, start_vertex, start_leg;
        
        if (measure_green_function) {
            if (gf_start_position(start_site, start_level, start_vertex, start_leg)) {
                // only measure green fucntion
                gf_measure_free_site(start_site);
                return 0;
            }
        } else {
            if (nnonzero == 0)
                return 0;
            
            start_vertex = op_indices[worker.mrandom_int(nnonzero)];
            start_leg = worker.mrandom_int(2 * UNIT_SIZE);
            start_level = 0; // unused

            Operator& start_op = opstring[start_vertex];
            lat_unit_type const& lat_unit = lat_units[start_op.unit_ref];
            start_site = lat_unit.sites[start_leg % UNIT_SIZE];
        }

        vertex_type const& vertex = model.vertex(opstring[start_vertex].vertex_index);
        unsigned start_site_type = lattice.sitei2alps_type(start_site);

        unsigned head_op = worker.mrandom_int(2);
        unsigned nbstates = this->nbstates[start_site_type];
        // here we can have: if (nbstates == 2) head_op = 1 - vertex.state[start_leg];
        if (head_op == 0) {
            if (vertex.state[start_leg] == 0)
                return 0;
        } else {
            if (vertex.state[start_leg] == nbstates - 1)
                return 0;
        }
        
        unsigned max_len = 0;
        if (abort_factor > 0.0) {
            opstring_copy = opstring;
            max_len = static_cast<unsigned>(abort_factor * nops);
        }
        
        double worm_weight = 0.0;
        
        if (measure_green_function)
            worm_weight = gf_start(head_op, start_site, start_site_type,
                start_leg, vertex);

        boost::uint64_t len = 0;

        unsigned cur_leg = start_leg;
        unsigned cur_vertex = start_vertex;
        while (1) {
            ++len;
            
            if (max_len > 0 && len > max_len) {
                std::cerr << "aborting worm update\n";
                opstring.swap(opstring_copy);
                break;
            }
            
            Operator& cur_op = opstring[cur_vertex];
            unsigned& vertex_index = cur_op.vertex_index;
                
            unsigned k = worm_prob_type::w_index(vertex_index, cur_leg, head_op);
            std::vector<wprob_type> const& wprob = wprobs[k];

            if (wprob.size() == 0)
                throw std::runtime_error("Internal error: size of prob table is zero.");

            typename std::vector<wprob_type>::const_iterator wp = wprob.begin();
            typename std::vector<wprob_type>::const_iterator last = wprob.end() - 1;

            if (wp != last) {
                double p = worker.mrandom_real();

                for (; wp != last; ++wp)
                    if (p < wp->prob) break;
            }

            head_op = wp->ex_op;
            vertex_index = wp->vertex_index;

            if (cur_vertex == start_vertex && wp->ex_leg == start_leg)
                break;

            // next vertex and entrance leg
            unsigned prev_vertex = cur_vertex;
            cur_vertex = cur_op.linked[wp->ex_leg] / (2 * UNIT_SIZE);
            cur_leg = cur_op.linked[wp->ex_leg] % (2 * UNIT_SIZE);

            if (cur_vertex == start_vertex && cur_leg == start_leg)
                break;
                
            if (measure_green_function)
                gf_next(start_site, start_level, cur_vertex, cur_leg,
                    prev_vertex, worm_weight);
        }
                
        return len;
    }
    
    bool gf_start_position(unsigned& start_site, unsigned& start_level,
        unsigned& start_vertex, unsigned& start_leg)
    {
        start_site = worker.mrandom_int(nsites);
        
        if (first[start_site] == MAX_NUMBER)
            return true;
        
        start_level = worker.mrandom_int(nops);
        
        unsigned cur_vertex;
        
        if (start_level < nops / 2) {
            cur_vertex = start_vertex = first[start_site] / (2 * UNIT_SIZE);
            start_leg = first[start_site] % (2 * UNIT_SIZE);
            while (cur_vertex < start_level) {
                unsigned next = opstring[cur_vertex].linked[start_leg + UNIT_SIZE];
                cur_vertex = next / (2 * UNIT_SIZE);
                start_leg = next % (2 * UNIT_SIZE);
                if (cur_vertex == start_vertex) break;
            }
        } else {
            cur_vertex = start_vertex = last[start_site] / (2 * UNIT_SIZE);
            start_leg = last[start_site] % (2 * UNIT_SIZE);
            while (cur_vertex >= start_level) {
                unsigned next = opstring[cur_vertex].linked[start_leg - UNIT_SIZE];
                cur_vertex = next / (2 * UNIT_SIZE);
                start_leg = next % (2 * UNIT_SIZE);
                if (cur_vertex == start_vertex) break;
            }
        }
        
        if (worker.mrandom_real() < 0.5) {
            unsigned prev = opstring[cur_vertex].linked[start_leg];
            cur_vertex = prev / (2 * UNIT_SIZE);
            start_leg = prev % (2 * UNIT_SIZE);
        }
        
        start_vertex = cur_vertex;
        
        return false;
    }
    
    void gf_measure_free_site(unsigned start_site)
    {
        unsigned stype = lattice.sitei2alps_type(start_site);
        unsigned state = this->state[start_site];
        double mer = model.raising_matrix_elements()[stype][state];
        double mel = model.lowering_matrix_elements()[stype][state];

        double me = 0.5 * (mer * mer + mel * mel);

        if (!worker.do_measurement_origin())
            green[lattice.distance(start_site, start_site)] += me;
        else if (start_site == worker.measurement_origin())
            green[start_site] += me;
    }
    
    double gf_start(unsigned head_op, unsigned start_site,
        unsigned start_site_type, unsigned start_leg,
        vertex_type const& vertex)
    {
        std::vector<std::vector<double> > const* mes;
        if (head_op == 0)
            mes = &model.lowering_matrix_elements();
        else
            mes = &model.raising_matrix_elements();
        double worm_weight = (*mes)[start_site_type][vertex.state[start_leg]];
        
        worm_weight *= worm_weight;

        if (!worker.do_measurement_origin())
            green[lattice.distance(start_site, start_site)] += worm_weight;
        else if (start_site == worker.measurement_origin())
            green[start_site] += worm_weight;
            
        return worm_weight;
    }
    
    void gf_next(unsigned start_site, unsigned start_level,
        unsigned cur_vertex, unsigned cur_leg, unsigned prev_vertex,
        double const& worm_weight)
    {
        lat_unit_type const& lat_unit =
                lat_units[opstring[cur_vertex].unit_ref];
        unsigned site = lat_unit.sites[cur_leg % UNIT_SIZE];
        
        bool crossed = false;
        if (cur_leg < UNIT_SIZE) {
            if (cur_vertex <= prev_vertex) {
                if (cur_vertex >= start_level || prev_vertex < start_level)
                    crossed = true;
            } else if (cur_vertex >= start_level && prev_vertex < start_level)
                crossed = true;
        } else {
            if (cur_vertex >= prev_vertex) {
                if (cur_vertex < start_level || prev_vertex >= start_level)
                    crossed = true;
            } else if (cur_vertex < start_level && prev_vertex >= start_level)
                crossed = true;
        }
        
        if (crossed) {
            if (!worker.do_measurement_origin())
                green[lattice.distance(start_site, site)] += worm_weight;
            else if (start_site == worker.measurement_origin())
                green[site] += worm_weight;
        }
    }
        
    void update_free_states()
    {
        for (unsigned i = 0; i < nsites; ++i)
            state[i] = worker.mrandom_int(nbstates[lattice.sitei2alps_type(i)]);
    }
    
    void opstring2state()
    {
        for (unsigned i = 0; i < nsites; ++i) {
            unsigned first_pos = first[i];

            if (first_pos != MAX_NUMBER) {
                unsigned cur_vertex = first_pos / (2 * UNIT_SIZE);
                unsigned cur_leg = first_pos % (2 * UNIT_SIZE);
                state[i] = model.vertex_state(opstring[cur_vertex].vertex_index)[cur_leg];
            } else
                // update "free" states
                state[i] = worker.mrandom_int(nbstates[lattice.sitei2alps_type(i)]);
        }
    }
    
    void adjust_opstring_len()
    {
        if (nnonzero < unsigned(0.8 * nops))
            return;
            
        unsigned nops1 = unsigned(1.1 * nops);
        if (nops1 == nops) ++nops1;

        opstring.reserve(nops1);
        opstring.resize(nops1);
        op_indices.reserve(nops1);
        op_indices.resize(nops1);

        unsigned j = 0;
        for (unsigned i = 0; i < nops; ++i)
            if (opstring[i].vertex_index != IDENTITY) {
                opstring[j] = opstring[i];
                ++j;
            }
        for (unsigned i = nnonzero; i < nops1; ++i)
            opstring[i].vertex_index = IDENTITY;

        j = nnonzero - 1;
        for (unsigned i = nops1 - 1; j < i; --i)
            if (j > 0 && nnonzero > (unsigned) worker.mrandom_int(nops1)) {
                opstring[i] = opstring[j];
                opstring[j].vertex_index = IDENTITY;
                --j;
            }

        nops = nops1;
    }
};

#endif
