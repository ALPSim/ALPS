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

#ifndef __SSE_WORM_PROB__
#define __SSE_WORM_PROB__

#include <cmath>
#include <stack>
#include <vector>
#include <algorithm>
#include <numeric>

#include "lp_sse.h"
#include "lattice.h"
#include "model.h"

template<typename L, typename M>
class Worm_prob {
public:
    struct Wprob {
        double prob;             // probability to update
        unsigned ex_leg;         // exit leg
        unsigned ex_op;          // new operator head
        unsigned vertex_index;   // new vertex index

        bool operator() (Wprob const&  wp1, Wprob const& wp2) const
        {
            return wp1.prob > wp2.prob;
        }
    };
    
    typedef L lattice_type;
    typedef M model_type;
    
    typedef typename lattice_type::lat_unit_type lat_unit_type;
    typedef typename lattice_type::lat_unit_sites_type lat_unit_sites_type;
    typedef typename model_type::vertex_type vertex_type;
    typedef typename model_type::vertex_state_type vertex_state_type;
    
    typedef Wprob wprob_type;
    typedef std::vector<std::vector<Wprob> > wprob_table_type;
    
    Worm_prob(alps::Parameters const& params,
            lattice_type const& lattice, model_type const& model) :
        lattice(lattice),
        model(model),
        nbstates(model.nbstates())
    {
        worm_weights = !params.value_or_default("NO_WORMWEIGHT", false);
        wormtype = params.value_or_default("WHICH_LOOP_TYPE", "minbounce");        
//        force_symmetry_constraint = params.defined("FORCE_SYMMETRY_CONSTRAINT");
    }
    
    static inline unsigned w_index(unsigned vi, unsigned en_leg, unsigned en_op)
    {
        return en_op + NOPERATORS * (en_leg + 2 * UNIT_SIZE * vi);
    }
    
    void calc_probabilities(wprob_table_type& wprobs, int p0, int p1)
    {
        phys_ops[0] = p0;
        phys_ops[1] = p1;
        
        nvertices = model.nvertices();
        
        unsigned ne = nvertices * 2 * UNIT_SIZE * NOPERATORS;
        wprobs.resize(ne);
        entrances_done.resize(ne);
        
        Matrix_type type;
        
        if (wormtype.compare("heatbath") == 0)
            type = HEATBATH;
        else if (wormtype.compare("locopt") == 0)
            type = LOCOPT;
        else if (wormtype.compare("minbounce") == 0)
            type = MINBOUNCE;
        else
            throw std::runtime_error("WHICH_LOOP_TYPE should be "
                    "minbounce, localopt or heatbath.");
                    
        calc_probabilities(wprobs, type);
    }
private:
    Worm_prob();
    
    static const unsigned NOPERATORS = 2;
    static const unsigned EMPTY;
    
    enum Matrix_type {HEATBATH = 0, LOCOPT, MINBOUNCE};
    
    struct Exit {
        double w;
        unsigned vi;
        unsigned leg;
        unsigned op;
        
        bool operator< (Exit const& y) const { return w < y.w; }
    };
    
    lattice_type const& lattice;
    model_type const& model;
    std::vector<unsigned> const& nbstates;
    
    std::vector<bool> entrances_done;
    
    bool worm_weights;
//    bool force_symmetry_constraint;
    
    std::string wormtype;
    
    int phys_ops[NOPERATORS];
    
    unsigned nvertices;
    
    void calc_probabilities(wprob_table_type& wprobs, Matrix_type type)
    {
        for (unsigned vi = 0; vi < nvertices; ++vi) {
            vertex_type const& vertex = model.vertex(vi);
            
            double w = model.c(vertex.unit_type) - vertex.me;
            if (vertex.diagonal && w < 1e-15)
                // zero weight
                continue;
            
            for (unsigned en_leg = 0; en_leg < 2 * UNIT_SIZE; ++en_leg)
            for (unsigned en_op = 0; en_op < NOPERATORS; ++en_op)
                if (!entrances_done[w_index(vi, en_leg, en_op)])
                    do_entrance(wprobs, type, vertex, en_leg, en_op);
        }
    }
    
    void do_entrance(wprob_table_type& wprobs, Matrix_type type,
        vertex_type const& vertex, unsigned en_leg, unsigned en_op)
    {
        lat_unit_sites_type sites;
        lattice.lat_unit_type2sitei(vertex.unit_type, sites);

        vertex_state_type const& vstate = vertex.state;
        
        unsigned nbstates2 =
            nbstates[lattice.sitei2alps_type(sites[en_leg % UNIT_SIZE])];
        if (!check_state(vstate[en_leg], phys_ops[en_op], nbstates2))
            return;
        
        std::vector<Exit> exits;
        
        for (unsigned ex_leg = 0; ex_leg < 2 * UNIT_SIZE; ++ex_leg)
        for (unsigned ex_op = 0; ex_op < NOPERATORS; ++ex_op) {
            unsigned nbstates2 =
                nbstates[lattice.sitei2alps_type(sites[ex_leg % UNIT_SIZE])];
            if (en_leg != ex_leg) {
                if (!check_state(vstate[ex_leg], phys_ops[ex_op], nbstates2))
                    continue;
            } else {
                if (!check_state(vstate[ex_leg], phys_ops[en_op] + phys_ops[ex_op], nbstates2))
                    continue;
            }
            
            // construct new vertex state
            vertex_state_type nstate;
            for (unsigned i = 0; i < 2 * UNIT_SIZE; ++i)
                nstate[i] = vstate[i];
            nstate[en_leg] += phys_ops[en_op];
            nstate[ex_leg] += phys_ops[ex_op];

            // find new vertex
            unsigned vi2 = model.find_vertex(nstate, vertex.unit_type);
            if (vi2 == model_type::INVALID_VERTEX)
                continue;
                
            vertex_type const& vertex2 = model.vertex(vi2);
            
            // get weight
            double w2;
            if (vertex2.diagonal)
                w2 = model.c(vertex2.unit_type) - vertex2.me;
            else
                // off diagonal matrix elements here should be always positive
                w2 = std::fabs(vertex2.me);
            if (worm_weights) {
                unsigned state = vertex2.state[ex_leg];
                unsigned stype = lattice.sitei2alps_type(sites[ex_leg % UNIT_SIZE]);
                if (ex_op == 1)
                    w2 *= model.lowering_matrix_elements()[stype][state];
                else
                    w2 *= model.raising_matrix_elements()[stype][state];
            }
            if (w2 < 1e-15)
                // zero weight
                continue;
            
            Exit exit = {w2, vi2, ex_leg, ex_op};
            exits.push_back(exit);
        }
        
        if (type == LOCOPT)
            std::sort(exits.begin(), exits.end());
        
        std::vector<double> ws(exits.size());
        for (unsigned i = 0; i < exits.size(); ++i)
            ws[i] = exits[i].w;
            
        std::vector<double> matrix(ws.size() * ws.size());

        switch (type) {
        case HEATBATH:
            calc_heatbath_matrix(ws, matrix);
            break;
        case LOCOPT:
            calc_localopt_matrix(ws, matrix);
            break;
        case MINBOUNCE:
            calc_minbounce_matrix(ws, matrix);
            break;
        }
        
        // map transition matrix to probability tables
        for (unsigned i = 0; i < exits.size(); ++i) {
            unsigned count = 0;
            for (unsigned j = 0; j < exits.size(); ++j) {
                double prob = matrix[exits.size() * i + j];
                if (prob > 1e-14)
                    // nonzero matrix element
                    ++count;
            }
            
            if (count == 0)
                throw std::runtime_error("Invalid transition matrix.");
            
            unsigned wi = w_index(exits[i].vi, exits[i].leg, 1 - exits[i].op);
            
            entrances_done[wi] = true;
            wprobs[wi].resize(count);

            count = 0;
            for (unsigned j = 0; j < exits.size(); ++j) {
                double prob = matrix[exits.size() * i + j];
                if (prob > 1e-14) {
                    // nonzero matrix element
                    Wprob wprob;

                    wprob.prob = prob;
                    wprob.ex_leg = exits[j].leg;
                    wprob.ex_op = exits[j].op;
                    wprob.vertex_index = exits[j].vi;

                    wprobs[wi][count++] = wprob;
                }
            }

            // order wprobs according to probability values
            std::sort(wprobs[wi].begin(), wprobs[wi].end(), wprobs[wi][0]);

            double psum = wprobs[wi][0].prob;
            for (unsigned k = 1; k < wprobs[wi].size(); ++k) {
                wprobs[wi][k].prob += psum;
                psum = wprobs[wi][k].prob;
            }

            if (std::fabs(psum - 1.0) > 1e-12)
                std::cerr << "Sum of probabilities must be equal to 1.";
        }
    }
    
    bool check_state(int state, int phys_op, int nbstates)
    {
        return state + phys_op >= 0 && state + phys_op < nbstates;
    }
        
    void calc_heatbath_matrix(
        std::vector<double> const& ws, std::vector<double>& matrix) const
    {
        double iwsum = 1.0 / std::accumulate(ws.begin(), ws.end(), 0.0);
     
        unsigned size = ws.size() * ws.size();
        for (unsigned i = 0; i < size; ++i)
            matrix[i] = iwsum * ws[i % ws.size()];
    }
    
    void calc_localopt_matrix(
        std::vector<double>& ws, std::vector<double>& matrix) const
    {
        double wsum = std::accumulate(ws.begin(), ws.end(), 0.0);
        
        double mul = 1.0;
        std::vector<double> ys(ws.size());
        
        for (unsigned i = 0; i < ws.size(); ++i) {
            ys[i] = mul;
            if (i < ws.size() - 1) {
                wsum -= ws[i];
                ys[i] = mul * ws[i] / wsum;
                mul -= ys[i];
            }
            
            double iw = 1.0 / ws[i];
            
            for (unsigned j = 0; j < ws.size(); ++j) {
                double& me = matrix[ws.size() * i + j];

                if (j < i)
                    me = ys[j];
                else if (j > i)
                    me = iw * ws[j] * ys[i];
                else if (i < ws.size() - 1)
                    me = 0.0;
                else
                    me = ys[i];
            }
        }
    }
    
    void calc_minbounce_matrix(
        std::vector<double> const& ws, std::vector<double>& matrix) const
    {
        LP_Solve_SSE lp_solver;

        lp_solver.init(ws);
        lp_solver.calc_solution(matrix);
    }
};

template<typename L, typename M> const unsigned Worm_prob<L, M>::EMPTY = 0xffffffff;

#endif
