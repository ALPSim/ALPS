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

#include "lattice.h"
#include "model.h"
#include "lpsolver.h"

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
		force_symmetry_constraint = params.defined("FORCE_SYMMETRY_CONSTRAINT");
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
		unsigned nt = ne * 2 * UNIT_SIZE * NOPERATORS * NOPERATORS;

		wprobs.resize(nt);
		entrances_done.resize(ne);
		trans_table.resize(nt);
		
		for (unsigned i = 0; i < nt; ++i)
			trans_table[i].id = EMPTY;
		
		build_trans_table();
		
		if (wormtype.compare("heatbath") == 0) {
			calc_heatbath_probabilities();
		} else if (wormtype.compare("locopt") == 0) {
			calc_localopt_probabilities();
		} else if (wormtype.compare("minbounce") == 0) {
			calc_minbounce_probabilities();
		} else
			throw std::runtime_error("WHICH_LOOP_TYPE should be "
					"minbounce, localopt or heatbath.");
		
		build_prob_table(wprobs);
	}
private:
	Worm_prob();
	
	static const unsigned NOPERATORS = 2;
	static const unsigned EMPTY;
	
	struct Transition {
		double prob;
		double me;               // matrix element of the old vertex
		unsigned reverse;        // index of the reverse transition table entry
		unsigned vertex_index;   // new vertex index
		unsigned id;
	};
	
	struct Entrance {
		unsigned vi;
		unsigned leg;
		unsigned op;
	};
	
	struct wi_s {
		double w;
		unsigned i;
		bool operator< (wi_s const& y) const { return w < y.w; }
	};
	
	lattice_type const& lattice;
	model_type const& model;
	std::vector<unsigned> const& nbstates;
	
	std::vector<bool> entrances_done;
	std::vector<Transition> trans_table;
	std::vector<std::vector<unsigned> > connected_trans;
	
	bool worm_weights;
	bool force_symmetry_constraint;
	
	std::string wormtype;
	
	int phys_ops[NOPERATORS];
	
	unsigned id;
	unsigned nvertices;
	
	bool check_state(int state, int phys_op, int nbstates)
	{
		return state + phys_op >= 0 && state + phys_op < nbstates;
	}
	
	void build_trans_table()
	{
		Entrance en;
		
		id = 0;
		
		for (en.vi = 0; en.vi < nvertices; ++en.vi) {
			vertex_type const& vertex = model.vertex(en.vi);
			
			lat_unit_sites_type sites;
			lattice.lat_unit_type2sitei(vertex.unit_type, sites);
		
			for (en.leg = 0; en.leg < 2 * UNIT_SIZE; ++en.leg)
			for (en.op = 0; en.op < NOPERATORS; ++en.op) {
				unsigned nbstates2 =
					nbstates[lattice.sitei2alps_type(sites[en.leg % UNIT_SIZE])];
				if (!check_state(vertex.state[en.leg], phys_ops[en.op], nbstates2))
					continue;
						
				std::vector<unsigned> connected;	
				do_entrance(en, connected);
				
				if (connected.size() != 0)
					connected_trans.push_back(connected);			
			}
		}
	}
	
	void do_entrance(Entrance const& en, std::vector<unsigned>& connected)
	{
		unsigned wi = w_index(en.vi, en.leg, en.op);
		
		if (entrances_done[wi])
			return;
			
		std::stack<Entrance> stack;
		
		int phys_en_op = phys_ops[en.op];

		vertex_state_type nstate;
		vertex_type const& vertex = model.vertex(en.vi);
		vertex_state_type const& vstate = vertex.state;
		
		lat_unit_sites_type sites;
		lattice.lat_unit_type2sitei(vertex.unit_type, sites);
		
		for (unsigned ex_leg = 0; ex_leg < 2 * UNIT_SIZE; ++ex_leg)
		for (unsigned ex_op = 0; ex_op < NOPERATORS; ++ex_op) {
			int phys_ex_op = phys_ops[ex_op];
			
			unsigned nbstates2 =
				nbstates[lattice.sitei2alps_type(sites[ex_leg % UNIT_SIZE])];
			if (en.leg != ex_leg) {
				if (!check_state(vstate[ex_leg], phys_ex_op, nbstates2))
					continue;
			} else {
				if (!check_state(vstate[ex_leg], phys_en_op + phys_ex_op, nbstates2))
					continue;
			}
			
			// construct new vertex state
			for (unsigned i = 0; i < 2 * UNIT_SIZE; ++i)
				nstate[i] = vstate[i];
			nstate[en.leg] += phys_en_op;
			nstate[ex_leg] += phys_ex_op;
			
			// find new vertex
			unsigned vi2 = model.find_vertex(nstate, vertex.unit_type);
			if (vi2 == model_type::INVALID_VERTEX)
				continue;
				
			unsigned ti = t_index(
							en.vi, en.leg, en.op, ex_leg, ex_op);
							
			connected.push_back(ti);
			
			Transition& trans = trans_table[ti];
			trans.prob = 0.0;
			
			if (vertex.diagonal)
				trans.me = model.c(vertex.unit_type) - vertex.me;
			else
				// off diagonal matrix elements here should be always positive
				trans.me = std::fabs(vertex.me);
								
			if (worm_weights) {
				unsigned state = vstate[en.leg];
				unsigned stype = lattice.sitei2alps_type(sites[en.leg % UNIT_SIZE]);
				if (en.op == 0)
					trans.me *= model.lowering_matrix_elements()[stype][state];
				else
					trans.me *= model.raising_matrix_elements()[stype][state];
			}
			
			trans.id = id++;
			trans.reverse = t_index(
							vi2, ex_leg, 1 - ex_op, en.leg, 1 - en.op);
			trans.vertex_index = vi2;
			
			Entrance en2 = {vi2, ex_leg, 1 - ex_op};
			stack.push(en2);
		}
		
		entrances_done[wi] = true;
		
		while (!stack.empty()) {
			do_entrance(stack.top(), connected);
			stack.pop();
		}
	}
	
	void calc_heatbath_probabilities()
	{
		for (unsigned vi = 0; vi < nvertices; ++vi)
		for (unsigned en_leg = 0; en_leg < 2 * UNIT_SIZE; ++en_leg)
		for (unsigned en_op = 0; en_op < NOPERATORS; ++en_op) {

			double psum = 0.0;
			for (unsigned ex_leg = 0; ex_leg < 2 * UNIT_SIZE; ++ex_leg)
			for (unsigned ex_op = 0; ex_op < NOPERATORS; ++ex_op) {
				unsigned ti = t_index(vi, en_leg, en_op, ex_leg, ex_op);
				if (trans_table[ti].id != EMPTY)
					psum += trans_table[trans_table[ti].reverse].me;
			}
						
			if (psum == 0.0)
				continue;

			for (unsigned ex_leg = 0; ex_leg < 2 * UNIT_SIZE; ++ex_leg)
			for (unsigned ex_op = 0; ex_op < NOPERATORS; ++ex_op) {
				unsigned ti = t_index(vi, en_leg, en_op, ex_leg, ex_op);
				if (trans_table[ti].id != EMPTY) {
					unsigned ti2 = trans_table[ti].reverse;
					trans_table[ti].prob = trans_table[ti2].me / psum;
				}
			}
		}
	}
	
	void calc_localopt_probabilities()
	{
		// fix me: do this in different way; use connected_trans
		
		// this works only with UNIT_SIZE == 2 and offdiagonal bond operators
		// that conserve quantum numbers
		
		for (unsigned vi = 0; vi < nvertices; ++vi)
		for (unsigned en_leg = 0; en_leg < 2 * 2; ++en_leg)
		for (unsigned en_op = 0; en_op < NOPERATORS; ++en_op) {
			
			double w[4] = {0.0, 0.0, 0.0, 0.0};
			double ti_matrix[4][4];
			for (unsigned i = 0; i < 4; ++i)
				for (unsigned j = 0; j < 4; ++j)
					ti_matrix[i][j] = EMPTY;
			
			bool valid = false;
			for (unsigned ex_leg = 0; ex_leg < 2 * 2; ++ex_leg)
			for (unsigned ex_op = 0; ex_op < NOPERATORS; ++ex_op) {
				unsigned ti = t_index(vi, en_leg, en_op, ex_leg, ex_op);
				
				if (trans_table[ti].id != EMPTY) {
					valid = true;

					w[ex_leg] = trans_table[trans_table[ti].reverse].me;
					ti_matrix[en_leg][ex_leg] = ti;
						
					if (en_leg == ex_leg) continue;
						
					unsigned rvi = trans_table[ti].vertex_index;
					for (unsigned j = 0; j < 4; ++j)
						ti_matrix[ex_leg][j] =
								reverse_index(rvi, ex_leg, 1 - ex_op, j);
						
					continue;
				}
			}
			
			if (!valid)
				continue;
						
			double matrix[4][4];
			calc_localopt_matrix(w, matrix);
			
			for (unsigned i = 0; i < 4; ++i)
				for (unsigned j = 0; j < 4; ++j) {
					unsigned ti = ti_matrix[i][j];
					if (ti != EMPTY)
						trans_table[ti].prob = matrix[i][j];
				}
		}
	}
	
	unsigned reverse_index(
		unsigned vi, unsigned en_leg, unsigned en_op, unsigned ex_leg)
	{
		for (unsigned ex_op = 0; ex_op < NOPERATORS; ++ex_op) {
			unsigned ti = t_index(vi, en_leg, en_op, ex_leg, ex_op);
			if (trans_table[ti].id != EMPTY && trans_table[ti].me != 0.0)
				return ti;
		}

		return EMPTY;
	}
	
	void calc_localopt_matrix(double (&w)[4], double (&matrix)[4][4]) const
	{	
		std::vector<wi_s> wi(4);
		for (unsigned i = 0; i < 4; ++i) {
			wi[i].w = w[i];
			wi[i].i = i;
		}
		std::sort(wi.begin(), wi.end());
		
		double tmp_matrix[4][4];
		
		double x0 = 1.0 / (wi[1].w + wi[2].w + wi[3].w);
		double y0 = x0 * wi[0].w;
		double x1 = (1.0 - y0) / (wi[2].w + wi[3].w);
		double y1 = x1 * wi[1].w;
		double x2 = (1.0 - y0 - y1) / wi[3].w;
		double y2 = x2 * wi[2].w;
		double y3 = (1.0 - y0 - y1 - y2);
		
	    tmp_matrix[0][0] = 0.0;
	    tmp_matrix[0][1] = wi[1].w * x0;
	    tmp_matrix[0][2] = wi[2].w * x0;
	    tmp_matrix[0][3] = wi[3].w * x0;

	    tmp_matrix[1][0] = y0;
	    tmp_matrix[1][1] = 0.0;
	    tmp_matrix[1][2] = wi[2].w * x1;
	    tmp_matrix[1][3] = wi[3].w * x1;

	    tmp_matrix[2][0] = y0;
		tmp_matrix[2][1] = y1;
	    tmp_matrix[2][2] = 0.;
	    tmp_matrix[2][3] = wi[3].w * x2;

	    tmp_matrix[3][0] = y0;
	    tmp_matrix[3][1] = y1;
	    tmp_matrix[3][2] = y2;
	    tmp_matrix[3][3] = y3;
		
		for (unsigned i = 0; i < 4; ++i)
			for (unsigned j = 0; j < 4; ++j)
				matrix[wi[i].i][wi[j].i] = tmp_matrix[i][j];
	}
	
	void calc_minbounce_probabilities()
	{
		unsigned ne = nvertices * 2 * UNIT_SIZE * NOPERATORS;
		unsigned nt = ne * 2 * UNIT_SIZE * NOPERATORS * NOPERATORS;
		
		std::vector<unsigned> indices;
		indices.resize(nt, EMPTY);
		
		LP_solver lpsolver;
		 
		for (unsigned i = 0; i < connected_trans.size(); ++i) {
			unsigned size = connected_trans[i].size();
			
			assert(size > 0);
			
			LP_solver::Equation constr1(size);
			LP_solver::Equation constr2(size);
			LP_solver::Equation constr3(size);
			LP_solver::Equation obj(size);

			lpsolver.create(size);
			
			for (unsigned j = 0; j < size; ++j)
				indices[connected_trans[i][j]] = j;

			for (unsigned j = 0; j < size; ++j) {
				unsigned ti = connected_trans[i][j];
				unsigned vi, en_leg, en_op, ex_leg, ex_op;
				unpack_t_index(ti, vi, en_leg, en_op, ex_leg, ex_op);
				
				for (unsigned ex_leg = 0; ex_leg < 2 * UNIT_SIZE; ++ex_leg)
				for (unsigned ex_op = 0; ex_op < NOPERATORS; ++ex_op) {
					unsigned ti = t_index(vi, en_leg, en_op, ex_leg, ex_op);

					Transition& trans = trans_table[ti];

					if (trans.id == EMPTY)
						continue;

					// bounds on probabilities
					lpsolver.set_low_bound(indices[ti], 0.0);
					lpsolver.set_up_bound(indices[ti], 1.0);

					// sum of exit probabilities must be equal to 1
					constr1.set_at(indices[ti], 1.0);

					if (en_leg == ex_leg) {
						assert(ti == trans.reverse);

						// objective --- minimize the sum of bounce probabilities
						obj.set_at(indices[ti], 1.0);
						continue;
					}

					if (trans.id > trans_table[trans.reverse].id)
						// do not do it twice
						continue;

					// detailed balance
					constr2.set_at(indices[ti], trans.me);
					constr2.set_at(indices[trans.reverse],
											-trans_table[trans.reverse].me);
					lpsolver.add_constraint(constr2, LP_EQUAL, 0.0);
				}
				
				// sum of exit probabilities must be equal to 1
				lpsolver.add_constraint(constr1, LP_EQUAL, 1.0);
				
				if (force_symmetry_constraint)
					for (unsigned k = 0; k < size; ++k) {
						unsigned ti2 = connected_trans[i][k];
						unsigned vi2, en_leg2, en_op2, ex_leg2, ex_op2;
						unpack_t_index(ti2, vi2, en_leg2, en_op2, ex_leg2, ex_op2);
						
						if (trans_table[ti].vertex_index == trans_table[ti2].vertex_index
							&& en_leg != en_leg2 && ex_leg == ex_leg2
							&& ex_leg != en_leg && ex_leg != en_leg2
							&& trans_table[ti].me == trans_table[ti2].me) {
							constr3.set_at(indices[ti], 1.0);
							constr3.set_at(indices[ti2], -1.0);
							lpsolver.add_constraint(constr3, LP_EQUAL, 0.0);
						}
					}
			}
			
			lpsolver.set_min_objective(obj);

			if (!lpsolver.solve())
				std::cout << "lp_solve found no solution\n";

			std::vector<double> solution;
			lpsolver.get_solution(solution);

			for (unsigned j = 0; j < size; ++j)
				trans_table[connected_trans[i][j]].prob = solution[j];
		}
	}
	
	void build_prob_table(wprob_table_type& wprobs)
	{
		for (unsigned vi = 0; vi < nvertices; ++vi)
		for (unsigned en_leg = 0; en_leg < 2 * UNIT_SIZE; ++en_leg)
		for (unsigned en_op = 0; en_op < NOPERATORS; ++en_op) {
		
			unsigned count = 0;
			for (unsigned ex_leg = 0; ex_leg < 2 * UNIT_SIZE; ++ex_leg)
			for (unsigned ex_op = 0; ex_op < NOPERATORS; ++ex_op) {
				unsigned ti = t_index(vi, en_leg, en_op, ex_leg, ex_op);
				Transition& trans = trans_table[ti];
					
				if (trans.prob > 0.0)
					++count;
			}
				
			if (count == 0)
				// no valid updates
				continue;
					
			unsigned wi = w_index(vi, en_leg, en_op);
			wprobs[wi].resize(count);

			count = 0;
			for (unsigned ex_leg = 0; ex_leg < 2 * UNIT_SIZE; ++ex_leg)
			for (unsigned ex_op = 0; ex_op < NOPERATORS; ++ex_op) {
				unsigned ti = t_index(vi, en_leg, en_op, ex_leg, ex_op);
				Transition& trans = trans_table[ti];

				if (trans.prob > 0.0) {
					Wprob wprob;
					
					wprob.prob = trans.prob;
					wprob.ex_leg = ex_leg;
					wprob.ex_op = ex_op;
					wprob.vertex_index = trans.vertex_index;

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
			
			if (std::fabs(psum - 1.0) > 1e-13)
				std::cerr << "Sum of probabilities must be equal to 1.\n";
		}
	}
	
	inline unsigned t_index(unsigned vi, unsigned en_leg, unsigned en_op,
			unsigned ex_leg, unsigned ex_op)
	{
		return ex_op + NOPERATORS * (ex_leg + 2 * UNIT_SIZE * (en_op +
			NOPERATORS * (en_leg + 2 * UNIT_SIZE * vi)));
	}

	inline void unpack_t_index(unsigned index, unsigned& vi, unsigned& en_leg,
			unsigned& en_op, unsigned& ex_leg, unsigned& ex_op)
	{
		ex_op = index % NOPERATORS;
		index = index / NOPERATORS;

		ex_leg = index % (2 * UNIT_SIZE);
		index = index / (2 * UNIT_SIZE);

		en_op = index % NOPERATORS;
		index = index / NOPERATORS;

		en_leg = index % (2 * UNIT_SIZE);
		vi = index / (2 * UNIT_SIZE);
	}
};

template<typename L, typename M> const unsigned Worm_prob<L, M>::EMPTY = 0xffffffff;

#endif
