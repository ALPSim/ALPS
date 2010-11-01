/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2006 by Sergei Isakov <isakov@itp.phys.ethz.ch>
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

#ifndef __LPSOLVER_H__
#define __LPSOLVER_H__

#include <vector>
#include <stack>

#include "lpkit.h"

enum LP_types {LP_EQUAL = EQ};

class LP_solver {
public:
	class Equation {
	public:
		Equation(unsigned size)
		{
			row = new double[size + 1];
			std::memset(row, 0, sizeof(double) * (size + 1));
		}

		~Equation()
		{
			delete [] row;
		}

		void set_at(unsigned index, double value)
		{
			row[index + 1] = value;
			undo.push(index);
		}
	private:
		Equation();
		Equation(const Equation&);

		lprec *lp;
		double *row;
		std::stack<unsigned> undo;

		friend class LP_solver;

		void undo_row()
		{
			while(!undo.empty()) {
				unsigned i = undo.top();
				row[i + 1] = 0.0;
				undo.pop();
			}
		}

		double* get_row()
		{
			return row;
		}
	};
	
	LP_solver() : lp(NULL), size(0) {}
	
	~LP_solver()
	{
		if (lp != NULL)
			delete_lp(lp);
	}
	
	void create(int size)
	{
		if (lp != NULL)
			delete_lp(lp);
		
		this->size = size;
		
		lp = make_lp(0, size);
		if (lp == NULL)
			throw std::runtime_error("Cannot initialize lpkit.");
			
		set_verbose(lp, IMPORTANT);
	}
	
	bool add_constraint(Equation& eq, unsigned type, double value)
	{
		::add_constraint(lp, eq.get_row(), type, value);
		eq.undo_row();
		return true;
	}
	
	bool set_min_objective(Equation& eq)
	{
		set_obj_fn(lp, eq.get_row());
		set_minim(lp);
		eq.undo_row();
		return true;
	}
	
	bool set_max_objective(Equation& eq)
	{
		set_obj_fn(lp, eq.get_row());
		set_maxim(lp);
		eq.undo_row();
		return true;
	}
	
	bool set_low_bound(unsigned index, double value)
	{
		set_lowbo(lp, index + 1, value);
		return true;
	}
	
	bool set_up_bound(unsigned index, double value)
	{
		set_upbo(lp, index + 1, value);
		return true;
	}
	
	bool solve()
	{
		if (::solve(lp) != 0)
			return false;
		else
			return true;
	}
	
	double get_solution(std::vector<double>& sol)
	{
		sol.resize(size);
		get_variables(lp, &sol[0]);

		return get_objective(lp);
	}
private:
	LP_solver(const LP_solver&);
	
	lprec *lp;
	unsigned size;
};

#endif
