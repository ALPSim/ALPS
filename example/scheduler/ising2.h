/***************************************************************************
* PALM++/scheduler example
*
* scheduler/ising.h an example Ising model simulation
*
* $Id$
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
**************************************************************************/

//=======================================================================
// This file defines the simulation specific classes for a simple
// simulation of a one-dimensional Ising model
//=======================================================================

#include <alps/scheduler.h>
#include <alps/lattice.h>

typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS,
			      // vertex property
			      boost::property<alps::vertex_type_t,int>,
			      // edge property
                              boost::property<alps::edge_type_t,int>
			      > graph_type;

class IsingSimulation2 : public alps::scheduler::LatticeMCRun<graph_type>
{
public:
  static void print_copyright(std::ostream&);

  IsingSimulation2(const alps::ProcessList&,const alps::Parameters&,int);
  void save(alps::ODump&) const;
  void load(alps::IDump&);
  void dostep();
  bool is_thermalized() const;
  double work_done() const;
  bool change_parameter(const std::string& name, const alps::StringValue& value);
  
private:
  double beta;                      // the inverse temperature
  uint32_t sweeps;                  // the number of sweeps done
  uint32_t thermalization_sweeps;   // the number of sweeps to be done for equilibration
  uint32_t total_sweeps;            // the total number of sweeps to be done after equilibration
  std::vector<int> spins;           // the vector to store the spins
};
