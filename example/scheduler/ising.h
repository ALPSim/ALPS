/***************************************************************************
* PALM++/scheduler example
*
* scheduler/ising.h an example Ising model simulation
*
* $Id$
*
* Copyright (C) 1994-2002 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

class IsingSimulation : public alps::scheduler::MCRun
{
public:
  IsingSimulation(const alps::ProcessList&,const alps::Parameters&,int);
  void save(alps::ODump&) const;
  void load(alps::IDump&);
  void dostep();
  bool is_thermalized() const;
  double work_done() const;
  bool change_parameter(const std::string& name, const alps::StringValue& value);
  
private:
  uint32_t length;                  // the system size
  double beta;                      // the inverse temperatre
  uint32_t sweeps;                  // the number of sweeps done
  uint32_t thermalization_sweeps;   // the number of sweeps to be done for equilibration
  uint32_t total_sweeps;            // the total number of sweeps to be done after equilibration
  std::vector<int> spins;           // the vector to store the spins
};

typedef alps::scheduler::SimpleFactory<IsingSimulation> IsingFactory;
