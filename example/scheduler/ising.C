/***************************************************************************
* PALM++/scheduler example
*
* scheduler/ising.C an example Ising model simulation
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
// This file implements the simulation specific classes for a simple
// simulation of a one-dimensional Ising model
//=======================================================================

#include "ising.h"
#include <cmath>

IsingSimulation::IsingSimulation(const alps::ProcessList& where,const alps::Parameters& p,int node)
  : alps::scheduler::MCRun(where,p,node),
    length(static_cast<uint32_t>(parms["L"])),
    beta(1./static_cast<double>(parms["T"])),
    sweeps(0),
    thermalization_sweeps(static_cast<uint32_t>(parms["THERMALIZATION"])),
    total_sweeps(static_cast<uint32_t>(parms["SWEEPS"])),
    spins(length)
{
  // initialize random spin configuration
  for(int i=0;i<length;i++)
    spins[i]=(random_real() <0.5 ? 1 : -1);

  // create measurement objects
  measurements << alps::RealObservable("Energy");
  measurements << alps::RealObservable("Magnetization");
  measurements << alps::RealObservable("Magnetization^2");
  measurements << alps::RealObservable("Magnetization^4");
#ifdef ALPS_HAVE_VALARRAY
  measurements << alps::RealVectorObservable("Correlations");
#endif
}

void IsingSimulation::load(alps::IDump& dump)
{
  dump >> sweeps;
  if(!where.empty()) // skip reading the spins if we are just evaluating
    dump >> spins; 
}

void IsingSimulation::save(alps::ODump& dump) const
{
  dump << sweeps << spins;
}

bool IsingSimulation::change_parameter(const std::string& name, const alps::StringValue& value)
{
  if(name=="SWEEPS")
    total_sweeps=static_cast<uint32_t>(value);
  else if (name=="THERMALIZATION" && !is_thermalized())
    thermalization_sweeps=static_cast<uint32_t>(value);
  else
    return false; // cannot do it
  return true; // could do it
}


bool IsingSimulation::is_thermalized() const
{
  return (sweeps >= thermalization_sweeps);
}

double IsingSimulation::work_done() const
{
  return (is_thermalized() ? (sweeps-thermalization_sweeps)/double(total_sweeps) :0.);
}

void IsingSimulation::dostep()
{
  // increment sweep count
  sweeps++;
  
  // perform updates
  for (int j=0;j<length;j++)  {
      // choose a random site and determine the neighbors
      int i = int(double(length)*random_real());
      int right=(i+1 < length ? i+1 : 0);
      int left=( i-1 < 0 ? length-1 : i-1);

      // calculate change in the weight of the configuration
      double p=exp(2.*beta*spins[i]*(spins[right]+spins[left]));

      // Metropolis updating: accept if random number is smaller than p
      if (p>=1. || random_real() < p)
	spins[i]=-spins[i];
    }
    
  // perform measurements
  double tmag=0;
  double ten=0;
  static std::valarray<double> corr;
  corr.resize(length);
  corr=0.;
  for (int i=0;i<length;i++) {
    int right=(i +1 < length ? i+1 : 0);
    tmag += spins[i];
    ten += -spins[i]*spins[right];
    for (int d=0;d<length;d++)
      corr[d]+=spins[i]*spins[(i+d)%length];
  }
 
  // normalize measurements and add them to the observables
  corr /= double(length);
  ten/=length;
  tmag/=length;
  measurements.get<alps::RealObservable>("Energy") << ten;
  measurements.get<alps::RealObservable>("Magnetization") << tmag;
  measurements.get<alps::RealObservable>("Magnetization^2") << tmag*tmag;
  measurements.get<alps::RealObservable>("Magnetization^4") << tmag*tmag*tmag*tmag;
#ifdef ALPS_HAVE_VALARRAY
  measurements.get<alps::RealVectorObservable>("Correlations") << corr;
#endif
}
