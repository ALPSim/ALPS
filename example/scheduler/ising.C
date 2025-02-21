/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2009 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

/* $Id$ */

//=======================================================================
// This file implements the simulation specific classes for a simple
// simulation of a one-dimensional Ising model
//=======================================================================

#include "ising.h"
#include <alps/alea.h>
#include <cmath>

void IsingSimulation::print_copyright(std::ostream& out)
{
  out << "Ising simulation example program using the ALPS Monte Carlo library\n"
      << "  copyright(c) 1994-2006 by Matthias Troyer <troyer@comp-phys.org>\n\n";
}

IsingSimulation::IsingSimulation(const alps::ProcessList& where,const alps::Parameters& p,int node)
  : alps::scheduler::MCRun(where,p,node),
    length(static_cast<alps::uint32_t>(p["L"])),
    beta(1./static_cast<double>(p["T"])),
    sweeps(0),
    thermalization_sweeps(static_cast<alps::uint32_t>(p["THERMALIZATION"])),
    total_sweeps(static_cast<alps::uint32_t>(p["SWEEPS"])),
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
  measurements << alps::RealVectorObservable("Correlations");
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
    total_sweeps=static_cast<alps::uint32_t>(value);
  else if (name=="THERMALIZATION" && !is_thermalized())
    thermalization_sweeps=static_cast<alps::uint32_t>(value);
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
  measurements["Energy"] << ten;
  measurements["Magnetization"] << tmag;
  measurements["Magnetization^2"] << tmag*tmag;
  measurements["Magnetization^4"] << tmag*tmag*tmag*tmag;
  measurements["Correlations"] << corr;
}
