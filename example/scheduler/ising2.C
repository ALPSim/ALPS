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

//=======================================================================
// This file implements the simulation specific classes for a simple
// simulation of a one-dimensional Ising model
//=======================================================================

#include "ising2.h"
#include <alps/alea/detailedbinning.h>
#include <cmath>

using namespace std;
using namespace alps;

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using boost::adjacent_vertices;
using boost::num_vertices;
using boost::vertex;
using boost::vertices;
using boost::edges;
using boost::out_edges;
using boost::source;
using boost::target;
#endif

void IsingSimulation2::print_copyright(std::ostream& out)
{
  out << "Ising simulation example program using the ALPS lattice library\n"
      << "  copyright(c) 1994-2003 by Matthias Troyer <troyer@comp-phys.org>\n\n";
}


IsingSimulation2::IsingSimulation2(const alps::ProcessList& where,const alps::Parameters& p,int node)
  : alps::scheduler::LatticeMCRun<graph_type>(where,p,node),
    beta(1./static_cast<double>(p["T"])),
    sweeps(0),
    thermalization_sweeps(static_cast<alps::uint32_t>(p["THERMALIZATION"])),
    total_sweeps(static_cast<alps::uint32_t>(p["SWEEPS"]))
{
  if (inhomogeneous())
    boost::throw_exception(std::runtime_error("Disordered lattices not supported by the Ising example program.\n"));

  spins.resize(num_sites()); // number of vertices = number of lattice sites
  // initialize random spin configuration
  for(int i=0;i<spins.size();i++)
    spins[i]=(random_real() <0.5 ? 1 : -1);

  // create measurement objects
  measurements << alps::RealObservable("Energy");
  measurements << alps::RealObservable("Magnetization");
}

void IsingSimulation2::load(alps::IDump& dump)
{
  dump >> sweeps;
  if(!where.empty()) // skip reading the spins if we are just evaluating
    dump >> spins; 
}

void IsingSimulation2::save(alps::ODump& dump) const
{
  dump << sweeps << spins;
}

bool IsingSimulation2::change_parameter(const std::string& name, const alps::StringValue& value)
{
  if(name=="SWEEPS")
    total_sweeps=static_cast<alps::uint32_t>(value);
  else if (name=="THERMALIZATION" && !is_thermalized())
    thermalization_sweeps=static_cast<alps::uint32_t>(value);
  else
    return false; // cannot do it
  return true; // could do it
}


bool IsingSimulation2::is_thermalized() const
{
  return (sweeps >= thermalization_sweeps);
}

double IsingSimulation2::work_done() const
{
  return (is_thermalized() ? (sweeps-thermalization_sweeps)/double(total_sweeps) :0.);
}

void IsingSimulation2::dostep()
{  
  // increment sweep count
  sweeps++;
  
  // perform updates
  for (int j=0;j<spins.size();j++)  {
    // choose a random site and determine the neighbors
    site_descriptor s = site(random_int(0,num_sites()-1));

    // Metropolis updates
    neighbor_iterator it,end;
    double de=0;
    for (boost::tie(it,end)=neighbors(s);it!=end; ++it)
      de += spins[s]*spins[*it];
    if (de <0. || random_real() < exp(-2.*beta*de))
      spins[s]=-spins[s];
    }
    
  // perform measurements
  double tmag=0;
  double ten=0;
  site_iterator s,s_end;
  for (boost::tie(s,s_end)=sites(); s!=s_end;++s)
    tmag += spins[*s];

  bond_iterator b,b_end;
  for (boost::tie(b,b_end)=bonds(); b!=b_end;++b)
    ten -= spins[source(*b)]*spins[target(*b)];
  
  measurements["Energy"] << ten/spins.size();
  measurements["Magnetization"] << tmag/spins.size();
}
