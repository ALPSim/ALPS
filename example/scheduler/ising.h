/***************************************************************************
* PALM++/scheduler example
*
* scheduler/ising.h an example Ising model simulation
*
* $Id$
*
* Copyright (C) 1994-2002 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*
* This software is part of the ALPS library, published under the 
* ALPS Library License; you can use, redistribute it and/or modify 
* it under the terms of the License, either version 1 or (at your option) 
* any later version.
*
* You should have received a copy of the ALPS Library License along with 
* the ALPS Library; see the file License.txt. If not, the license is also 
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
**************************************************************************/

//=======================================================================
// This file defines the simulation specific classes for a simple
// simulation of a one-dimensional Ising model
//=======================================================================

#include <alps/scheduler.h>

class IsingSimulation : public alps::scheduler::MCRun
{
public:
  static void print_copyright(std::ostream&);
  
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

typedef alps::scheduler::SimpleMCFactory<IsingSimulation> IsingFactory;
