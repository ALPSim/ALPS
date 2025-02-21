/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2005 by Fabien Alet <alet@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>,
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

#include <alps/osiris/comm.h>
#include "SSE.hpp"

void SSE::print_copyright(std::ostream& out)
{
  out << "Quantum Monte simulations using the generalized directed loop algorithm v. 1.1\n"
      << "  available from http://alps.comp-phys.org/\n"
      << "  copyright (c) 2001-2005 by Fabien Alet <alet@comp-phys.org>,\n"
      << "                             Synge Todo <wistaria@comp-phys.org>,\n"
      << "                             and Matthias Troyer <troyer@comp-phys.org>\n"
      << "  see F. Alet, S. Wessel and M. Troyer, Phys. Rev. E 71, 036706 (2005) for details.\n\n";
}

SSE::SSE(const alps::ProcessList& w, const alps::Parameters& myparms,int n)
  : QMCRun<>(w,myparms,n),
    // Simulation parameters
  epsilon(parms.value_or_default("EPSILON", 0.)),
  cutoff_L(parms.value_or_default("INITIAL_CUTOFF",10)),
  number_of_worms_per_sweep(parms.value_or_default("NUMBER_OF_WORMS_PER_SWEEP",1)),
  nb_steps(parms.defined("SWEEPS") ? alps::uint64_t(parms.required_value("SWEEPS")) 
    :  (parms.defined("MCS") ? alps::uint64_t(parms["MCS"]) : alps::uint64_t(parms["Steps"]))),
  nb_thermalisation_steps(parms.defined("THERMALIZATION") ? alps::uint32_t(parms["THERMALIZATION"]) 
    :  (parms.defined("thermalization") ? alps::uint32_t(parms["thermalization"]) : nb_steps/10)),
  each_measurement(parms.value_or_default("SKIP",1)),
  steps_done_total(0),
  measurements_done(0),
  WHICH_LOOP_TYPE(parms.value_or_default("WHICH_LOOP_TYPE","minbounce")),
  NO_WORMWEIGHT(parms.value_or_default("NO_WORMWEIGHT", false)),
  nb_worms_thermalization(0),
  count_worms_thermalization(0),
  site_state(num_sites(),0),
  operator_string(cutoff_L),
  op_indices(cutoff_L),
  current_number_of_non_identity(0),
  boundary_crossing(alps::get_or_default(alps::boundary_crossing_t(),graph(),alps::boundary_crossing())),
  worm_abort(parms.value_or_default("WORM_ABORT", 0.))
{
  if (worm_abort)
    measure_green_function_=false;
  // Initialize simulation
  create_observables();
}

// Check if simulation is finished
double SSE::work_done() const
{
  return (is_thermalized() ? (steps_done_total-nb_thermalisation_steps)/double(nb_steps) :0.);
}

// It is desirable to be able to increase the number of SWEEPS
bool SSE::change_parameter(const std::string& name, const StringValue& value) {
  alps::uint64_t new_sweeps = 0;
  if(name=="SWEEPS")
    new_sweeps = alps::uint64_t(value);
  if(name=="MCS")
    new_sweeps = alps::uint64_t(value);
  if(name=="Steps")
    new_sweeps = alps::uint64_t(value);
  // Is it sensoble to do it ?
  if(new_sweeps > 0) {
    nb_steps = new_sweeps;
    return true;
   }
  // Otherwise we cannot do it
  return false;
}


// Do a MC step - and do measurements (if thermalized)
void SSE::dostep()
{
  do_update();
  steps_done_total++;
  
  //  if (!(steps_done_total%1000)) cout << steps_done_total << endl;
  if (is_thermalized() && (++measurements_done==each_measurement)) {
    measurements_done=0;
    do_measurements();
  }
}

// Check if simulation os thermalized
bool SSE::is_thermalized() const
{
  return (steps_done_total >= nb_thermalisation_steps); 
}

// save simulation
void SSE::save(alps::ODump& dump) const
{
  dump <<  steps_done_total << site_state << operator_string << cutoff_L <<  current_number_of_non_identity << number_of_worms_per_sweep << nb_worms_thermalization << count_worms_thermalization;
}

// load simulation
void SSE::load(alps::IDump& dump)
{
  if(dump.version() >= 302)
    dump >>  steps_done_total;
  else {
    // data type have changed from 32 to 64 Bit between version 301 and 302
    alps::uint32_t steps_done_total_tmp;
    dump >>  steps_done_total_tmp;
    // perform the conversion which may be necessary
    steps_done_total = steps_done_total_tmp;
   }
  if(!where.empty()) {
    dump >>  site_state >> operator_string >>  cutoff_L >> current_number_of_non_identity >> number_of_worms_per_sweep >>  nb_worms_thermalization >> count_worms_thermalization;
    op_indices.resize(cutoff_L);
  } else
    measurements.compact();
}

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
   return alps::scheduler::start(argc,argv,alps::scheduler::SimpleMCFactory<SSE>());
#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (std::exception& exc) {
    std::cerr << exc.what() << "\n";
      alps::comm_exit(true);
      return -1;
    }
  catch (...) {
    std::cerr << "Fatal Error: Unknown Exception!\n";
    return -2;
  }
#endif
}
