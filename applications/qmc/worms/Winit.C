/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2005 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>
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

/* $Id$ */

#include "WRun.h"
#include <iomanip>
//#include <alps/osiris/std/set.h>

void WRun::print_copyright(std::ostream& out)
{
  out << "Worm algorithm quantum Monte Carlo simulation v2.0\n"
      << "  copyright (c) 1997-2007 by Simon Trebst <trebst@comp-phys.org>\n"
      << "                          and Matthias Troyer <troyer@comp-phys.org>\n"
      << " for details see the publication:\n"     
      << " A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).\n\n";
}



//- Constructor/setup run -------------------------------------------------

WRun::WRun(const alps::ProcessList& w, const alps::Parameters& myparms,int n)
  : QMCRun<>(w,myparms,n,myparms.defined("NUMBER_OF_PARTICLES")),
    canonical(parms.defined("NUMBER_OF_PARTICLES")),
    corrections_upwards(0),
    corrections_downwards(0),
    preadjustment_done(false),
    adjustment_done(false),
    adjust_parameter(parms.defined("ADJUST") ? static_cast<std::string>(parms["ADJUST"]) : std::string("mu")),
    worm_head(2,wormhead_type(graph(),kinks)),
    stat(4),
    steps(0),
    min_number(parms.value_or_default("MIN_NUMBER",0)),
    max_number(parms.value_or_default("MAX_NUMBER",min_number-1)),
    eta(beta*double(parms.value_or_default("eta", 1.))),
    thermal_sweeps(parms.required_value("THERMALIZATION")),
    skip_measurements(parms.value_or_default("SKIP",1)),
    measurements_done(skip_measurements),
    have_worm(false),
    chain_kappa(parms.value_or_default("CHAIN_KAPPA",false)),
    num_kinks(num_sites()),
    worms_per_kink(parms.value_or_default("WORMS_PER_KINK",1)),
    worms_per_update(1.),
    log_numeric_limits_double(std::log(std::numeric_limits<double>::max())),
    Sign(1),
    nonlocal(parms.value_or_default("NONLOCAL",true)),
    use_1D_stiffness(parms.value_or_default("USE_1D_STIFFNESS",false)),  //@#$br
    chain_number(num_sites()),
    bond_type(alps::get_or_default(alps::bond_type_t(),graph(),0)),
    boundary_crossing(alps::get_or_default(alps::boundary_crossing_t(),graph(),alps::boundary_crossing()))
{

  // no site compressibility measurements yet
  measure_site_compressibility_=false;
  
  if(eta<0)
    boost::throw_exception(std::out_of_range("negative eta is illegal"));
  if(thermal_sweeps<0)
    boost::throw_exception( std::out_of_range("negative thermalization is illegal"));
  
  // kinks
  if(!where.empty()) {
    int vol = num_sites();
    kinks.resize(vol);
    initial_state_.resize(vol,state_type());
  }
  
  // model
  initialize_site_states();
  initialize_hamiltonian();        
  // print_hamiltonian();

  // subintervals for non-local interactions
  subinterval.clear();
  subinterval_valid = false;
  current_head_num = 0;

  create_observables();
}

//- Dumps -----------------------------------------------------------------

void WRun::save(alps::ODump& dump) const
{
  dump << 100 << num_kinks << worms_per_update << initial_state_ << steps << kinks << Sign << last_id_;
}

void WRun::load(alps::IDump& dump)
{
  int version;
  dump >> version >> num_kinks >> worms_per_update;
  if(!where.empty())
    dump >> initial_state_ >> steps >> kinks >> Sign >> last_id_;
} 
