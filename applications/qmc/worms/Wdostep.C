/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2004 by Matthias Troyer <troyer@comp-phys.org>,
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

using namespace alps;
//#define TIMINGS

double WRun::work_done() const
{
  return (is_thermalized() ? (steps-thermal_sweeps)/double(parms.required_value("SWEEPS")) :0.);
}  

void WRun::start()
{
  green.resize(1+num_sites());
  green=0.;
  measurements_done=skip_measurements;
}

void WRun::dostep()
{
#ifdef TIMINGS
  double tt;
#endif
#ifdef TIMINGS
  tt=-dclock();
#endif
  stat=0;
  if (is_thermalized())
    for (int i=0;i<worms_per_update;++i)
      make_worm();
  else {
    int worm_num=0;
    for (long long length=0;length<=worms_per_kink*num_kinks;++worm_num) {
      length+=make_worm();
    }
    worms_per_update=0.99*worms_per_update+0.01*worm_num;
  }

#ifdef TIMINGS
  tt += dclock();
  std::cerr << "Worm time: " << tt << " seconds.\n";
  tt=-dclock();
#endif
#ifdef CHECK_OFTEN
    check_spins();
#ifdef TIMINGS
  tt += dclock();
  std::cerr << "Check time: " << tt << " seconds\n";
  tt=-dclock();
#endif
#endif

if (canonical&&!adjustment_done&&steps>25) {
  adjustment();
}

if (canonical) {
  if (static_cast<int>(parms["NUMBER_OF_PARTICLES"])==get_particle_number()) 
    measure();
}
else 
  measure();

  steps++;
#ifdef TIMINGS
  tt += dclock();
  std::cerr << "Meas time: " << tt << " seconds\n";
#endif
}   // WRun::dostep

bool WRun::is_thermalized() const
{
  return (steps >= thermal_sweeps);
} 


