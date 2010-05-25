/****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Applications
*
* Copyright (C) 2006-2010 by Lode Pollet <lpollet@physics.harvard.edu>,
*                            Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
*
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

/* $Id: lowa.measurement.h 3520 2010-03-03 16:55:00Z tamama $ */


/*
 *
 * 1) Code modification      -- partly done
 * 2) Replacing raw pointers -- not done yet
 *
 */


#ifndef ALPS_APPLICATIONS_LOWA_DOSTEP_H
#define ALPS_APPLICATIONS_LOWA_DOSTEP_H


namespace alps {
namespace applications {


void Lowa::dostep() {
  /*
   * Lode's comment: Lowa::dostep() is the core unit. It generates a new diagonal configuration.
   *
   * Tama's comment:
   *   1) worm_from_insertion_to_delete() generates a new closed worldline configuration, where diagonal statistics, eg. particle number ... , can be taken.
   *   2) worm_from_insertion_to_delete() returns true if we can now do saving.
   */

  bool is_now_okay_to_save = worm_from_insertion_to_deletion();
  if (is_now_okay_to_save)
  {
    std::cout << "\nSaving...\t";
    MCstep_total += MCstep;
    MCstep        = 0;

    update_obs();
    save_lowa();

    std::cout << "\t...done\n";
  }
  else {
    ++counter_MEASURE;
    ++counter_MEASURE_GREEN;
    ++counter_TEST;


    if (counter_MEASURE >= Nmeasure) {

      if ((Ncan < 0) || (get_Npart() == Ncan)) {
        ++sweeps;
        time (&times2);
        
        std::cout << "Measuring in real space...\tN = " << get_Nmeaspart() << " ( " << get_Npart() << " ),\t Total elapsed MCstep : " << MCstep_total+MCstep << "\t worm insert ratio : " << worm_insertion_acceptance_ratio() << "\t Average nr MCstep per second  : " << (MCstep_total+MCstep - MCold) / ( times2 - times1) << " Sweeps : " << sweeps << "\t";
        
        MCold = MCstep_total + MCstep;
        times1 = times2;

        update_obs();
        if (is_thermalized())  {  take_diagonal_measurements();  }

        std::cout << "...done\n";

      }

      counter_MEASURE = 0;

    }


    if (counter_MEASURE_GREEN >= Nmeasure_green) {

      if (measure_time_series_density_matrix)   {  ++sweeps_green;  }

      std::cout << "Measuring in momentum space...\t Sweeps (green) : " << sweeps_green << "\t";            

      update_off_diag_obs();
      if (is_thermalized())  {  take_offdiagonal_measurements();  }
      if (measure_time_series_density_matrix)  {  reset_av_dnsmat();  }

      std::cout << "...done\n";

      counter_MEASURE_GREEN = 0;

    }


    if (counter_TEST >= Ntest) {
       std::cout << "\n Testing...\t";
       if (!test_conf()) {std::cout << "Something is wrong here... Press exit the program...\n"; std::cin.get();};
       std::cout << "\t...OK\n";
       counter_TEST = 0;
    }
  }
}


}
}

#endif
