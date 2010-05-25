/*****************************************************************************
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

/* $Id: lowa.propagation.h 3520 2010-03-03 16:55:00Z tamama $ */


/*
 *
 * 1) Code modification      -- partly done
 * 2) Replacing raw pointers -- not done yet
 *
 */


/*
 *
 * This file describes how the worm head propagates right after the insertion of the worm pair and right before its deletion.
 *
 */ 


#ifndef ALPS_APPLICATIONS_LOWA_PROPAGATION_H
#define ALPS_APPLICATIONS_LOWA_PROPAGATION_H



void Lowa::worm_propagation(Lowa::time_type& p)
{

  time_type El, Er, dE;
  time_type segm_tau;

  // we distinguish between moving to the right and moving to the left
  if (is_worm_moving_forward) {

    // find out the energies before and after the worm
    if (is_worm_rising) {
      dE = U(worm_head.to()) * (site_it[worm_head.to()]->before() - 1) + _mu_eff[worm_head.to()];
      dE = ( dE > 0 ? _Eoffset : -dE + _Eoffset);
    }
    else {
      dE =  -U(worm_head.to()) * (site_it[worm_head.to()]->before())- _mu_eff[worm_head.to()];
      dE = (dE > 0 ? _Eoffset : -dE + _Eoffset);
    }
#ifdef DEBUGMODE
    std::cout << "\n moving right p = " << p << "\tdE " << dE << "\tdens " << site_it[worm_head.to()]->before() << "\t_mu " << _mu_eff[worm_head.to()] << "\tworm_it " << *worm_it << "\tsite_it " << *site_it[worm_head.to()];
    std::cout << "\n Worm : " << worm_head.to() << "\t" << worm_head.time();
    if (worm_it != site_it[worm_head.to()]) std::cout << "\nWRONG assignment of worm_it";
#endif

    if (is_kink_encountered_forward(p, dE, segm_tau) ) {
#ifdef DEBUGMODE
      std::cout << "\n moving right p = " << p << "\tdE " << dE << "\tsegm_tau " << segm_tau;
#endif
      if ((is_time0_cyclically_between_time1_and_time2(worm_tail.time(), worm_head.time(), worm_it->time())) && (worm_tail.time()!=(worm_it->time() ))) {
#ifdef DEBUGMODE
    std::cout << "\nMoving upwards, measuring density matrix " << worm_tail.time() << "\t" << worm_head.time() << "\t" << worm_it->time();
#endif

    // measure density matrix
    worm_head.set_time(worm_tail.time());
    //dns[dist(worm_head.to(), worm_tail.to())] += _state_enclosed_within_wormpair_at_creation;
    if (worm_head.to() != worm_tail.to()) {
#ifdef MEASURE_DENSITY_MATRIX
          double cosphi=1.; for (int k=0; k < dim; k++) cosphi *= std::cos(phase[k]* (distsq(worm_head.to(),k) - distsq(worm_tail.to(),k)));
#ifdef HARDCORE
#ifdef TRAPPEDSYSTEM
          av_dnsmat[dist(worm_head.to(), worm_tail.to())] += 0.25*_state_enclosed_within_wormpair_at_creation * cosphi;
      av_dnsmat_inf[dist(worm_head.to(), worm_tail.to())] += 0.25*_state_enclosed_within_wormpair_at_creation;
#else
          av_dnsmat[dist(worm_head.to(), worm_tail.to())] += 0.5*_state_enclosed_within_wormpair_at_creation * cosphi;
      av_dnsmat_inf[dist(worm_head.to(), worm_tail.to())] += 0.5*_state_enclosed_within_wormpair_at_creation;
#endif
#else
#ifdef TRAPPEDSYSTEM       
          av_dnsmat[dist(worm_head.to(), worm_tail.to())] += 0.5*_state_enclosed_within_wormpair_at_creation * cosphi;
          av_dnsmat_inf[dist(worm_head.to(), worm_tail.to())] += 0.5*_state_enclosed_within_wormpair_at_creation;
#else
          av_dnsmat[dist(worm_head.to(), worm_tail.to())] += 1.*_state_enclosed_within_wormpair_at_creation * cosphi;
      av_dnsmat_inf[dist(worm_head.to(), worm_tail.to())] += 1.*_state_enclosed_within_wormpair_at_creation;
#endif 
#endif
#endif
        }
    p = 0.;
    return;
      }

      worm_head.set_time(worm_it->time());  //worm jumps to next kink

      if (worm_it->from() == worm_it->to()) // the time of the worm head is the same as the time of the worm tail
      {
#ifdef DEBUGMODE
        std::cout << "\nMoving upwards, time head = time tail" << "\t" << worm_head << "\tworm_it " << *worm_it;
#endif
        if (worm_head.time() == worm_tail.time()) { // worm bites in tail (same sites)
#ifdef DEBUGMODE
      std::cout << "\n ...removing worm";
#endif
      //dns[0] = worm_it->after();  
#ifdef MEASURE_DENSITY_MATRIX
#ifdef HARDCORE
#ifdef TRAPPEDSYSTEM
          av_dnsmat[0] += worm_it->after();
      av_dnsmat_inf[0] += worm_it->after();
#else
          av_dnsmat[0] += worm_it->after();
      av_dnsmat_inf[0] += worm_it->after();
#endif
#else
#ifdef TRAPPEDSYSTEM           
          av_dnsmat[0] += worm_it->after();
      av_dnsmat_inf[0] += worm_it->after();
#else
          av_dnsmat[0] += worm_it->after();
      av_dnsmat_inf[0] += worm_it->after();
#endif 
#endif
#endif


      reset_assoc_upon_delete(worm_head.to(), site_it[worm_head.to()]);
          site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].erase(site_it[worm_head.to()]);
          if (site_it[worm_head.to()] == _kinks_on_lattice[worm_head.to()].end()) site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].begin();
      worm_it = site_it[worm_head.to()];
          is_worm_diagonal = true;

          p = 0.;
          return;
        }
        else {  // worm head passes dummy
#ifdef DEBUGMODE
      std::cout << "\n ...passing dummy" << worm_head.to() << "\t" << worm_tail.to() << "\trising : "<< is_worm_rising;
#endif

      if (is_worm_rising) {
        site_it[worm_head.to()]->set_before(site_it[worm_head.to()]->before()-1);
        site_it[worm_head.to()]->set_after( site_it[worm_head.to()]->after() -1);
      }
      else {
        site_it[worm_head.to()]->set_before(site_it[worm_head.to()]->before()+1);
        site_it[worm_head.to()]->set_after( site_it[worm_head.to()]->after() +1);
      }
#ifdef DEBUGMODE
      std::cout << "\n ...passed dummy" << *worm_it << "\t" << *site_it[worm_head.to()];
#endif
          move_kinks_iterator_cyclically_forward(site_it[worm_head.to()],worm_head.to());
      worm_it = site_it[worm_head.to()];
#ifdef DEBUGMODE
      std::cout << "\n ...passed dummy" << *worm_it;
#endif
          p -= segm_tau * dE;
          return;
        } // else ... worm passes tail
      } // time of worm head equals time of worm tail  
      else if (worm_it->from() == worm_head.to())
      {
#ifdef DEBUGMODE
        std::cout << "\nMoving upwards, next.from = head.to";
#endif
        if (!is_worm_rising) {  //pass the interaction with probability 1
      // no iterators change, worm remains on same site
      site_it[worm_head.to()]->set_before(site_it[worm_head.to()]->before()+1);
      site_it[worm_head.to()]->set_after( site_it[worm_head.to()]->after() +1);
          move_kinks_iterator_cyclically_forward(site_it[worm_head.to()],worm_head.to());
      worm_it = site_it[worm_head.to()];
          p -= segm_tau*dE;
          return;
        }
        else // try to delete/relink interaction
        {
      for (site_type j = 0; j < zcoord; j++) {
        kinks_iterator_type it = site_it[worm_it->from()]->assoc(j);
        if (worm_it->time() == it->time()) {
          site_it[worm_it->to()] = it;
          break;
        }
      }
      //std::cout << "\n setup_local_transition_weight_upon_deletion iterator for to" << site_it[worm_it->to()]->time() << "\t" << &*site_it[worm_it->to()];
          setup_local_transition_weight_upon_deletion(worm_it->to()); // this should also set the iterators right
          dim_type i3 = select_index_via_heatbatch_algorithm();
          if (i3 == zcoord) { //annihilate interaction
#ifdef DEBUGMODE
            std::cout << "\n  annihilate " << i3;
#endif
        site_type ib = worm_it->to();

        // delete interaction on the from site
        reset_assoc_upon_delete(worm_it->from(), site_it[worm_it->from()]);
        reset_assoc_upon_delete(worm_it->to(), site_it[worm_it->to()]);
        site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == _kinks_on_lattice[worm_head.to()].end()) site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].begin();

        // delete interaction on the to site
            worm_head.set_to(ib);
            site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == _kinks_on_lattice[worm_head.to()].end()) site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].begin();
        worm_it = site_it[worm_head.to()]; // worm jumps
            p -= segm_tau * dE;
        --_nrvertex;
        //print_configuration(std::cout);
            return;
          }
          else //relink
          {
#ifdef DEBUGMODE
            std::cout << "\n relink or bounce " << i3 << "\t" << nb(worm_it->to(),i3);
#endif
        // define the new kink
        site_type ito = worm_it->to();
            site_type ib = nb(ito,i3);


        site_it[ito]->set_from(ib);

        // delete the old kink on the from site
        reset_assoc_upon_delete(worm_it->from(), site_it[worm_it->from()]);
        site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == _kinks_on_lattice[worm_head.to()].end()) site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].begin();

        // insert new kink
        int n0 = site_it[ib]->before();
        kink_type new_elem(zcoord, n0+1,n0, ib, ito, worm_head.time() );
        site_it[ib] = _kinks_on_lattice[ib].insert(site_it[ib], new_elem);
        //std::cout << "\nafter insert, before reset_assoc_upon!!! " << site_it[ib]->time() << "\t" << site_it[ito]->time();
        reset_assoc_upon_insert(ib, site_it[ib]);

        worm_head.set_to(ib); // worm jumps
            move_kinks_iterator_cyclically_backward(site_it[ib],ib);  
        worm_it = site_it[ib];
        is_worm_moving_forward = !is_worm_moving_forward;
        is_worm_rising         = !is_worm_rising;
            p -= segm_tau * dE;
            return;
          }
        } // ... else delete/relink interaction
      }
      else if (worm_it->to() == worm_head.to()) {
#ifdef DEBUGMODE
    std::cout << "\nMoving upwards next.to = head.to";
#endif
        if (is_worm_rising) { // pass interaction with probability 1
      // no iterators change, worm remains on same site
      site_it[worm_head.to()]->set_before(site_it[worm_head.to()]->before()-1);
      site_it[worm_head.to()]->set_after( site_it[worm_head.to()]->after() -1);
          move_kinks_iterator_cyclically_forward(site_it[worm_head.to()],worm_head.to());
      worm_it = site_it[worm_head.to()];
          p -= segm_tau * dE;
          return;
        }
        else { // try to delete/relink interaction
      for (site_type j = 0; j < zcoord; j++) {
        kinks_iterator_type it = site_it[worm_it->to()]->assoc(j);
        if (worm_it->time() == it->time()) {
          site_it[worm_it->from()] = it;
          break;
        }
      }
      setup_local_transition_weight_upon_deletion(worm_it->from());
          dim_type i3 = select_index_via_heatbatch_algorithm();

          if (i3 == zcoord) { //annihilate interaction

            site_type ib = worm_it->from();

#ifdef DEBUGMODE
        std::cout << "\n Annihilate " << i3 << std::endl;
#endif

        // delete kink on the to site
        reset_assoc_upon_delete(worm_it->to(), site_it[worm_it->to()]);
        reset_assoc_upon_delete(worm_it->from(), site_it[worm_it->from()]);
        site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == _kinks_on_lattice[worm_head.to()].end()) site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].begin();

        // delete kink on the from site
            worm_head.set_to(ib);
            site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == _kinks_on_lattice[worm_head.to()].end()) site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].begin();
        worm_it = site_it[worm_head.to()]; // worm jumps

            p -= segm_tau * dE;
        --_nrvertex;
            return;
          }
          else {// relink the interaction
#ifdef DEBUGMODE
        std::cout << "\n Relink : " << i3 << std::endl;
#endif
        // define the new kink
        site_type ifrom = worm_it->from();
            site_type ib = nb(ifrom,i3);

        site_it[ifrom]->set_to(ib);

        // delete the old kink on the to site
        reset_assoc_upon_delete(worm_it->to(), site_it[worm_it->to()]);
        site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == _kinks_on_lattice[worm_head.to()].end()) site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].begin();

        // insert new kink
        int n0 = site_it[ib]->before();
        kink_type new_elem(zcoord, n0-1,n0, ifrom, ib, worm_head.time() );
        site_it[ib] = _kinks_on_lattice[ib].insert(site_it[ib], new_elem);
        reset_assoc_upon_insert(ib, site_it[ib]);

        worm_head.set_to(ib); // worm jumps
            move_kinks_iterator_cyclically_backward(site_it[ib],ib);  
        worm_it = site_it[ib];
        is_worm_moving_forward = !is_worm_moving_forward;
        is_worm_rising         = !is_worm_rising;
            p -= segm_tau*dE;
        return;
          }
        }  //....try to delete/relink interaction
      }
      else // worm passes the interaction (no sites in common), but diagonal energy might change
      {
#ifdef DEBUGMODE
        std::cout << "\nMoving upwards, worm passes remote interaction. This should not occur at all!!";
    //print_configuration(std::cout);
    std::cin.get();
#endif
      }
    } // is_kink_encountered_forward

    else { // no interaction reached, try to insert new interaction

#ifdef DEBUGMODE
      std::cout << "\n moving right p = " << p << "\t dE " << dE << "\t segm_tau" << segm_tau;
      std::cout << "\nMoving to the right, inserting interaction " << worm_head.to() << "\t" << is_worm_rising ;
#endif
      time_type newtime = worm_head.time() + p/dE;
      while (newtime > beta) {
         newtime -= beta;
      }

      if (  is_time0_cyclically_between_time1_and_time2(worm_tail.time(), worm_head.time(), newtime) ) { 
#ifdef DEBUGMODE
    std::cout << "\nMoving to the right, not inserting interaction, but measuring the density matrix ";
#endif

    // measure density matrix

    worm_head.set_time(worm_tail.time());
    //dns[dist(worm_head.to(), worm_tail.to())] += _state_enclosed_within_wormpair_at_creation;
    if (worm_head.to() != worm_tail.to()) {
#ifdef MEASURE_DENSITY_MATRIX
          double cosphi=1.; for (int k=0; k < dim; k++) cosphi *= std::cos(phase[k]* (distsq(worm_head.to(),k) - distsq(worm_tail.to(),k)));
#ifdef HARDCORE
#ifdef TRAPPEDSYSTEM
          av_dnsmat[dist(worm_head.to(), worm_tail.to())] += 0.25*_state_enclosed_within_wormpair_at_creation * cosphi;
      av_dnsmat_inf[dist(worm_head.to(), worm_tail.to())] += 0.25*_state_enclosed_within_wormpair_at_creation;
#else
          av_dnsmat[dist(worm_head.to(), worm_tail.to())] += 0.5*_state_enclosed_within_wormpair_at_creation * cosphi;
      av_dnsmat_inf[dist(worm_head.to(), worm_tail.to())] += 0.5*_state_enclosed_within_wormpair_at_creation;
#endif
#else
#ifdef TRAPPEDSYSTEM           
          av_dnsmat[dist(worm_head.to(), worm_tail.to())] += 0.5*_state_enclosed_within_wormpair_at_creation * cosphi;
      av_dnsmat_inf[dist(worm_head.to(), worm_tail.to())] += 0.5*_state_enclosed_within_wormpair_at_creation;
#else
          av_dnsmat[dist(worm_head.to(), worm_tail.to())] += 1.*_state_enclosed_within_wormpair_at_creation * cosphi;
      av_dnsmat_inf[dist(worm_head.to(), worm_tail.to())] += 1.*_state_enclosed_within_wormpair_at_creation;
#endif 
#endif
#endif
    }
    p = 0.;
    //print_configuration(std::cout);
    return;
      }

      setup_local_transition_weight_upon_insertion(dE, worm_head.to(), newtime); // this should also set the iterators right...
      worm_head.set_time(newtime);

      dim_type i3 = select_index_via_heatbatch_algorithm();
#ifdef DEBUGMODE
      std::cout << "\n ..." << i3 << "\t" << nb(worm_head.to(),i3);
#endif
      if ( i3 < zcoord)
      {
        site_type ib = nb(worm_head.to(),i3);

        kink_type new_elem(zcoord);
        new_elem.set_time(worm_head.time());
        if (is_worm_rising) {
      new_elem.set_to(worm_head.to());
      new_elem.set_from(ib);
      new_elem.set_after(site_it[worm_head.to()]->before());
      new_elem.set_before(site_it[worm_head.to()]->before()-1);
      site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].insert(site_it[worm_head.to()], new_elem);
      new_elem.set_after(site_it[ib]->before()-1);
      new_elem.set_before(site_it[ib]->before());
      site_it[ib] = _kinks_on_lattice[ib].insert(site_it[ib], new_elem);
        }
    else {
          new_elem.set_to(ib);
          new_elem.set_from(worm_head.to());
      new_elem.set_after(site_it[worm_head.to()]->before());
      new_elem.set_before(site_it[worm_head.to()]->before()+1);
      site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].insert(site_it[worm_head.to()], new_elem);
      new_elem.set_after(site_it[ib]->before()+1);
      new_elem.set_before(site_it[ib]->before());
      site_it[ib] = _kinks_on_lattice[ib].insert(site_it[ib], new_elem);
    }
    reset_assoc_upon_insert(worm_head.to(), site_it[worm_head.to()]);
    reset_assoc_upon_insert(ib, site_it[ib]);
        move_kinks_iterator_cyclically_forward(site_it[ib],ib);

    // worm jumps
        worm_head.set_to(ib);
    worm_it = site_it[ib];
        p = 0.;
    ++_nrvertex;
        return;
      }
      else // bounce back
      {
#ifdef DEBUGMODE
        std::cout << "\nInserting interaction failed, bounstd::cing back, moving upwards";
#endif

        move_kinks_iterator_cyclically_backward(site_it[worm_head.to()],worm_head.to());
    worm_it = site_it[worm_head.to()];
    is_worm_moving_forward = !is_worm_moving_forward;
    is_worm_rising         = !is_worm_rising;
        p = 0.;
        return;
      }
    }
   } // is_worm_moving_forward 
   
//*************************************************************************************
  else { // !is_worm_moving_forward; this is mirror symmetry of the case is_worm_moving_forward. it's all about getting the +1 and -1 right
    if (is_worm_rising) {
      dE = U(worm_head.to()) * (site_it[worm_head.to()]->after() - 1) + _mu_eff[worm_head.to()];
      dE = (dE > 0 ? _Eoffset : -dE + _Eoffset);
    }
    else {
      dE = -U(worm_head.to()) * (site_it[worm_head.to()]->after()) - _mu_eff[worm_head.to()];
      dE = (dE > 0 ? _Eoffset : -dE + _Eoffset);
    }
#ifdef DEBUGMODE
    std::cout << "\n moving left p = " << p << "\tdE " << dE << "\tdens " << site_it[worm_head.to()]->after() << "\t_mu " << _mu_eff[worm_head.to()] << "\tworm_it " << *worm_it << "\tsite_it " << *site_it[worm_head.to()];
    std::cout << "\n Worm : " << worm_head.to() << "\t" << worm_head.time();
    //char chr; std::cin>> chr;
    if (worm_it != site_it[worm_head.to()]) std::cout << "\nWRONG assignment of worm_it";
#endif

    if (is_kink_encountered_backward(p, dE, segm_tau)) {

      if ( ( is_time0_cyclically_between_time1_and_time2(worm_tail.time(), worm_it->time(), worm_head.time() )) && ( worm_tail.time()  !=  worm_it->time() ) && (worm_tail.time() != worm_head.time()) ) { 
    // measure density matrix
#ifdef DEBUGMODE
    std::cout << "\n measuring density matrix";
#endif

    worm_head.set_time(worm_tail.time());
    //dns[dist(worm_head.to(), worm_tail.to())] += _state_enclosed_within_wormpair_at_creation;
    if (worm_head.to() != worm_tail.to()) {
#ifdef MEASURE_DENSITY_MATRIX
          double cosphi=1.; for (int k=0; k < dim; k++) cosphi *= std::cos(phase[k]* (distsq(worm_head.to(),k) - distsq(worm_tail.to(),k)));
#ifdef HARDCORE
#ifdef TRAPPEDSYSTEM
          av_dnsmat[dist(worm_head.to(), worm_tail.to())] += 0.25*_state_enclosed_within_wormpair_at_creation * cosphi;
      av_dnsmat_inf[dist(worm_head.to(), worm_tail.to())] += 0.25*_state_enclosed_within_wormpair_at_creation;
#else
          av_dnsmat[dist(worm_head.to(), worm_tail.to())] += 0.5*_state_enclosed_within_wormpair_at_creation * cosphi;
      av_dnsmat_inf[dist(worm_head.to(), worm_tail.to())] += 0.5*_state_enclosed_within_wormpair_at_creation;
#endif
#else
#ifdef TRAPPEDSYSTEM           
      av_dnsmat[dist(worm_head.to(), worm_tail.to())] += 0.5*_state_enclosed_within_wormpair_at_creation * cosphi;
      av_dnsmat_inf[dist(worm_head.to(), worm_tail.to())] += 0.5*_state_enclosed_within_wormpair_at_creation;
#else
          av_dnsmat[dist(worm_head.to(), worm_tail.to())] += 1.*_state_enclosed_within_wormpair_at_creation * cosphi;
      av_dnsmat_inf[dist(worm_head.to(), worm_tail.to())] += 1.*_state_enclosed_within_wormpair_at_creation;
#endif 
#endif
#endif
    }
    p = 0.;
    return;
      }

      worm_head.set_time(worm_it->time());  //worm jumps to next kink
      if (worm_it->from() == worm_it->to()) { // the time of the worm head is the same as the time of the worm tail
#ifdef DEBUGMODE
        std::cout << "\nMoving downwards, time head equals time tail";
#endif
        if (worm_head.time() == worm_tail.time()) // worm bites in tail (same sites)
        {
#ifdef DEBUGMODE
          std::cout << "\n  remove tail in moving downwards";
#endif
          //dns[0] = worm_it->before(); 
#ifdef MEASURE_DENSITY_MATRIX
#ifdef HARDCORE
#ifdef TRAPPEDSYSTEM
          av_dnsmat[0] += worm_it->before();
      av_dnsmat_inf[0] += worm_it->before();
#else
          av_dnsmat[0] += worm_it->before();
      av_dnsmat_inf[0] += worm_it->before();
#endif
#else
#ifdef TRAPPEDSYSTEM           
          av_dnsmat[0] += worm_it->before();
      av_dnsmat_inf[0] += worm_it->before();
#else
          av_dnsmat[0] += worm_it->before();
      av_dnsmat_inf[0] += worm_it->before();
#endif 
#endif
#endif

      reset_assoc_upon_delete(worm_head.to(), site_it[worm_head.to()]);

      site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].erase(site_it[worm_head.to()]);
      if (site_it[worm_head.to()] == _kinks_on_lattice[worm_head.to()].end()) site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].begin();
      worm_it = site_it[worm_head.to()];
          is_worm_moving_forward = true;
          is_worm_diagonal       = true;
          p = 0.;
          return;
        }
    else {  // worm head passes tail on other site, update green function (different sites)
#ifdef DEBUGMODE
      std::cout << "\n  passing dummy in moving downwards";
#endif

      if (is_worm_rising) {
        site_it[worm_head.to()]->set_before(site_it[worm_head.to()]->before()-1);
        site_it[worm_head.to()]->set_after(site_it[worm_head.to()]->after() -1);
      }
      else {
        site_it[worm_head.to()]->set_before(site_it[worm_head.to()]->before()+1);
        site_it[worm_head.to()]->set_after(site_it[worm_head.to()]->after() +1);
      }
          move_kinks_iterator_cyclically_backward(site_it[worm_head.to()],worm_head.to());      
      worm_it = site_it[worm_head.to()];
      p -= segm_tau * dE;
          //print_configuration(std::cout);
      return;
    } // else ... worm passes tail
      } // time of worm head equals time of worm tail
      else if (worm_it->from() == worm_head.to()) {
#ifdef DEBUGMODE
    std::cout << "\nMoving downwards, prev.from = head.to " << is_worm_rising;
#endif
        if (is_worm_rising) { // pass interaction with probability 1
      site_it[worm_head.to()]->set_after( site_it[worm_head.to()]->after()  - 1);
      site_it[worm_head.to()]->set_before(site_it[worm_head.to()]->before() - 1);
          move_kinks_iterator_cyclically_backward(site_it[worm_head.to()],worm_head.to());
      worm_it = site_it[worm_head.to()];
          p -= segm_tau * dE;
          return;
        }
        else { // try to delete/relink interaction
      for (site_type j = 0; j < zcoord; j++) {  
        kinks_iterator_type it = site_it[worm_it->from()]->assoc(j);
        if (worm_it->time() == it->time()) {
          site_it[worm_it->to()] = it;
          break;
        }
      }
      setup_local_transition_weight_upon_deletion(worm_it->to());        // This should also put the iterators right
          dim_type i3 = select_index_via_heatbatch_algorithm();

          if (i3 == zcoord) {  
        //annihilate interaction
#ifdef DEBUGMODE
            std::cout << "\n  annihilate " << i3;
#endif

        site_type ib = worm_it->to();


        // delete interaction on the from site
        reset_assoc_upon_delete(worm_it->from(), site_it[worm_it->from()]);
        reset_assoc_upon_delete(worm_it->to(), site_it[worm_it->to()]);
        site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == _kinks_on_lattice[worm_head.to()].end()) site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].begin();
            move_kinks_iterator_cyclically_backward(site_it[worm_head.to()],worm_head.to());
        // delete interaction on the to site
            worm_head.set_to(ib);
            site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == _kinks_on_lattice[worm_head.to()].end()) site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].begin();
            move_kinks_iterator_cyclically_backward(site_it[worm_head.to()],worm_head.to());        
        worm_it = site_it[worm_head.to()]; // worm jumps
            p -= segm_tau * dE;
        --_nrvertex;
        //print_configuration(std::cout);
            return;
          }
          else { //relink
#ifdef DEBUGMODE
        std::cout << "\n relink or bounce " << i3 << "\t" << nb(worm_it->to(),i3);
#endif
        // define the new kink
        site_type ito = worm_it->to();
            site_type ib  = nb(ito,i3);

        site_it[ito]->set_from(ib);

        // delete the old kink on the from site
        reset_assoc_upon_delete(worm_it->from(), site_it[worm_it->from()]);
        site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == _kinks_on_lattice[worm_head.to()].end()) site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].begin();

        // insert new kink
        int n0 = site_it[ib]->before(); 
        kink_type new_elem(zcoord,n0,n0-1, ib, ito, worm_head.time() );
        site_it[ib] = _kinks_on_lattice[ib].insert(site_it[ib], new_elem);
        reset_assoc_upon_insert(ib, site_it[ib]);

        worm_head.set_to(ib); // worm jumps
            move_kinks_iterator_cyclically_forward(site_it[ib],ib);
        worm_it = site_it[ib];
        is_worm_moving_forward = !is_worm_moving_forward;
        is_worm_rising         = !is_worm_rising;
            p -= segm_tau * dE;
            return;
          }
        }  //....try to delete/relink interaction
      }
      else if (worm_it->to() == worm_head.to()) {
#ifdef DEBUGMODE
        std::cout << "\nMoving downwards, prev.to = head.to";
#endif
        if (!is_worm_rising) {  //pass the interaction with probability 1
      site_it[worm_head.to()]->set_after( site_it[worm_head.to()]->after()  + 1);
      site_it[worm_head.to()]->set_before(site_it[worm_head.to()]->before() + 1);
          move_kinks_iterator_cyclically_backward(site_it[worm_head.to()],worm_head.to());
      worm_it = site_it[worm_head.to()];
          p -= segm_tau * dE;
          return;
        }
        else { // try to delete/relink interaction
      for (site_type j = 0; j < zcoord; j++) {  
        kinks_iterator_type it = site_it[worm_it->to()]->assoc(j);
        if (worm_it->time() == it->time()) {
          site_it[worm_it->from()] = it;
          break;
        }
      }
      //std::cout << "\n setup_local_transition_weight_upon_deletion iterator for from" << site_it[worm_it->from()]->time() << "\t" << &*site_it[worm_it->from()];
      setup_local_transition_weight_upon_deletion(worm_it->from());
          dim_type i3 = select_index_via_heatbatch_algorithm();
          if (i3 == zcoord) { //annihilate interaction
#ifdef DEBUGMODE
            std::cout << "\n  annihilate " << i3;
#endif
        site_type ib = worm_it->from();

        // delete interaction on the from site
        reset_assoc_upon_delete(worm_it->to(),   site_it[worm_it->to()]  );
        reset_assoc_upon_delete(worm_it->from(), site_it[worm_it->from()]);
        site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == _kinks_on_lattice[worm_head.to()].end()) site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].begin();
            move_kinks_iterator_cyclically_backward(site_it[worm_head.to()],worm_head.to());
        // delete interaction on the to site
            worm_head.set_to(ib);
            site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == _kinks_on_lattice[worm_head.to()].end()) site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].begin();
            move_kinks_iterator_cyclically_backward(site_it[worm_head.to()],worm_head.to());
        worm_it = site_it[worm_head.to()];
            p -= segm_tau * dE;
        --_nrvertex;
            return;
          }
          else { // relink the interaction
#ifdef DEBUGMODE
        std::cout << "\n relink or bounce " << i3 << "\t" << nb(worm_it->from(),i3);
#endif
        // define the new kink
        site_type ifrom = worm_it->from();
            site_type ib  = nb(ifrom,i3);

        site_it[ifrom]->set_to(ib);

        // delete the old kink on the to site
        reset_assoc_upon_delete(worm_it->to(), site_it[worm_it->to()]);
        site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == _kinks_on_lattice[worm_head.to()].end()) site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].begin();

        // insert new kink
        int n0 = site_it[ib]->before();
        kink_type new_elem(zcoord,n0,n0+1, ifrom, ib, worm_head.time() );
        site_it[ib] = _kinks_on_lattice[ib].insert(site_it[ib], new_elem);
        reset_assoc_upon_insert(ib, site_it[ib]);

        worm_head.set_to(ib); // worm jumps
            move_kinks_iterator_cyclically_forward(site_it[ib],ib);
        worm_it = site_it[ib];
        is_worm_moving_forward = !is_worm_moving_forward;
        is_worm_rising         = !is_worm_rising;
            p -= segm_tau * dE;
            return;
          }
        } // ... else delete/relink interaction
      }
      else {// worm passes the interaction (no sites in common), but diagonal energy might change
#ifdef DEBUGMODE
        std::cout << "\nMoving downwards, worm passes remote interaction. This should not occur at all!!";
    char ch; std::cin >> ch;
#endif
      }
    } // is_kink_encountered_backward
    else { // no interaction reached, try to insert new interaction
      time_type newtime = worm_head.time() - p/dE;
      while (newtime <= 0.) {
         newtime += beta;
      }
#ifdef DEBUGMODE
      std::cout << "\nMoving to the left, inserting interaction " << worm_head.to() << "\t" << is_worm_rising << "\t " << newtime;
#endif
      if (  (worm_tail.time() != worm_head.time() ) && (is_time0_cyclically_between_time1_and_time2(worm_tail.time(), newtime, worm_head.time() )) ) { 
    // measure density matrix
#ifdef DEBUGMODE
    std::cout << "\nMoving to the left, measuring density matrix ";
#endif

    worm_head.set_time(worm_tail.time());
    //dns[dist(worm_head.to(), worm_tail.to())] += _state_enclosed_within_wormpair_at_creation;
    if (worm_head.to() != worm_tail.to()) {
#ifdef MEASURE_DENSITY_MATRIX
          double cosphi=1.; for (int k=0; k < dim; k++) cosphi *= std::cos(phase[k]* (distsq(worm_head.to(),k) - distsq(worm_tail.to(),k)));
#ifdef HARDCORE
#ifdef TRAPPEDSYSTEM
          av_dnsmat[dist(worm_head.to(), worm_tail.to())] += 0.25*_state_enclosed_within_wormpair_at_creation * cosphi;
      av_dnsmat_inf[dist(worm_head.to(), worm_tail.to())] += 0.25*_state_enclosed_within_wormpair_at_creation;
#else
          av_dnsmat[dist(worm_head.to(), worm_tail.to())] += 0.5*_state_enclosed_within_wormpair_at_creation * cosphi;
      av_dnsmat_inf[dist(worm_head.to(), worm_tail.to())] += 0.5*_state_enclosed_within_wormpair_at_creation;
#endif
#else
#ifdef TRAPPEDSYSTEM           
          av_dnsmat[dist(worm_head.to(), worm_tail.to())] += 0.5*_state_enclosed_within_wormpair_at_creation * cosphi;
      av_dnsmat_inf[dist(worm_head.to(), worm_tail.to())] += 0.5*_state_enclosed_within_wormpair_at_creation;
#else
          av_dnsmat[dist(worm_head.to(), worm_tail.to())] += 1.*_state_enclosed_within_wormpair_at_creation * cosphi;
      av_dnsmat_inf[dist(worm_head.to(), worm_tail.to())] += 1.*_state_enclosed_within_wormpair_at_creation;
#endif 
#endif
#endif
    }

    p = 0.;
    return;
      }

      move_kinks_iterator_cyclically_forward(site_it[worm_head.to()],worm_head.to());  // then the iterator is the same as for right_insert
      worm_head.set_time(newtime);
      setup_local_transition_weight_upon_insertion(dE, worm_head.to(), newtime); // this should also set the iterators right...

      dim_type i3 = select_index_via_heatbatch_algorithm();

      if ( i3 < zcoord)
      {
        site_type ib = nb(worm_head.to(),i3);
#ifdef DEBUGMODE
    std::cout << "\n ..." << i3 << "\t" << ib << "\t" << worm_head.to();
#endif
        kink_type new_elem(zcoord);
        new_elem.set_time(worm_head.time());
        if (is_worm_rising) {
      new_elem.set_to(ib);
          new_elem.set_from(worm_head.to());
      new_elem.set_after(site_it[worm_head.to()]->before() );
      new_elem.set_before(site_it[worm_head.to()]->before() + 1  );
      site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].insert(site_it[worm_head.to()], new_elem);

      new_elem.set_after(site_it[ib]->before());
      new_elem.set_before(site_it[ib]->before()-1);
      site_it[ib] = _kinks_on_lattice[ib].insert(site_it[ib], new_elem);
        }
        else {
          new_elem.set_to(worm_head.to());
      new_elem.set_from(ib);
      new_elem.set_after(  site_it[worm_head.to()]->before()  );
      new_elem.set_before( site_it[worm_head.to()]->before() - 1    );
      site_it[worm_head.to()] = _kinks_on_lattice[worm_head.to()].insert(site_it[worm_head.to()], new_elem);
      
      new_elem.set_after( site_it[ib]->before()  );
      new_elem.set_before(site_it[ib]->before()+1);
      site_it[ib] = _kinks_on_lattice[ib].insert(site_it[ib], new_elem);
    }
        reset_assoc_upon_insert(worm_head.to(), site_it[worm_head.to()]);
    reset_assoc_upon_insert(ib, site_it[ib]);
        move_kinks_iterator_cyclically_backward(site_it[worm_head.to()],worm_head.to());   // restore 
        move_kinks_iterator_cyclically_backward(site_it[ib],ib);  

    
    // worm jumps
        worm_head.set_to(ib);
    worm_it = site_it[ib];
        p = 0.;
    ++_nrvertex;
    //print_configuration(std::cout);
        return;
      }
      else // bounce back
      {
#ifdef DEBUGMODE
        std::cout << "\nInserting interaction failed, bounstd::cing back, moving upwards";
#endif
    worm_it = site_it[worm_head.to()];
    is_worm_moving_forward = !is_worm_moving_forward;
    is_worm_rising         = !is_worm_rising;
        p = 0.;
        return;
      }
    }
  } // is_worm_moving_forward < 0  
}



inline Lowa::dim_type Lowa::select_index_via_heatbatch_algorithm()
{
  obs_type q = (random_real()*total_local_transition_weight);
  for (dim_type index=0; index <= zcoord; ++index)  {  if (q < cummulative_local_transition_weight[index])  {  return index;  }  }
}



inline void Lowa::setup_local_transition_weight_upon_deletion(Lowa::site_type const cur_site)
// Note: cur_site is here NOT the worm site, it is the other bond site
// Note: only works for 3D currently
{
  fock_basis_type cur_dens = (is_worm_moving_forward ? site_it[cur_site]->before() : site_it[cur_site]->after());
  /*
   * find the densities on the neigboring sites
   * site_it[cur_site] should point to the kink which we want to delete
   * all neighboring iterators should be put right as well
   */
  for (dim_type i=0; i < zcoord; ++i)
  {
    site_type s = nb(cur_site,i);
    site_it[s] = site_it[cur_site]->assoc(i);
    nbs[i] = site_it[s]->before();
    if (s == worm_head.to()) {
      nbs[i] = (is_worm_moving_forward ? site_it[worm_head.to()]->after() : site_it[worm_head.to()]->before() );
      // when moving upwards, we need the density just below the worm; which is the same as after the interaction.
    }
    //std::cout << "\n calc_del_w : " << cur_site << "\t "<< &*site_it[cur_site] << "\t" << s << "\t" << site_it[s]->time();
  }

  if (is_worm_rising)
  { 
    local_transition_weight[0] = (nbs[0] == _statemax ? 0 : (tx_plus(cur_site) * (nbs[0]+1)));
    local_transition_weight[1] = (nbs[1] == _statemax ? 0 : (ty_plus(cur_site) * (nbs[1]+1)));
    local_transition_weight[2] = (nbs[2] == _statemax ? 0 : (tz_plus(cur_site) * (nbs[2]+1)));
    local_transition_weight[3] = (nbs[3] == _statemax ? 0 : (tx_minus(cur_site) * (nbs[3]+1)));
    local_transition_weight[4] = (nbs[4] == _statemax ? 0 : (ty_minus(cur_site) * (nbs[4]+1)));
    local_transition_weight[5] = (nbs[5] == _statemax ? 0 : (tz_minus(cur_site) * (nbs[5]+1)));

    obs_type Ebis = U(cur_site) * cur_dens + _mu_eff[cur_site];
    local_transition_weight[6] = (Ebis > 0 ? (_Eoffset+Ebis) : _Eoffset);
  }
  else
  {
    local_transition_weight[0] = tx_plus(cur_site) * nbs[0];
    local_transition_weight[1] = ty_plus(cur_site) * nbs[1];
    local_transition_weight[2] = tz_plus(cur_site) * nbs[2];
    local_transition_weight[3] = tx_minus(cur_site) * nbs[3];
    local_transition_weight[4] = ty_minus(cur_site) * nbs[4];
    local_transition_weight[5] = tz_minus(cur_site) * nbs[5];

    obs_type Ebis = -U(cur_site) * (cur_dens-1) - _mu_eff[cur_site];
    local_transition_weight[zcoord] = (Ebis > 0 ? (_Eoffset + Ebis) : _Eoffset);
  }

  // heatbath local_transition_weights
  // this seems enough... please change to locally optimal ones if you are not happy with autocorr
  total_local_transition_weight = std::accumulate(&local_transition_weight[0], &local_transition_weight[zcoord+1], 0.);  
  cummulative_local_transition_weight[0] = local_transition_weight[0];
  for (dim_type i=1; i <= zcoord; ++i)
  {
    cummulative_local_transition_weight[i] = local_transition_weight[i] + cummulative_local_transition_weight[i-1];
  }
}


inline void Lowa::setup_local_transition_weight_upon_insertion(Lowa::obs_type const E, Lowa::site_type const cur_site, Lowa::time_type const instime)
// Note: only works for 3D currently
{
  local_transition_weight[6] = E; //  when bouncing
  
  for (dim_type i = 0; i < zcoord; ++i)
  {
    site_type s = nb(cur_site,i);
    site_it[s] = site_it[cur_site]->assoc(i);
    kinks_iterator_type previt = site_it[s];
    //std::cout << "\n calc_ins_local_transition_weight : " << cur_site << "\t" << &*site_it[cur_site] << "\t" << s << *site_it[s];

    move_kinks_iterator_cyclically_backward(previt,s);
    while (!is_time0_cyclically_between_time1_and_time2(instime, previt->time(), site_it[s]->time() ) ) {
      site_it[s] = previt;
      move_kinks_iterator_cyclically_backward(previt,s);
    }
    nbs[i] = site_it[s]->before();
  }

  if (is_worm_rising)
  {
    local_transition_weight[0] = tx_plus(cur_site) * nbs[0];
    local_transition_weight[1] = ty_plus(cur_site) * nbs[1];
    local_transition_weight[2] = tz_plus(cur_site) * nbs[2];
    local_transition_weight[3] = tx_minus(cur_site) * nbs[3];
    local_transition_weight[4] = ty_minus(cur_site) * nbs[4];
    local_transition_weight[5] = tz_minus(cur_site) * nbs[5];
  }
  else
  {
    local_transition_weight[0] = (nbs[0] == _statemax ? 0 : (tx_plus(cur_site) * (nbs[0]+1)));
    local_transition_weight[1] = (nbs[1] == _statemax ? 0 : (ty_plus(cur_site) * (nbs[1]+1)));
    local_transition_weight[2] = (nbs[2] == _statemax ? 0 : (tz_plus(cur_site) * (nbs[2]+1)));
    local_transition_weight[3] = (nbs[3] == _statemax ? 0 : (tx_minus(cur_site) * (nbs[3]+1)));
    local_transition_weight[4] = (nbs[4] == _statemax ? 0 : (ty_minus(cur_site) * (nbs[4]+1)));
    local_transition_weight[5] = (nbs[5] == _statemax ? 0 : (tz_minus(cur_site) * (nbs[5]+1)));
  }

  total_local_transition_weight = std::accumulate(&local_transition_weight[0], &local_transition_weight[zcoord+1], 0.);

  cummulative_local_transition_weight[0] = local_transition_weight[0];
  for (dim_type i=1; i <= zcoord; ++i)
  {
    cummulative_local_transition_weight[i] = local_transition_weight[i] + cummulative_local_transition_weight[i-1];
  }
}



inline bool Lowa::is_kink_encountered_forward(Lowa::time_type& p, Lowa::obs_type const dE, Lowa::time_type& t2)
{
  // see if the next interaction can be reached when moving to the right
  t2 = worm_it->time() - worm_head.time();          // t2 = time between worm head and the next interaction
  if (t2 <= 0.) t2 += beta;
  if (dE == 0.) return (true);
  time_type t1 = p/dE;
  while (std::abs(t1 - t2) < tol) {
    p = -std::log(random_real());
    t1 = p/dE;
  }                                                 // p is automatically chosen again if it is not suitable, ie. t1 = p/dE
  if (t1 > t2) return (true);
  return (false);
}



inline bool Lowa::is_kink_encountered_backward(Lowa::time_type& p, Lowa::obs_type const dE, Lowa::time_type& t2)
{
  // see if the next interaction can be reached when moving to the left
  t2 = worm_head.time() - worm_it->time();
  if (t2 <= 0.) t2 += beta;
  if (dE == 0.) return (true);
  time_type t1 = p/dE;
  while (std::abs(t1 - t2) < tol) {
    p = -std::log(random_real());
    t1 = p/dE;
  }
  if (t1 > t2) return (true);   // take care when t1 and t2 get close...
  return (false);
}



#endif 
