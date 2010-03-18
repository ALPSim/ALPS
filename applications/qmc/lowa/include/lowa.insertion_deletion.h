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

/* $Id: lowa.insertion_deletion.h 3520 2010-03-03 16:55:00Z tamama $ */


/*
 *
 * 1) Code modification      -- partly done
 * 2) Replacing raw pointers -- not done yet
 *
 */


#ifndef ALPS_APPLICATIONS_LOWA_INSERTION_DELETION_H
#define ALPS_APPLICATIONS_LOWA_INSERTION_DELETION_H


bool Lowa::worm_from_insertion_to_deletion()
{
  if (is_worm_diagonal)  // preparing to insert worm-pair
  { 
    ++_Z_dnsmat;

    // proposing wormpair insertion
    ++MCstep;
    ++no_of_proposed_worm_insertions;

    is_worm_moving_forward = ((random_real() < 0.5) ? false : true);
    site_type wormpair_site   = static_cast<site_type>(random_real()*_N);   // choosing site randomly
    time_type wormpair_time   = random_real()*beta;
  
    // find suitable place to make a worm pair
    if (_kinks_on_lattice[wormpair_site].size() != 1) // dummy kink always there
    {
      kinks_iterator_type prevworm_it;
      prevworm_it = ( site_it[wormpair_site] == _kinks_on_lattice[wormpair_site].begin() ? _kinks_on_lattice[wormpair_site].end() : site_it[wormpair_site] );
      --prevworm_it;
      /*
       * Lode's comment:
       *  you cannot create a worm at the same place where the previous one was removed...
       *  you really need to go around in space-time. maybe for insulating phases
       *  the next couple of lines are time consuming
       */
      while (!is_time0_cyclically_between_time1_and_time2(wormpair_time, prevworm_it->time(), site_it[wormpair_site]->time()))
      {
        prevworm_it = site_it[wormpair_site];
        move_kinks_iterator_cyclically_forward(site_it[wormpair_site],wormpair_site);
      }
    }
    worm_it = site_it[wormpair_site];
    fock_basis_type state_enclosed_within_wormpair_at_creation = site_it[wormpair_site]->before();

    bool l_t, mkworm = true;
    // detailed balance slightly differs for hard-core bosons
    // see also density matrix for that
#ifdef HARDCORE
    l_t = true;
    is_worm_rising = (state_enclosed_within_wormpair_at_creation == 0 ? false : true);
#else
    if ((state_enclosed_within_wormpair_at_creation > 0) && (state_enclosed_within_wormpair_at_creation < _statemax)) {
      random_real() < 0.5 ? is_worm_rising = 1 : is_worm_rising = 0;
    }
    else 
    {
      if (state_enclosed_within_wormpair_at_creation == 0) 
      {
        if (random_real() < 0.5)  {  is_worm_rising = false;  }
        else                      {  l_t = true; mkworm = false; }  
      }
      else { // startocc equals _statemax, worm should only be created with prob 1/2
        if (random_real() < 0.5)  {  is_worm_rising = true;  }
        else                      {  l_t = false; mkworm = false; }
      }
    }
  
    if (!mkworm) // density matrix is a bit tricky
    {
#ifdef MEASURE_DENSITY_MATRIX
#ifdef HARDCORE
#ifdef TRAPPEDSYSTEM
      av_dnsmat[0]     += (l_t ? 0. : _statemax);
      av_dnsmat_inf[0] += (l_t ? 0. : _statemax);
#else
      av_dnsmat[0]     += (l_t ? 0. : _statemax);
      av_dnsmat_inf[0] += (l_t ? 0. : _statemax);
#endif
#else
#ifdef TRAPPEDSYSTEM
      av_dnsmat[0]     += (l_t ? 0. : _statemax);
      av_dnsmat_inf[0] += (l_t ? 0. : _statemax); //dns[0];
#else
      av_dnsmat[0]     += (l_t ? 0. : _statemax);
      av_dnsmat_inf[0] += (l_t ? 0. : _statemax); //dns[0];
#endif
#endif
#endif
      //std::cout << "\n No worm inserted, mkworm = 0";
  
      return false;
    }
  
#endif
  
    ++no_of_accepted_worm_insertions;
    new_measurement = true;
  
    // we can now really create the worms
    worm_tail.set_from(wormpair_site);
    worm_tail.set_to(wormpair_site);
    worm_tail.set_time(wormpair_time);
    worm_head = worm_tail;
  
    is_worm_diagonal = false;
  
    // the worm tail is added to the interaction std::list.
    // the worm head is not, but we always need to know
    // where it is and what the density just before/after
    // the worm head is
    if (is_worm_moving_forward) {
      if ( is_worm_rising ) {
        worm_tail.set_before(state_enclosed_within_wormpair_at_creation);
        worm_tail.set_after(--state_enclosed_within_wormpair_at_creation);
        worm_head.set_before(worm_tail.after());
        worm_head.set_after(worm_tail.before());
        _state_enclosed_within_wormpair_at_creation = state_enclosed_within_wormpair_at_creation;
      }
      else {
        worm_tail.set_before(state_enclosed_within_wormpair_at_creation);
        worm_tail.set_after(++state_enclosed_within_wormpair_at_creation);
        worm_head.set_before(worm_tail.after());
        worm_head.set_after(worm_tail.before());
        _state_enclosed_within_wormpair_at_creation = state_enclosed_within_wormpair_at_creation;
      }
      site_it[wormpair_site] = _kinks_on_lattice[wormpair_site].insert(site_it[wormpair_site], worm_tail);
      reset_assoc_upon_insert(wormpair_site, site_it[wormpair_site]);
      move_kinks_iterator_cyclically_forward(site_it[wormpair_site],wormpair_site);
      worm_it = site_it[wormpair_site];
    }
    else {
      if ( is_worm_rising ) {
        worm_tail.set_after(state_enclosed_within_wormpair_at_creation);
        worm_tail.set_before(--state_enclosed_within_wormpair_at_creation);
        worm_head.set_before(worm_tail.after());
        worm_head.set_after(worm_tail.before());
        _state_enclosed_within_wormpair_at_creation = state_enclosed_within_wormpair_at_creation;
      }
      else {
        worm_tail.set_after(state_enclosed_within_wormpair_at_creation);
        worm_tail.set_before(++state_enclosed_within_wormpair_at_creation);
        worm_head.set_before(worm_tail.after());
        worm_head.set_after(worm_tail.before());
        _state_enclosed_within_wormpair_at_creation = state_enclosed_within_wormpair_at_creation;
      }
      site_it[wormpair_site] = _kinks_on_lattice[wormpair_site].insert(site_it[wormpair_site], worm_tail);
      reset_assoc_upon_insert(wormpair_site, site_it[wormpair_site]);
      move_kinks_iterator_cyclically_backward(site_it[wormpair_site],wormpair_site);
      worm_it = site_it[wormpair_site];
    }
#ifdef DEBUGMODE
    std::cout << "\nInitial head : " << worm_head;
    std::cout << "\nInitial tail : " << worm_tail;
    std::cout << "\nInitialized worm..." << worm_tail << std::endl;
    std::cout << "\nis_worm_rising " << is_worm_rising << "\tdir " << is_worm_moving_forward;
    //std::cin.get();

    //print_configuration(std::cout);
#endif
  }

  do {

    time_type pexp = -std::log(random_real());
#ifdef DEBUGMODE
    std::cout << "\nExponential jump time : " << pexp  << " dir " << is_worm_moving_forward<< "\t" << is_worm_rising << std::endl;
#endif
    while (pexp!=0.) {  //exponential time jump
      ++MCstep;
      worm_propagation(pexp);
    }
    if (MCstep >= Nsave)  {  return true;  }

  } while (!is_worm_diagonal);

#ifdef DEBUGMODE
  std::cout << "\n*** DIAGONAL *** ";
  print_configuration(std::cout);
  std::cout << "\n*****************";
  if (!test_conf()) {std::cin.get();};
#endif
  return false;
}



inline bool Lowa::is_time0_cyclically_between_time1_and_time2(Lowa::time_type const t0, Lowa::time_type const t1, Lowa::time_type const t2) const
{ // clockwise-cyclic order, mod beta
  if (t2 > t1) {  if ((t2 >= t0) && (t0 >  t1))  {  return true;  }  }
  else         {  if ((t0 >  t1) || (t0 <= t2))  {  return true;  }  }
  return false;
}



/*
 * Tama's comment : 
 *  1) Making the iterator of a local-site kink_list behave cylindrically in imaginary time
 *  2) tau is beta-periodic in imaginary time for finite temperatue bosonic system : refer to many body texts
 */
  
inline void Lowa::move_kinks_iterator_cyclically_forward(Lowa::kinks_iterator_type& base_iterator, Lowa::site_type const site)
{
  ++base_iterator;
  if (base_iterator == _kinks_on_lattice[site].end())   {  base_iterator = _kinks_on_lattice[site].begin();  }
}

inline void Lowa::move_kinks_iterator_cyclically_backward(Lowa::kinks_iterator_type& base_iterator, Lowa::site_type const site)
{
  if (base_iterator == _kinks_on_lattice[site].begin()) {  base_iterator = _kinks_on_lattice[site].end();    }
  --base_iterator;
}



void Lowa::reset_assoc_upon_insert(Lowa::site_type cursite, Lowa::kinks_iterator_type it)
{
  // go one up on the operator string for the current kink and look to what is pointing
  kinks_iterator_type ito = it;
  move_kinks_iterator_cyclically_forward(ito,cursite);

  for (site_type j = 0; j < zcoord; j++) {
    site_type s = nb(cursite, j);
    kinks_iterator_type itl;

    itl = ito->assoc(j);
    kinks_iterator_type itp = itl;
    move_kinks_iterator_cyclically_backward(itp,s);
    it->set_assoc(j, itl);
    while (!is_time0_cyclically_between_time1_and_time2(it->time(), itp->time(), itl->time())) {
      itl = itp;
      move_kinks_iterator_cyclically_backward(itp,s);
      it->set_assoc(j, itl);
    }
  }

  // interactions on the neighbors below the time of it might also change...
  for (site_type j = 0; j < zcoord; j++) {
    site_type oppdir = (j + dim) % zcoord;
    site_type s = nb(cursite, j);
    kinks_iterator_type itl;
    itl = it->assoc(j);
    kinks_iterator_type itp = itl;
    move_kinks_iterator_cyclically_backward(itp,s);
    kinks_iterator_type itw = itp->assoc(oppdir);
    if (itl->time() == it->time()) itl->set_assoc(oppdir, it);
    while ((itp->time() != itw->time())  && ( is_time0_cyclically_between_time1_and_time2(it->time(), itp->time(), itw->time()) )) {
      itp->set_assoc(oppdir, it);
      move_kinks_iterator_cyclically_backward(itp,s);
      itw = itp->assoc(oppdir);
    }
  }
}



void Lowa::reset_assoc_upon_delete(Lowa::site_type cursite, Lowa::kinks_iterator_type it)
{
  if (_kinks_on_lattice[cursite].size() == 1) {
     for (site_type j = 0; j < zcoord; j++) {
      site_type oppdir = (j + dim) % zcoord;
      site_type s = nb(cursite, j);
          for (kinks_iterator_type itt=_kinks_on_lattice[s].begin(); itt != _kinks_on_lattice[s].end(); ++itt)
            itt->set_assoc(oppdir, _kinks_on_lattice[cursite].begin());
          }
  }
  else {
    kinks_iterator_type newpoint=it;
    move_kinks_iterator_cyclically_forward(newpoint,cursite);
    kinks_iterator_type itp = it;
    move_kinks_iterator_cyclically_backward(itp,cursite);
    // on nb sites, go down until you find an kink that points at the current kink
    for (site_type j = 0; j < zcoord; j++) {
      site_type oppdir = (j + dim) % zcoord;
      site_type s = nb(cursite, j);
      kinks_iterator_type it_end=it->assoc(j);
      kinks_iterator_type it_begin=itp->assoc(j);
      kinks_iterator_type itt=it_begin;
      if ((itt == it_end) && (itt == dummy_kink_it[s])) {
        for (kinks_iterator_type itt=_kinks_on_lattice[s].begin(); itt != _kinks_on_lattice[s].end(); ++itt) {
          if (itt->assoc(oppdir) == it) itt->set_assoc(oppdir, newpoint);
        }
      }
      else {
        while (itt != it_end) {
          if (itt->assoc(oppdir) == it) {
            itt->set_assoc(oppdir, newpoint);
          }
          move_kinks_iterator_cyclically_forward(itt,s);
        }
        if (it_end->assoc(oppdir) == it) {
          it_end->set_assoc(oppdir, newpoint);
        }
      } //else
    }
  }
}






#endif
