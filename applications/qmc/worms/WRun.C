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

//- Worm generation -------------------------------------------------------

int64_t WRun::make_worm()
{
  int64_t length=0;
  if(create_worm())  {
    stat[0]+=1.;
    do {
      length++;
      if(measure_green_function_)
        measure_green();

      // can/should we annihilate the two kinks?
      if(worm_head[0].site()==worm_head[1].site() && random_real()<P_REMOVE) {
        stat[1]+=1.;
        stat[3]+=1.;
        if(annihilate_worm())
          break;
        else
          continue;
      }

      // choose updating method randomly
      int update_method=random_int(3);

      /*
      int head_num = random_int(2);
      wormhead_type& worm_end  = worm_head[head_num];
      wormhead_type& worm_tail = worm_head[1-head_num];
      
      if(head_num != current_head_num) 
        subinterval_valid = false;
      current_head_num = head_num;
      */

      wormhead_type& worm_end  = worm_head[0];
      wormhead_type& worm_tail = worm_head[1];
      subinterval_valid = false;

      stat[2]+=1;
      stat[3]+=1;
      switch(update_method) {
      case 0: 
        // update: insert jump
        if(num_neighbors(worm_end.site())) {
          int nb  = random_int(num_neighbors(worm_end.site()));
          int dir = 2*random_int(2)-1;
          // std::cerr << "calling insert jump" << std::endl;
          insert_jump(worm_end,worm_tail,dir,nb);
        }
        break;
      case 1:
        // update: remove jump
        if(num_neighbors(worm_end.site())) {
          int nb  = random_int(num_neighbors(worm_end.site()));
          int dir = 2*random_int(2)-1;
          // std::cerr << "calling remove jump" << std::endl;
          remove_jump(worm_end,worm_tail,dir,nb);
        }
        break;
      case 2:
        // update: shift kink
        // std::cerr << "calling shift kink" << std::endl;
        shift_kink(worm_end,worm_tail);
        break;
      default:
        boost::throw_exception(std::logic_error("invalid update_method"));
        break;
      }
#ifdef CHECK_OFTEN
          check_spins();
#endif
    } while(true); // break above by annihilating kinks
  }
  return length;
}   // WRun::make_worm


inline double WRun::Delta_H0(const site_descriptor& s, const state_type& state1, const state_type& state2, 
                             std::vector<state_type>& adjacent_state) {
  //
  //  determines the energy difference of state1 and state2 on site s taking into account 
  //  the states of all neighbors as given in adjacent_state.
  //  This function doublecounts the bond energy terms.
  //

  // determine onsite energy
  double Result = onsite_energy(state1, s) - onsite_energy(state2, s);


  if (nonlocal) {
    // determine neighbor energy
    std::vector<state_type>::iterator adj_state_it = adjacent_state.begin();
    neighbor_bond_iterator nbi, nbi_end;
    for(boost::tie(nbi, nbi_end) = neighbor_bonds(s); nbi != nbi_end; ++nbi) {
      Result += neighbor_energy(state1, *adj_state_it, *nbi);
      Result -= neighbor_energy(state2, *adj_state_it, *nbi);
      adj_state_it++;
    }  
  }
  return Result;
}   // WRun::Delta_H0

inline double WRun::Delta_H0(const site_descriptor& s1, const state_type& initial_state1, const state_type& final_state1, std::vector<state_type>& adj_state1,
                      const site_descriptor& s2, const state_type& initial_state2, const state_type& final_state2, std::vector<state_type>& adj_state2,
                      const bond_descriptor& this_bond) {
  //
  //  determines the energy of state1 on site s1 and state2 on site s2 at time t1 taking into account 
  //  the states of all neighbors at that particular time, thereby avoiding double counting of the
  //  bond between s1 and s2.
  //

  // determine onsite energy
  double initial_energy = onsite_energy(initial_state1, s1) + onsite_energy(initial_state2, s2);
  double final_energy   = onsite_energy(final_state1, s1) + onsite_energy(final_state2, s2);

  if (nonlocal) {
    // determine neighbor energy of states on site s1, avoid double counting
    neighbor_bond_iterator nbi, nbi_end;

    std::vector<state_type>::iterator adj_state_it = adj_state1.begin();
    for(boost::tie(nbi, nbi_end) = neighbor_bonds(s1); nbi != nbi_end; ++nbi) {
      if( target(*nbi) != s2 ) {
        initial_energy += neighbor_energy(initial_state1, *adj_state_it, *nbi);
        final_energy   += neighbor_energy(final_state1, *adj_state_it, *nbi);
      }
      adj_state_it++;
    }  

    // determine neighbor energy of states on  site s2, avoid double counting
    adj_state_it = adj_state2.begin();
    for(boost::tie(nbi, nbi_end) = neighbor_bonds(s2); nbi != nbi_end; ++nbi) {
      if( target(*nbi) != s1 ) {
        initial_energy += neighbor_energy(initial_state2, *adj_state_it, *nbi);     
        final_energy   += neighbor_energy(final_state2, *adj_state_it, *nbi);     
      }
      adj_state_it++;
    }  

    // add edge between the two sites
    initial_energy += neighbor_energy(initial_state1, initial_state2, this_bond);
    final_energy   += neighbor_energy(final_state1, final_state2, this_bond);
  }

  return initial_energy - final_energy;
}   // WRun::Delta_H0

//- Subinterval algorithms ------------------------------------------------

void WRun::traverse_subintervals(const time_struct& t1, const time_struct& t2, const double& time_interval, 
                                 const site_descriptor& s,
                                 const state_type& state1, const state_type& state2, double cutoff=0.0) {
  // static variables
  static std::vector<cyclic_iterator> adjacent_kink(8);
  static std::vector<state_type> adjacent_state(8);
  static std::vector<int> closest_neighbor(2);
  static subinterval_info this_interval;

  std::vector<cyclic_iterator>::iterator adj_kink_it;

  // integrated weight, time and modification factor for stat. weight
  double iWeight = 0.0;
  double iTime   = 0.0;
  double factor  = 1.0;

  // initialize subinterval information / static variables
  subinterval.clear();
  adjacent_kink.clear(); 
  adjacent_state.clear(); 
  closest_neighbor.clear();

  // subinterval times
  time_struct start_time = t1;
  time_struct end_time   = t2;

  // set up valid neighbors
  for(int nb=0; nb<num_neighbors(s); nb++) {
    cyclic_iterator nn = first_kink(neighbor(s, nb)); 
    if( nn.valid() ) { 
      adjacent_kink.push_back( adjacent(nn, t1) );
      adjacent_state.push_back( adjacent_kink.back()->state() );
      adjacent_kink.back()++;
    }
    else {
      adjacent_kink.push_back( nn );
      adjacent_state.push_back( initial_state(neighbor(s, nb)) );
    }
  }

  double delta_e = Delta_H0(s, state1, state2, adjacent_state);

  // iterate all subintervals
  do {
    iTime = time_interval;
    end_time = t2;
    closest_neighbor.clear();

    // find closest kinks
    int nb = 0;
    for(adj_kink_it = adjacent_kink.begin(); adj_kink_it != adjacent_kink.end(); adj_kink_it++) {
      if( adj_kink_it->valid() ) {
        if((*adj_kink_it)->time() == end_time) 
          closest_neighbor.push_back(nb);
        else 
          if((*adj_kink_it)->time()-t1 < iTime) {
            end_time = (*adj_kink_it)->time();
            iTime = (*adj_kink_it)->time()-t1;
            
            closest_neighbor.clear();
            closest_neighbor.push_back(nb);
          }
      }
      nb++;
    }   
    
    // insert subinterval
    if( iTime<=time_interval ) {
      double delta_t = iTime - (subinterval.size() ?  subinterval.back().integrated_time : 0.);            

      iWeight += factor * integrated_weight(beta*delta_e, delta_t);
   
      // add information about current subinterval
      this_interval.start_time        = start_time;
      this_interval.end_time          = end_time;
      this_interval.delta_t              =        delta_t;
      this_interval.delta_e           = delta_e;
      this_interval.integrated_time   = iTime;
      this_interval.integrated_weight = iWeight;
      subinterval.push_back(this_interval);
    }

    if( iTime<time_interval && (!cutoff || iWeight<cutoff) ) {
      // update subinterval times
      start_time = end_time;

      // update prefactor for integrated weight
      factor *= exp(-beta*this_interval.delta_e*this_interval.delta_t);

      // update adjacent kinks/states of closest neighbors
      for(int i=0; i<closest_neighbor.size(); ++i) {
        int cn = closest_neighbor[i];

        // update energy difference
        Update_Delta_H0(delta_e, state1, state2, adjacent_state[cn], adjacent_kink[cn]->state(),
                        *(neighbor_bonds(s).first+cn));

        // update state
        adjacent_state[cn] = adjacent_kink[cn]->state();

        // update kink/validity
        if( (adjacent_kink[cn]+1)->time()-t1 > (adjacent_kink[cn]->time()-t1) )
          adjacent_kink[cn]++; 
        else 
          adjacent_kink[cn].invalidate();
      }
    }
  } while( iTime<time_interval && (!cutoff || iWeight<cutoff) );

  if(!subinterval.size())
    boost::throw_exception(std::logic_error("No subinterval found in traverse_subintervals"));

  subinterval_valid = true;
}   // WRun::traverse_subintervals

void WRun::traverse_subintervals(wormhead_type& head) {
  cyclic_iterator h = head.kink();
  cyclic_iterator start = h-1;
  cyclic_iterator end   = head.kink()+1;

  state_type state1 = start->state();
  state_type state2 = h->state();

  assert( state1!=state2 );

  traverse_subintervals(start->time(), end->time(), end->time()-start->time(), head.site(), state1, state2);
}   // WRun::traverse_subintervals

//- Worm creation probability ---------------------------------------------

double WRun::worm_creation_probability(const time_struct& t1, const double& time_interval, const site_descriptor& s,
                                       const state_type& state1, const state_type& state2, double cutoff=0.0) {

  traverse_subintervals(t1, t1+time_interval, time_interval, s, state1, state2, cutoff);
  subinterval_valid = false;
  return subinterval.back().integrated_weight;
}   // WRun::worm_creation_probability

//------------------------------------------------------------------------

std::pair<time_struct, time_struct> WRun::adjacent_subinterval(const site_descriptor& s, const time_struct& t1,
                                                                          std::vector<state_type>& adjacent_state) {
  /*
    determines the next time step where there is a kink on site s or any neighboring site.
    Function used by insert/remove jump.
  */
  time_struct left_time  = t1;
  time_struct right_time = t1;

  adjacent_state.clear();

  // check this site
  cyclic_iterator nn = first_kink(s);
  if( nn.valid() ) {
    cyclic_iterator adjacent_kink = adjacent(nn, t1);
    left_time  = adjacent_kink->time()==t1 ? (adjacent_kink-1)->time() : adjacent_kink->time();
    right_time = (adjacent_kink+1)->time();
  }


  // check valid neighbors
  for(int nb=0; nb<num_neighbors(s); nb++) {
    cyclic_iterator nn = first_kink(neighbor(s, nb));
    if( nn.valid() ) {
      //cyclic_iterator first=nn;
      cyclic_iterator adjacent_kink = adjacent(nn, t1);
      time_struct left_t2  = adjacent_kink->time()==t1 ? (adjacent_kink-1)->time() : adjacent_kink->time();
      time_struct right_t2 = (adjacent_kink+1)->time();

      if( right_t2-t1<right_time-t1 ) 
        right_time = right_t2;
      if( t1-left_t2<t1-left_time )
        left_time = left_t2;
 
      adjacent_state.push_back( adjacent_kink->state() );
    }
    else 
    adjacent_state.push_back( initial_state(neighbor(s, nb)) );
  }

  return std::make_pair(left_time, right_time);
}   // WRun::adjacent_subinterval_time

//- Update routines -------------------------------------------------------

void WRun::shift_kink(wormhead_type& head,wormhead_type&)
{
  cyclic_iterator h = head.kink();
  cyclic_iterator start = h-1;
  cyclic_iterator end   = h+1;

  time_struct newtime;
 

  do { 
  if (!nonlocal) {
    double time_interval  = end->time() - start->time();
    double delta_e = onsite_energy(start->state(),head.site()) -onsite_energy(h->state(),head.site());
    newtime = start->time() + 
      finite_exponential_random<engine_type>(*engine_ptr,-beta*delta_e,time_interval)();
  }
  else {
    double random_x = random();

    if(!subinterval_valid) 
      traverse_subintervals(head);
    double i_weight = subinterval.back().integrated_weight;

    int isub = 0;
    while(random_x > subinterval[isub].integrated_weight/i_weight) 
      ++isub;

    double subinterval_weight = subinterval[isub].integrated_weight-(isub ? subinterval[isub-1].integrated_weight : 0.);
    double random_y = (random_x*i_weight-(isub ? subinterval[isub-1].integrated_weight :0.))/subinterval_weight;

    newtime = subinterval[isub].start_time 
      + new_finite_exponential_random(random_y, -beta*subinterval[isub].delta_e, subinterval[isub].delta_t)();
  }
// avoid problems with hitting the end or start of the interval
  } while (newtime == end->time() || newtime == start->time());

  head = move_kink(head.site(),h,newtime);
#ifdef USE_VECTOR
  other_head.invalidate(head.site());
#endif
}   // WRun::shift_kink

//-------------------------------------------------------------------------

void WRun::insert_jump(wormhead_type& head, wormhead_type&,int dir,int nb)
{
  // creation: up                down
  // s2   s3     s2  s3        s2   s3     s2 s3
  // s2   s3  -> s1  s4        s1   s3  -> s2 s4
  // s1   s3     s1  s3        s1   s3     s1 s3

  double time;
  static std::vector<state_type> adjacent_state1;
  static std::vector<state_type> adjacent_state2;

  bond_descriptor bond = *(neighbor_bonds(head.site()).first+nb);
  site_descriptor newsite = target(bond);
  cyclic_iterator h = head.kink();
  cyclic_iterator start2 = h->adjacent(first_kink(newsite));
  cyclic_iterator end2=(start2.valid() ? start2+1 : start2);
  time_struct th = h->time();

  state_type state1 = (h-1)->state();
  state_type state2 = h->state();
  state_type state3 = start2.valid() ? start2->state() : initial_state_[newsite];
  state_type state4;

  // determine matrix element
  bool do_create = (dir>0 ? (state1<state2) : (state2<state1));
  double q=hopping_matrix_element(/*(dir>0 ? state2 : state1),state3,*/bond) *
    creation_matrix_element(state3,do_create,newsite);
  state4 = (do_create ? create(state3) : annihilate(state3));
  if(q==0.)
    return;

  if (!nonlocal) {
    // find interval
    if(dir>0) {
      //cyclic_iterator start = h;
      cyclic_iterator end   = h+1;
      time=end->time()-h->time();
      if(end2.valid()) {
        double time2=end2->time()-h->time();
        if (time2<time)
          time=time2;
      }
    }
    else {
      cyclic_iterator start = h-1;
      ///cyclic_iterator end   = h;
      time=h->time()-start->time();
      if(start2.valid()) {
        double time2=h->time()-start2->time();
        if (time2<time)
          time=time2;
      }
    }
  }
  else {

    // find interval
    time_struct left_t1, right_t1; 
    boost::tie(left_t1, right_t1) = adjacent_subinterval(head.site(), head.time(), adjacent_state1);

    time_struct left_t2, right_t2;
    boost::tie(left_t2, right_t2) = adjacent_subinterval(newsite, head.time(), adjacent_state2);

    time_struct t1 = dir>0 ? right_t1 : left_t1;
    time_struct t2 = dir>0 ? right_t2 : left_t2;

    if(dir>0)
      time = (t1-th < t2-th) ? t1-th : t2-th;
    else
      time = (th-t1 < th-t2) ?th-t1 : th-t2;
  }

  assert(time>0 && time<=1.);

  if(time<=0 || time>1.) {
    std::cerr << "time = " << time << std::endl;
    check_spins();
    boost::throw_exception(std::logic_error("invalid time in insert_jump"));
  }
  
  double statistical_weight; 
  if (!nonlocal)
     statistical_weight=(onsite_energy(state1,head.site()) - onsite_energy(state2,head.site()))*dir +
                        onsite_energy(state4,newsite) - onsite_energy(state3,newsite);
  else {
    if(dir>0)
      statistical_weight 
        = Delta_H0(head.site(), state1, state2, adjacent_state1, 
                 newsite, state4, state3, adjacent_state2, bond);
    else
      statistical_weight 
        = Delta_H0(head.site(), state2, state1, adjacent_state1,
                 newsite, state4, state3, adjacent_state2, bond);
  }
  
  // determine acceptance rate
  double pacc = fabs(q)*integrated_weight(statistical_weight, beta*time);
  pacc *= double(num_neighbors(head.site())) / double(num_neighbors(newsite));
  
  // accept or reject the insertion
  if (pacc < 1.&& pacc < random_real())
    return;

  // determine jump time
  time_struct jump_time=th +dir*
     finite_exponential_random<engine_type>(*engine_ptr,-beta*statistical_weight,time)();

  // insert jump
  time_struct t=th;
  iterator it=move_kink(head.site(),h,jump_time);
  it->set_id(++last_id_);
#ifdef USE_VECTOR
  other_head.invalidate(head.site());
#endif
  head.set_site(newsite);
#ifdef USE_VECTOR
  other_head.invalidate(head.site());
#endif
  it=insert_kink(newsite,end2,kink_type(dir > 0 ? jump_time : t,state3));
  if (dir> 0)
    it->set_id(last_id_);
#ifndef USE_VECTOR
  if (dir<0) 
    head=it;
#endif
  it=insert_kink(newsite,it,kink_type(dir > 0 ? t :  jump_time,state4));
  if (dir>0) 
    head=it;
  else
    it->set_id(last_id_);

  num_kinks++;
  if(is_signed_) 
    Sign *= q<0 ? 1 : -1;

  subinterval_valid = false;
}   // WRun::insert_jump

//-------------------------------------------------------------------------

void WRun::remove_jump(wormhead_type& head, wormhead_type& , int dir, int nb)
{
  bond_descriptor bond = *(neighbor_bonds(head.site()).first+nb);
  site_descriptor newsite = target(bond);

  cyclic_iterator h = head.kink();
  cyclic_iterator k1=(dir>0 ? h+1 : h-1);
  cyclic_iterator k2=h->adjacent(first_kink(newsite));
  if (!k2.valid())
    return;
  if (dir>0)
    k2++;
  if(k1->time()!=k2->time()) 
    return;
  if (k1->id()!=k2->id()) // this is not the same bond
    return;

  state_type state1=(k1-1)->state();
  state_type state2=k1->state();
  state_type state3=(k2-1)->state();
  state_type state4=k2->state(); 
  
  if ((dir>0 && state2 != (h-1)->state()) || (dir<0 && state1 != h->state()))
    return; // wrong type of jump
  
  // meaning of 3 and 4 exchanged for dir<0 compared to adding of jump
  double q = hopping_matrix_element(/*state1,state3,*/bond)*
    creation_matrix_element(state1,state1<state2,head.site());
  if(fabs(q)<1e-10) {      
    std::cerr << " q = " << q << " state1 = " << int(state1) << " state2 = " << int(state2) 
              << " bond_type = " << bond_type[bond] << " id " << k1->id() << " " << k2->id() << std::endl;
    print_spins();
    boost::throw_exception(std::logic_error("zero matrix element in remove_jump"));
  }
  
  double time;
  static std::vector<state_type> adjacent_state1;
  static std::vector<state_type> adjacent_state2;

  if (!nonlocal) {
    double time1 = k1->time()-(k1-1)->time();
    double time2 = k2->time()-(k2-1)->time();
    double left_time = (time1<time2 ? time1 : time2);
  
    time1 = (k1+1)->time() - k1->time();
    time2 = (k2+1)->time() - k2->time();
    double right_time = (time1<time2 ? time1 : time2);
  
    time = left_time + right_time;
  } 
  else {
    time_struct jump_time = k1->time();

  
    time_struct left_t1, right_t1;
    time_struct left_t2, right_t2;
  
    boost::tie(left_t1, right_t1) = adjacent_subinterval(head.site(), jump_time, adjacent_state1);
  
    // still allowed to remove jump?
    if( dir>0 ? jump_time-head.time()>jump_time-left_t1 : head.time()-jump_time>right_t1-jump_time ) 
      return;
  
    boost::tie(left_t2, right_t2) = adjacent_subinterval(newsite, jump_time, adjacent_state2);
  
    // still allowed to remove jump?
    if( dir>0 ? jump_time-head.time()>jump_time-left_t2 : head.time()-jump_time>right_t2-jump_time )
      return;
  
    double left_time  = (jump_time-left_t1<jump_time-left_t2) ? jump_time-left_t1 : jump_time-left_t2;
    double right_time = (right_t1-jump_time<right_t2-jump_time) ? right_t1-jump_time : right_t2-jump_time;
    time = left_time + right_time;
  }
  
  double statistical_weight;
  if (!nonlocal)
   statistical_weight = -dir*(
                -onsite_energy(state1,head.site()) + onsite_energy(state2,head.site())
                -onsite_energy(state3,newsite) + onsite_energy(state4,newsite));
  else
    statistical_weight 
    = dir*Delta_H0(head.site(), state1, state2, adjacent_state1,
                   newsite, state3, state4, adjacent_state2, bond);
  
  // determine acceptance rate
  double pacc = 1./(fabs(q)*integrated_weight(statistical_weight, time*beta));
  pacc *= double(num_neighbors(head.site())) / double(num_neighbors(newsite));

  // accept or reject the insertion
  if (pacc< 1. && pacc <= random_real())
    return;
  
  // remove jump
  if (k1->time()>h->time()) {
    erase_kink(head.site(),k1);
    erase_kink(head.site(),h);
  }
  else {
    erase_kink(head.site(),h);
    erase_kink(head.site(),k1);
  }
  if (kinks[head.site()].empty())
    initial_state_[head.site()]=(dir > 0 ? state2 : state1);
#ifdef USE_VECTOR
  other_head.invalidate(head.site());
#endif
  head.set_site(neighbor(head.site(),nb));
#ifdef USE_VECTOR
  other_head.invalidate(head.site());
#endif
  head=move_kink(head.site(),k2,head.time());

  num_kinks--;
  if(is_signed_)
    Sign *= q<0 ? 1 : -1;

  subinterval_valid = false;
}   // WRun::remove_jump

//- Worm creation/annihilation --------------------------------------------

bool WRun::create_worm() 
{
  // choose random site
  site_descriptor s = site(random_int(num_sites()));

  // choose random time t1
  time_struct t1 = random_real();
  assert(t1>=0);

  // determine kink segment
  cyclic_iterator start = adjacent(first_kink(s), t1);
  if (start.valid() && start->time() == t1)
    return false;
  cyclic_iterator end   = (start.valid() ? start+1 : start);

  // determine states
  state_type state1 = (start.valid() ? start->state() : initial_state_[s]);
  bool do_create=(random_real()<0.5);
 if((do_create && state1==max_state()) || (!do_create && state1==min_state()))
    return false;
    
  state_type state2 = (do_create ? create(state1) : annihilate(state1));

  int num=kinks[s].size();

  double time = num ? end->time()-t1 : 1.;

  // calculate acceptance probability
  
  double pacc;
  
  if (nonlocal) {
    assert(!subinterval_valid);
    
    double cme  = creation_matrix_element(state1,do_create,s)/eta * beta;
    if(!cme) 
      return false;
    double wcp  = worm_creation_probability(t1, time, s, state2, state1, 1./cme);
    pacc = cme * wcp;
    subinterval_valid = false;
  }
  else {
    pacc = creation_matrix_element(state1,do_create,s)/eta;
	if (pacc==0.)
		  return false;
	double statistical_weight =
        onsite_energy(state2,s)-onsite_energy(state1,s);
    pacc *= integrated_weight(statistical_weight, beta*time);
  }

  // accept/reject 
  if(pacc<1. && pacc < random_real()) 
    return false;

  // choose time t2 - will be shifted afterwards
  time_struct t2 = t1+0.5*time;
  if (end.valid() && end->time()==t2)
    return false;

  // insert wormheads
  worm_head[0].set_site(s);
  worm_head[1].set_site(s);
  iterator k=insert_kink(s,end,kink_type(t2,state1));
#ifndef USE_VECTOR
  worm_head[1]=k;
#else
  worm_head[1].set_time(t2);
#endif
  worm_head[0]=insert_kink(s,k,kink_type(t1,state2));

  // shift time t2
  shift_kink(worm_head[1],worm_head[0]);
  subinterval_valid = false;

  return true;
}   // WRun::create_worm
 
//-------------------------------------------------------------------------

bool WRun::annihilate_worm()
{
  // determine site
  site_descriptor s = worm_head[0].site();

  int first=random_int(2);
  cyclic_iterator w0=worm_head[first].kink();
  cyclic_iterator w1=worm_head[1-first].kink();
  if(w0+1!=w1) {
    return false;
  }
  
  // determine states
  state_type state1=(w0-1)->state();
  state_type state2=w0->state();

  int num=kinks[s].size()-2;

  // determine kink segment
  cyclic_iterator end   = w1+1;

  time_struct t1 = w0->time();
  double time = num ? end->time()-t1 : 1.;

  double pacc;
  if (nonlocal) {
    cyclic_iterator start = w0-1;
    double wcp =  num > 0 ? worm_creation_probability(t1, time, s, state2, state1) :
      worm_creation_probability(start->time(), time, s, state2, state1);
    pacc=eta/(creation_matrix_element(state1,state2>state1,s) * beta * wcp);
  }
  else {
    double statistical_weight = onsite_energy(state2,s)- onsite_energy(state1,s);
    double wcp = integrated_weight(statistical_weight, beta*time);
    pacc=eta/(creation_matrix_element(state1,state2>state1,s) * wcp);
  }

  if(pacc<1 && pacc<random_real()) 
    return false; // rejected

  if( w0->time()>w1->time() ) {
    erase_kink(s,w0);
    erase_kink(s,w1);
  }
  else {
    erase_kink(s,w1);
    erase_kink(s,w0);
  }

  if(num==0) 
    initial_state_[s]=state1;
  
  have_worm = false;
  subinterval_valid = false;
  return true;
}   // WRun::annihilate_worm
