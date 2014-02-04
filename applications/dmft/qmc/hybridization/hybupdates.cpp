/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <egull@umich.edu>,
 *
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
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
#include <iomanip>
#include"hyb.hpp"
#include"hyblocal.hpp"

//this is the heart of the Monte Carlo procedure: we have the following updates:
//1: change the zero order state, swap an empty orbital versus a filled one
//2: shift an existing segment start- or end-point
//3: insert or remove a new segment
//4: insert or remove an anti-segment
//5: perform a segment flip between different orbtials
//see our review for details of these updates
void hybridization::update(){

  //one sweep is composed of N_MEAS Monte Carlo updates and one measurement (the latter only if thermalized)
  sweeps++;

  for(std::size_t i=0;i<N_meas;++i){
    double update_type=random();
    if(spin_flip){
      if(update_type < 0.1){
        change_zero_order_state_update();
      }else if(update_type < 0){
        shift_segment_update();
      }else if(update_type < 0.5){
        insert_remove_segment_update();
      }else if(update_type < 0.8){
        insert_remove_antisegment_update();
      }else{
        insert_remove_spin_flip_update();
      }
    }else{
      if(update_type < 0.1){
        change_zero_order_state_update();
      }else if(update_type < 0){
        shift_segment_update();
      }else if(update_type < 0.6){
        insert_remove_segment_update();
      }else{
        insert_remove_antisegment_update();
      }
    }

    if(is_thermalized()){
      measure_order();
      if(MEASURE_time){
        local_config.get_F_prefactor(F_prefactor);//compute segment overlaps in local config
        measure_G(F_prefactor);
      }
    }

  }//N_meas

  if(VERBOSE && sweeps%output_period==0 && crank==0) {
    //  if(VERBOSE && crank==0 && boost::chrono::steady_clock::now() - lasttime > delay) {
    //    lasttime = boost::chrono::steady_clock::now();
    int tot_acc=0,cur_prec = std::cout.precision();
    for (int i=0;i<nacc.size();i++) tot_acc += nacc[i];
    std::cout << std::endl << "|------ Simulation details (master only) after " << sweeps << " sweeps ------|" << std::endl;
    std::cout << "  Total acceptance rate = " << std::setprecision(2) << std::fixed;
    std::cout << (((double)tot_acc)/(sweeps*N_meas))*100 << "%" << std::endl;
    std::cout << "  Individual acceptance rates for update " << std::endl;
    for (int i=0;i<nacc.size();i++) {
      std::cout << "     " << update_type[i] << " = ";
      std::cout << std::setprecision(2) << std::fixed << (((double)nacc[i])/(sweeps*N_meas))*100 << "%";
      std::cout << " (proposal rate = ";
      std::cout << std::setprecision(2) << std::fixed << (((double)nprop[i])/(sweeps*N_meas))*100 << "%)" << std::endl;
    }
    std::cout << "|-----------------------------------------------------------------|" << std::endl;
    std::cout.unsetf(std::ios_base::fixed);
    std::cout.precision(cur_prec);
  }
}

void hybridization::change_zero_order_state_update(){
  //choose the orbital in which we do the update
  nprop[0]++;
  int orbital=(int)(random()*n_orbitals);
  
  //changing the zero order state only makes sense if we are at zero order.
  if(!local_config.order(orbital)==0) return;
  
  //propose to change orbital from occuppied to unoccuppied.
  if(local_config.zero_order_orbital_occupied(orbital)){
    double local_weight_change=1./local_config.local_weight_change(segment(0,beta), orbital, false);
    if(std::abs(local_weight_change)>random()){
      nacc[0]++;
      local_config.set_zero_order_orbital_occupied(orbital, false);
      if(local_weight_change<0)
        sign*=-1.;
    }
  }
  //propose to change from unoccuppied to occuppied.
  else{
    double local_weight_change=local_config.local_weight_change(segment(0,beta), orbital, false);
    //std::cout<<cmagenta<<"local weight change is: "<<local_weight_change<<cblack<<std::endl;
    if(std::abs(local_weight_change)>random()){
      nacc[0]++;
      local_config.set_zero_order_orbital_occupied(orbital, true);
      if(local_weight_change<0)
        sign*=-1.;
    }
  }
}
void hybridization::shift_segment_update(){
  ///TODO: implement this update!
}
void hybridization::insert_remove_segment_update(){
  //choose the orbital in which we do the update
  int orbital=(int)(random()*n_orbitals);
  if(random()<0.5){ insert_segment_update(orbital);}
  else            { remove_segment_update(orbital);}
}
void hybridization::insert_remove_antisegment_update(){
  //choose the orbital in which we do the update
  int orbital=(int)(random()*n_orbitals);
  if(random()<0.5){ insert_antisegment_update(orbital);}
  else            { remove_antisegment_update(orbital);}
}
void hybridization::insert_remove_spin_flip_update(){
  //choose the orbital in which we do the update
  int orbital=(int)(random()*n_orbitals);
  spin_flip_update(orbital);
}

void hybridization::insert_segment_update(int orbital){
  nprop[1]++;
  //std::cout<<clred<<"starting insertion update."<<cblack<<std::endl;
  if(local_config.order(orbital)==0 && local_config.zero_order_orbital_occupied(orbital)) return; //can't insert segment, orbital is fully occuppied.
  double t_start=random()*beta; //start time of a segment
  if(local_config.exists(t_start)){ /*std::cerr<<"rare event, duplicate: "<<t_start<<std::endl; */ return;} //time already exists.
  double t_next_segment_start=local_config.find_next_segment_start_distance(t_start,orbital);
  double t_next_segment_end=local_config.find_next_segment_end_distance(t_start,orbital);
  //std::cout<<"============================================"<<cblack<<std::endl;
  //std::cout<<clblue<<"orbital: "<<orbital<<" time is: "<<t_start<<" segment start distance: "<<t_next_segment_start<<" end distance: "<<t_next_segment_end<<cblack<<std::endl;
  //std::cout<<*this<<std::endl;
  //std::cout<<"============================================"<<cblack<<std::endl;
  if(t_next_segment_end < t_next_segment_start) return; //we're trying to create a segment on top of another segment. abort.
  
  //draw an end time
  double t_len = random()*t_next_segment_start;
  double t_end=t_start+t_len;
  if(t_end >= beta) t_end-=beta;
  if(local_config.exists(t_end)){ /*std::cerr<<"rare event, duplicate: "<<t_end<<std::endl; */return;} //time already exists.
  if(t_end<=t_start || t_end<=0.0){ /*std::cerr<<"rare event, zero length segment: "<<t_start<<" "<<t_end<<std::endl; */return;} //time already exists.
  
  //compute local weight of the new segment with t_start and t_end
  segment new_segment(t_start, t_end);
  double local_weight_change=local_config.local_weight_change(new_segment, orbital, false);
  
  //compute hybridization weight change
  double hybridization_weight_change=hyb_config.hyb_weight_change_insert(new_segment, orbital);
  
  //compute the proposal probability ratio
  double permutation_factor=beta*t_next_segment_start/(local_config.order(orbital)+1);
  
  //perform metropolis
  double weight_change=local_weight_change*hybridization_weight_change*permutation_factor;
  
  /*std::cout<<" new segment: "<<new_segment<<std::endl;
   std::cout<<clred<<"weight change: "<<weight_change<<" l: "<<local_weight_change<<" h: "<<hybridization_weight_change<<" p: "<<permutation_factor<<cblack<<std::endl;*/
  
  if(std::abs(weight_change)>random()){
    nacc[1]++;
    if(weight_change < 0) sign*=-1.;
    local_config.insert_segment(new_segment, orbital);
    hyb_config.insert_segment(new_segment, orbital);
  }
}
void hybridization::remove_segment_update(int orbital){
  nprop[2]++;
  //std::cout<<clblue<<"starting removal update."<<cblack<<std::endl;
  int k=local_config.order(orbital);
  
  if(k==0) return; //no point, this is an empty orbital
  
  int segment_nr=(int)(random()*k);
  
  segment segment_to_remove=local_config.get_segment(segment_nr, orbital);
  
  double local_weight_change=1./local_config.local_weight_change(segment_to_remove, orbital, false);
  
  //compute hybridization weight change
  double hybridization_weight_change=1.0/hyb_config.hyb_weight_change_remove(segment_to_remove, orbital);
  
  //compute the proposal probability ratio
  double t_next_segment_start=local_config.find_next_segment_start_distance(segment_to_remove.t_start_,orbital);
  double permutation_factor=local_config.order(orbital)/(beta*t_next_segment_start);
  
  //perform metropolis
  double weight_change=local_weight_change*hybridization_weight_change*permutation_factor;
  
  /*if(segment_nr==k-1){
   std::cout<<" segment_to_remove: "<<segment_to_remove<<std::endl;
   std::cout<<"t_next_segment_start: "<<t_next_segment_start<<std::endl;
   std::cout<<clblue<<"weight change: "<<weight_change<<" l: "<<local_weight_change<<" h: "<<hybridization_weight_change<<" p: "<<permutation_factor<<cblack<<std::endl;
   }*/
  if(std::abs(weight_change)>random()){
    nacc[2]++;
    if(weight_change < 0) sign*=-1.;
//      double fwo = full_weight();
    local_config.remove_segment(segment_to_remove, orbital);
    hyb_config.remove_segment(segment_to_remove, orbital);
//      double fwa = full_weight();
//      std::cout << clgreen<<"weight change removal: "<<fwa<<" control: "<<fwo*std::abs(weight_change)<<std::endl;
  }
}
void hybridization::insert_antisegment_update(int orbital){
  nprop[3]++;
  if(local_config.order(orbital)==0 && !local_config.zero_order_orbital_occupied(orbital)) return; //can't insert an antisegment, orbital is empty.
  double t_start=random()*beta; //start time of the anti segment
  if(local_config.exists(t_start)){ /*std::cerr<<"rare event, duplicate: "<<t_start<<std::endl; */return;} //time already exists.
  double t_next_segment_start=local_config.find_next_segment_start_distance(t_start,orbital);
  double t_next_segment_end=local_config.find_next_segment_end_distance(t_start,orbital);
  
  if(t_next_segment_start < t_next_segment_end) return; //we're trying to create an antisegment where there is no segment abort.
  
  //draw an end time
  double t_len = random()*t_next_segment_end;
  double t_end=t_start+t_len; //((t_len<0.1*beta)?t_len:0.1*beta); //random()*t_next_segment_end;
  if(t_end >= beta) t_end-=beta;
  if(local_config.exists(t_end)){ /*std::cerr<<"rare event, duplicate: "<<t_end<<std::endl; */return;} //time already exists.
  if(t_end<=t_start || t_end<=0.0){ /*std::cerr<<"rare event, zero length segment: "<<t_start<<" "<<t_end<<std::endl; */return;} //time already exists.
  
  //std::cout<<clgreen<<"antisegment insertion update: "<<std::endl<<cblack<<*this<<std::endl;
  //std::cout<<clgreen<<" antisegment start time: (cdagger): "<<t_start<<" end time (c): "<<t_end<<std::endl;
  
  //compute local weight of the removed segment with t_start and t_end
  segment new_segment(t_start, t_end);
  //std::cout<<clred<<"antisegment insert."<<std::endl;
  double local_weight_change=local_config.local_weight_change(new_segment, orbital, true);
  //std::cout<<clred<<"antisegment insert done."<<std::endl;
  
  //compute hybridization weight change
  segment new_antisegment(t_end,t_start);
  double hybridization_weight_change=hyb_config.hyb_weight_change_insert(new_antisegment, orbital);
  
  //compute the proposal probability ratio
  double permutation_factor=beta*t_next_segment_end/(local_config.order(orbital)+1);
  
  //perform metropolis
  double weight_change=local_weight_change*hybridization_weight_change*permutation_factor;
  
  //std::cout<<clred<<"weight change: "<<weight_change<<" l: "<<local_weight_change<<" h: "<<hybridization_weight_change<<" p: "<<permutation_factor<<cblack<<std::endl;
  
  if(std::abs(weight_change)>random()){
    nacc[3]++;
    //std::cout<<cred<<"accepting insert antisegment."<<cblack<<std::endl;
    if(weight_change < 0) sign*=-1.;
    local_config.insert_antisegment(new_antisegment, orbital);
    hyb_config.insert_antisegment(new_antisegment, orbital);
    //std::cout<<cred<<"done accepting insert antisegment."<<cblack<<std::endl;
  }
}

void hybridization::remove_antisegment_update(int orbital){
  nprop[4]++;
  int k=local_config.order(orbital);
  
  if(k==0) return; //no point, this is an empty orbital
  int segment_nr=(int)(random()*k);
  
  //try to merge segment k and segment k+1
  segment segment_earlier=local_config.get_segment(segment_nr, orbital);
  segment segment_later  =local_config.get_segment(segment_nr==k-1?0:segment_nr+1, orbital);
  
  //std::cout<<clcyan<<"antisegment removal update: "<<cblack<<*this<<std::endl;
  //std::cout<<clcyan<<" antisegment start time: (cdagger) "<<segment_earlier.t_end_<<" end time: (c): "<<segment_later.t_start_<<std::endl;
  //compute local weight of the antisegment. note that time direction here has to be forward
  segment segment_forward(segment_earlier.t_end_, segment_later.t_start_);
  //std::cout<<clgreen<<"antisegment remove."<<std::endl;
  double local_weight_change=1./local_config.local_weight_change(segment_forward, orbital, true);
  //std::cout<<clgreen<<"antisegment remove done."<<std::endl;
  
  //compute hybridization weight change
  segment antisegment(segment_later.t_start_, segment_earlier.t_end_);
  double hybridization_weight_change=hyb_config.hyb_weight_change_remove(antisegment, orbital);
  
  //compute the proposal probability ratio
  double t_next_segment_end=local_config.order(orbital)==1?beta:segment_later.t_end_-segment_earlier.t_end_; if(t_next_segment_end<0.) t_next_segment_end+=beta;
  double permutation_factor=local_config.order(orbital)/(beta*t_next_segment_end);
  //perform metropolis
  double weight_change=local_weight_change/hybridization_weight_change*permutation_factor;
  //std::cout<<clblue<<"weight change: "<<weight_change<<" l: "<<local_weight_change<<" h: "<<hybridization_weight_change<<" p: "<<permutation_factor<<cblack<<std::endl;
  
  if(std::abs(weight_change)>random()){
    nacc[4]++;
    //std::cout<<cred<<"accepting remove antisegment."<<cblack<<std::endl;
    if(weight_change < 0) sign*=-1.;
    local_config.remove_antisegment(antisegment, orbital);
    hyb_config.remove_antisegment(antisegment, orbital);
    //std::cout<<cred<<"done accepting remove antisegment."<<cblack<<std::endl;
  }
}

/********************************************************\
 New spin-flip updates
 Idea: take remove_segment_update and insert_segment_update and combine them so one segment is removed from on orbital (spin up or down) and inserted on the corresponding other_orbital (spin down or up), if the other_orbital is not filled.
 \********************************************************/
void hybridization::spin_flip_update(int orbital){
  
  nprop[5]++;
  
  int k=local_config.order(orbital);
  
  if(k==0) return; //no point, this is an empty orbital
  //    if (local_config.zero_order_orbital_occupied(orbital)) return;
  int other_orbital=(int)(random()*n_orbitals);
  int other_orbital_order=local_config.order(other_orbital);
  if (orbital == other_orbital) return;
  if (local_config.zero_order_orbital_occupied(other_orbital)) return;
  int segment_nr=(int)(random()*k);
  segment segment_to_flip=local_config.get_segment(segment_nr, orbital);
  
  
  if (local_config.has_overlap(segment_to_flip,other_orbital)) return;
  
  double t_next_segment_start=local_config.find_next_segment_start_distance(segment_to_flip.t_start_,other_orbital);
  double t_next_segment_end=local_config.find_next_segment_end_distance(segment_to_flip.t_start_,other_orbital);
  double seg_length = segment_to_flip.t_end_ - segment_to_flip.t_start_;
  if (seg_length<0.0) seg_length += beta;
  

  //Totally experimential - mistakes here?
  double t_start = segment_to_flip.t_start_,t_end=segment_to_flip.t_end_;
  if(t_end >= beta) t_end-=beta;
  if(t_end<=t_start || t_end<=0.0) { /*std::cerr<<"rare (impossible?) event: t_start = t_end."<<std::endl; */return; }
  
  //compute local weight change: As we intend to propose a flip, we can
  //safely ignore the intermediate state and directly compare the energies
  //of the two states involved
  // This energy is associated with the present segment. We give it a negative weight
  // because it is to be removed
  double de = -local_config.local_energy(segment_to_flip,orbital);
  //double full_weight_orig=full_weight();
  double hybridization_weight_change_1=1./hyb_config.hyb_weight_change_remove(segment_to_flip, orbital); // from line 187 - remove_segment_update
  double weight_change_1=hybridization_weight_change_1;
  if (weight_change_1<0) sign*=-1;
  local_config.remove_segment(segment_to_flip, orbital);
  hyb_config.remove_segment(segment_to_flip, orbital);
  //double full_weight_removed=full_weight();
  segment new_segment(t_start,t_end);
  de += local_config.local_energy(new_segment,other_orbital);
  double hybridization_weight_change_2=hyb_config.hyb_weight_change_insert(new_segment, other_orbital);   double weight_change_2=std::exp(de)*hybridization_weight_change_2;

  //the permutation factor has the probability of proposing this move (the old order in this orbital) divided by the probability of proposing the reverse move (the new order in the new orbital)
  double permutation_factor=k/(double)(other_orbital_order+1.);

  if(std::abs(weight_change_1*weight_change_2*permutation_factor)>random()){ //Accepted
    nacc[5]++;
    if(weight_change_2 < 0) sign*=-1.;
    local_config.insert_segment(new_segment, other_orbital);
    hyb_config.insert_segment(new_segment, other_orbital);
  } else { //Not accepted, thus restore old configuration
    double wc = hyb_config.hyb_weight_change_insert(new_segment,orbital);
    if (wc<0.0) sign *= -1;
    local_config.insert_segment(new_segment, orbital);
    hyb_config.insert_segment(new_segment, orbital);
  }
}
