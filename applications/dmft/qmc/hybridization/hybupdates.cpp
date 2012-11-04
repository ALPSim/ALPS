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

#include"hyb.hpp"
#include"hyblocal.hpp"

//this is the heard of the Monte Carlo procedure: we have the following updates:
//1: change the zero order state, swap an empty orbital versus a filled one.
//2: shift an existing segment start- or end-point
//3: insert a new segment
//4: insert or remove an anti-segment.
//see our review for details of these updates.
void hybridization::update(){
  for(std::size_t n=0; n<N_accu; ++n){
    for(std::size_t i=0;i<N_meas;++i){
      double update_type=random();
      if(update_type < 0.1){
        change_zero_order_state_update();
      }else if(update_type < 0){
        shift_segment_update();
      }else if(update_type < 0.6){
        insert_remove_segment_update();
      }else{
        insert_remove_antisegment_update();
      }
      sweeps++;

      //these are cheap measurements that should be done every time.
      if(is_thermalized()){
        measure_order();
        measure_G();
        meas_count++;
      }
    }//N_meas
    if(is_thermalized()){
      std::vector<std::map<double,double> > F_prefactor;
      if(MEASURE_freq || MEASURE_legendre || MEASURE_g2w || MEASURE_h2w)
        local_config.get_F_prefactor(F_prefactor);//compute segment overlaps in local config

      measure_Gw(F_prefactor);
      measure_Gl(F_prefactor);
      measure_sector_statistics();

      //measure 2-particle quantities
      measure_nn();
      measure_nnt();
      measure_nnw();

      if(MEASURE_g2w  || MEASURE_h2w) measure_G2w(F_prefactor);

      accu_count++;
    }
  }//N_accu
}

void hybridization::change_zero_order_state_update(){
  //choose the orbital in which we do the update
  int orbital=(int)(random()*n_orbitals);

  //changing the zero order state only makes sense if we are at zero order.
  if(!local_config.order(orbital)==0) return;
  
  //propose to change orbital from occuppied to unoccuppied.
  if(local_config.zero_order_orbital_occupied(orbital)){
    double local_weight_change=1./local_config.local_weight_change(segment(0,beta), orbital, false);
    if(std::abs(local_weight_change)>random()){
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
void hybridization::insert_segment_update(int orbital){
  //std::cout<<clred<<"starting insertion update."<<cblack<<std::endl;
  if(local_config.order(orbital)==0 && local_config.zero_order_orbital_occupied(orbital)) return; //can't insert segment, orbital is fully occuppied.
  double t_start=random()*beta; //start time of a segment
  if(local_config.exists(t_start)){ std::cerr<<"rare event, duplicate: "<<t_start<<std::endl; return;} //time already exists.
  double t_next_segment_start=local_config.find_next_segment_start_distance(t_start,orbital);
  double t_next_segment_end=local_config.find_next_segment_end_distance(t_start,orbital);
  //std::cout<<"============================================"<<cblack<<std::endl;
  //std::cout<<clblue<<"orbital: "<<orbital<<" time is: "<<t_start<<" segment start distance: "<<t_next_segment_start<<" end distance: "<<t_next_segment_end<<cblack<<std::endl;
  //std::cout<<*this<<std::endl;
  //std::cout<<"============================================"<<cblack<<std::endl;
  if(t_next_segment_end < t_next_segment_start) return; //we're trying to create a segment on top of another segment. abort.
  
  //draw an end time
  double t_end=t_start+random()*t_next_segment_start;
  if(t_end > beta) t_end-=beta;
  if(t_start == t_end || local_config.exists(t_end)){ std::cerr<<"rare event, duplicate: "<<t_end<<" "<<t_start<<std::endl; return;} //time already exists.
  
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
    if(weight_change < 0) sign*=-1.;
    local_config.insert_segment(new_segment, orbital);
    hyb_config.insert_segment(new_segment, orbital);
  }
}
void hybridization::remove_segment_update(int orbital){
  //std::cout<<clblue<<"starting removal update."<<cblack<<std::endl;
  int k=local_config.order(orbital);

  if(k==0) return; //no point, this is an empty orbital

  int segment_nr=(int)(random()*k);
  
  segment segment_to_remove=local_config.get_segment(segment_nr, orbital);

  double local_weight_change=1./local_config.local_weight_change(segment_to_remove, orbital, false);
  
  //compute hybridization weight change
  double hybridization_weight_change=1./hyb_config.hyb_weight_change_remove(segment_to_remove, orbital);
  
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
    if(weight_change < 0) sign*=-1.;
    local_config.remove_segment(segment_to_remove, orbital);
    hyb_config.remove_segment(segment_to_remove, orbital);
  }
}
void hybridization::insert_antisegment_update(int orbital){
  if(local_config.order(orbital)==0 && !local_config.zero_order_orbital_occupied(orbital)) return; //can't insert an antisegment, orbital is empty.
  double t_start=random()*beta; //start time of the anti segment
  if(local_config.exists(t_start)){ std::cerr<<"rare event, duplicate: "<<t_start<<std::endl; return;} //time already exists.
  double t_next_segment_start=local_config.find_next_segment_start_distance(t_start,orbital);
  double t_next_segment_end=local_config.find_next_segment_end_distance(t_start,orbital);

  if(t_next_segment_start < t_next_segment_end) return; //we're trying to create an antisegment where there is no segment abort.
  
  //draw an end time
  double t_end=t_start+random()*t_next_segment_end;
  if(t_end > beta) t_end-=beta;
  if(t_start == t_end || local_config.exists(t_end)){ std::cerr<<"rare event, duplicate: "<<t_end<<" "<<t_start<<std::endl; return;} //time already exists.
  
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
    //std::cout<<cred<<"accepting insert antisegment."<<cblack<<std::endl;
    if(weight_change < 0) sign*=-1.;
    local_config.insert_antisegment(new_antisegment, orbital);
    hyb_config.insert_antisegment(new_antisegment, orbital);
    //std::cout<<cred<<"done accepting insert antisegment."<<cblack<<std::endl;
  }
}
void hybridization::remove_antisegment_update(int orbital){
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
    //std::cout<<cred<<"accepting remove antisegment."<<cblack<<std::endl;
    if(weight_change < 0) sign*=-1.;
    try{
    local_config.remove_antisegment(antisegment, orbital);
    hyb_config.remove_antisegment(antisegment, orbital);
    }catch(std::exception &e){
      std::cerr<<cred<<"caught exception in remove antisegment!"<<cblack<<std::endl;
      std::cout<<" orbital is: "<<orbital<<std::endl;
      std::cout<<" segment no is: "<<segment_nr<<std::endl;
      std::cout<<" segment earlier is: "<<segment_earlier<<std::endl;
      std::cout<<" segment later is: "<<segment_later<<std::endl;
      std::cout<<" segment forward is: "<<segment_forward<<std::endl;
      std::cout<<" antisegment is: "<<antisegment<<std::endl;
      std::cout<<" t next segment end: "<<t_next_segment_end<<std::endl;
      exit(1);
    }
    //std::cout<<cred<<"done accepting remove antisegment."<<cblack<<std::endl;
  }
}
