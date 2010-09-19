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

void WRun::create_observables()
{
  create_common_observables();
  
  measurements << alps::RealVectorObservable("Statistics");
  measurements << alps::make_observable(alps::RealVectorObservable("Winding number^2"), is_signed_);

  if(use_1D_stiffness){   //@#$br: take care of D=1 ---> need a histogram of the winding numbers
     measurements << alps::make_observable(alps::RealVectorObservable("Winding number histogram"), is_signed_); 
  }
  
  // warn a user if using the D>1 stiffness estimator in 1D  
  //if(1==dimension() && !use_1D_stiffness){
  //  std::cout<<" *** Warning: using the <W^2> estimator for the stiffness \n *** Consider switching USE_1D_STIFFNESS=true \n" ;
  //}


 
  if (nonlocal) {
    measurements << alps::RealObservable("Ratio time intervals");
    measurements << alps::RealObservable("Ratio time intervals 2");
    measurements << alps::RealObservable("Time interval");
    measurements << alps::RealVectorObservable("Statistics time intervals");
  }
  
  site_iterator it,end;
  num_chains=0;
  if (chain_kappa) {
    measurements << alps::make_observable(alps::RealVectorObservable("Chain density"), is_signed_);
    measurements << alps::make_observable(alps::RealVectorObservable("Chain density^2"), is_signed_);
    if (!where.empty()) {
      for (boost::tie(it,end)=sites();it!=end;++it) 
        if (coordinate(*it).size()>=2) {
          chain_number[*it]=int(coordinate(*it)[1]);
          if (coordinate(*it)[1]>=num_chains)
            num_chains=int(coordinate(*it)[1]+1);
        }
    }
  }
}

int WRun::get_particle_number() {
  double n=0.;
  for (int i=0;i<num_sites();i++)
    n += QMCRun<>::diagonal_matrix_element["n"][site_type(i)][initial_state(i)];
  return static_cast<int>(n);
}

void WRun::adjustment() {

  std::cerr << "Starting adjustment\n";
  int number_of_bosons=static_cast<int>(parms["NUMBER_OF_PARTICLES"]);
  double correction=static_cast<double>(parms["CORRECTION"]);

  if ((corrections_upwards==0||corrections_downwards==0)&&steps%25==0) {
    // coarse adjustment
    if (number_of_bosons>get_particle_number()) {
      parms[adjust_parameter] = static_cast<double>(parms[adjust_parameter]) + 10.*correction;
      corrections_upwards++;
    }
    else {
      parms[adjust_parameter] = static_cast<double>(parms[adjust_parameter]) - 10.*correction;
      corrections_downwards++;
    }
    std::cerr << "Init\n";
    initialize_hamiltonian();
  }
  else
    if (!preadjustment_done&&steps%25==0)
      preadjustment_done = true;
  
  if (!adjustment_done&&preadjustment_done) {
    // fine adjustment
    nob.push_back(get_particle_number());
    if (steps%60==0) {
      double n=0.;
      // averaging of measured number of bosons over last 50 steps
      for (int i=10;i<nob.size();i++)
        n += nob[i];
      n /= (nob.size()-10);
      nob.clear();
      // check if average is close enough to number_of_bosons
      if (fabs(number_of_bosons-n)<number_of_bosons/200.) {
        adjustment_done = true;
        std::cout << "Adjusted parameter is " << static_cast<double>(parms[adjust_parameter]) << ".\n";
        std::cout << "Adjustment took " << steps << " steps.\n";
        if (steps>thermal_sweeps/4.) {
          std::cout << "Adjustment took longer than 25 % of thermalization phase.\n";
          exit(1);
        }
      }
      else {
        if (number_of_bosons>n)
          parms[adjust_parameter] = static_cast<double>(parms[adjust_parameter]) + correction;
        else
          parms[adjust_parameter] = static_cast<double>(parms[adjust_parameter]) - correction;
    std::cerr << "Init\n";
        initialize_hamiltonian();
      }
    }
  }
  std::cerr << "Done adjustment\n";
}

void WRun::make_meas()
{
  std::vector<std::vector<double> > const& matrix_element_n = QMCRun<>::diagonal_matrix_element["n"];

  // measure energy and particle numbers
  double energy=0.;
  std::vector<state_type> local(num_sites());
  std::valarray<double> chain_n(num_chains);

  for (int i=0;i<num_sites();i++) {
    // determine energy
    state_type s=initial_state(i);
    if (nonlocal)
      energy += H0(site(i)) - kinks[i].size()/(2.*beta);
    else
      energy += onsite_energy(s,i) -kinks[i].size()/(2.*beta);

    local[i]=s;
    if (chain_kappa)
      chain_n[chain_number[i]]+=matrix_element_n[site_type(i)][s];
  }

  if (!do_common_measurements(Sign,local)) {
    green = 0.;
    return;
  }

  // measure winding numbers
  bond_iterator bi, bi_end;

  std::vector<int> bond_dir(dimension());
  std::valarray<double> winding_number(0., dimension());
  for(boost::tie(bi, bi_end) = bonds(); bi != bi_end; ++bi) {
    int s1 = source(*bi);
    int s2 = target(*bi);
    vector_type v=bond_vector_relative(*bi);

    for(int d=0; d<dimension(); ++d)
      bond_dir[d] = boundary_crossing[*bi].crosses(d);
    
    cyclic_iterator i1 = first_kink(s1);
    cyclic_iterator i2 = first_kink(s2);
    if (i1.valid() && i2.valid()) {
      while(true) {
        if (i1->time() == i2->time()) {
          int hop_dir = (i1->state() > (i1-1)->state() ? -1 : 1);
          for(int d=0; d<dimension(); ++d) {
            winding_number[d]+=hop_dir*v[d];
          }
        }
        if (i1->time()<i2->time()) {
          ++i1;
          if (i1==first_kink(s1))
            break;
        }
        else {
          ++i2;
          if (i2==first_kink(s2))
            break;
        } 
      }
    }
  }
  std::valarray<double> winding_number2(0., dimension());
  winding_number2 = winding_number*winding_number;

  if(use_1D_stiffness){
     // @#$br : record the winding numbers = -1, 0, 1
     //         the vector entries are 0,1,2, respectively            
     std::valarray<double> winding_numberz(3);
     for(int i=0; i<3;++i){  winding_numberz[i]=0. ;}
     if( fabs( winding_number[0] ) <=1 ) {
          winding_numberz[(int)winding_number[0]+1]+=Sign ; 
          measurements["Winding number histogram"]  << winding_numberz ; //*Sign; 
     }
  }
  else{
     // determine stiffness
     double stiffness=0.;
     for(int d=0; d<dimension(); ++d) {
        stiffness += winding_number[d]*winding_number[d];
     }
     measurements["Stiffness"] << stiffness*Sign/(dimension()*beta);
  }

  double vol = num_sites();

  measurements["Statistics"]  << stat; 
  winding_number2 *= Sign;
  measurements["Winding number^2"] << winding_number2;
  measurements["Energy"] << energy*Sign;
  measurements["Energy Density"] << energy/vol*Sign;
  //measurements["Stiffness"] << stiffness*Sign/(dimension()*beta);

  if (measure_green_function_) {
    green /= 2.*double(skip_measurements);
    green *= Sign;
    measurements["Green's Function"] << green;
    green=0.;
  }
  if(chain_kappa) {
    chain_n/=double(num_sites()/num_chains);
    
    std::valarray<double> chain_n2(0., dimension());
    chain_n2 = chain_n2*chain_n2;

    chain_n  *= Sign;
    chain_n2 *= Sign;

    measurements["Chain density"] << chain_n;
    measurements["Chain density^2"] << chain_n2;
  }

}   // WRun::make_meas


void WRun::measure_green()
{
/*
  int delta=worm_head[1].site()-worm_head[0].site();
  if (delta<0)
    delta +=num_sites();
  for (int h=0;h<(delta==0 ? 2 :1);++h) {
    cyclic_iterator head=worm_head[1-h].kink();
    cyclic_iterator start=head-1;
    cyclic_iterator end=head+1;
    site_descriptor s=worm_head[1-h].site();
    bool do_create = (head->state() < start->state());

    if ((end->time()-start->time())>=(end->time()-worm_head[h]->time())) {
      // needs to be fixed ...
      double statistical_weight= (onsite_energy(start->state(),s)-onsite_energy(head->state(),s));
      double time_interval=end->time() - start->time();
      double time_delta = worm_head[h]->time() - start->time();
      if(time_delta==1.) 
        time_delta=0.;
      double P=exp(-beta*time_delta*statistical_weight)/integrated_weight(statistical_weight, beta*time_interval);
      green[(delta!=0 || do_create) ? delta+1 : delta]+=P;    
    }
  }
  */
}   // WRun::measure_green

