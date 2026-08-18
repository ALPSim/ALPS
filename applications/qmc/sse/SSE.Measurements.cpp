/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2005 by Fabien Alet <alet@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include "SSE.hpp"


double SSE::get_sign()
{
  double sign=1.;
  
  for ( vector<vertex_type>::iterator it=operator_string.begin();it!=operator_string.end();++it) { 
    if ((it->non_diagonal())) { // Off-Diagonal
      state_type* MP= it->leg;
      if (matrix_sign[bond_type[bond(it->bond_number)]][MP[0]][MP[1]][MP[2]][MP[3]]) { sign=-sign; /*cout << "sign becomes " << sign << endl;*/}
    }
  }
  //cout << "Sign is finally *** " << sign << endl;
  return sign;
}

void SSE::create_observables()
{
  initialize_simulation();

  create_common_observables();
  
  if (measure_green_function_) {
    green.resize(measurement_origin_ ? num_sites() : num_distances());
    green=0.;
  }
    
  if (super_type::measure_site_compressibility_) {
    localint.resize(num_sites());
    localint2.resize(num_sites());
  }

  measurements << alps::make_observable(RealObservable("n"),is_signed_);
  measurements << alps::make_observable(RealObservable("n^2"),is_signed_);
  measurements << alps::make_observable(RealObservable("n^3"),is_signed_);
#ifdef AUTO_TIMES
  measurements << alps::RealObservable("Worm Size");
#endif
}
   

void SSE::do_measurements()
{
  double NbSites=num_sites();

  double sign=1.;

  if (is_signed_)
    sign = get_sign();

  measurements["Energy"] << (-1./beta*current_number_of_non_identity+total_shift)*sign;
  measurements["Energy Density"] << (-1./beta*current_number_of_non_identity+total_shift)/NbSites*sign;

  localint = (localint2+localint*localint)*beta/double(operator_string.size())
              /double(operator_string.size()+1);
  do_common_measurements(sign,site_state,localint);
    
  // Green function measurement
  
  if (measure_green_function_) {
    if (measurement_origin_)
      for (int i=0;i<green.size();++i)
        green[i]*=double(num_sites())/double(number_of_worms_per_sweep);
    else
      for (int i=0;i<green.size();++i)
        green[i]*=double(num_sites())/double(number_of_worms_per_sweep*distance_mult[i]);
    measurements["Green's Function"] << green;
    green=0.;
  }
  
  //************************ Stiffness measurement ******************** 

  std::vector<int> bond_dir(dimension());
  std::valarray<double> winding_number(0., dimension());
  boost::multi_array<double,2> bond_type_windings_;
  if (measure_bond_type_stiffness_)
    bond_type_windings_.resize(boost::extents[num_bond_types_][dimension()]);

  for ( vector<vertex_type>::iterator it=operator_string.begin();it!=operator_string.end();++it)
    if(it->non_diagonal()) { 
      bond_descriptor mybond=bond(it->bond_number);              
      int hopping_sign = ((it->leg[0] < it->leg[2])  ? 1 : -1);
      vector_type v=bond_vector_relative(mybond);
      unsigned int bt = original_bond_type[bond_type(mybond)];
      for (int i=0;i<v.size() && i<dimension();++i) {
        winding_number[i]+=hopping_sign*v[i];
        if (measure_bond_type_stiffness_)
          bond_type_windings_[bt][i]+=hopping_sign*v[i];
      }
    }
  
  double winding=0.;
  for(int d=0; d<dimension(); ++d) { 
    winding+=winding_number[d]*winding_number[d]; 
  }
  measurements["Stiffness"] << winding/(beta*dimension())*sign;
  if (measure_bond_type_stiffness_) {
    std::valarray<double> bt_stiffness(0.,num_bond_types_);
    std::valarray<double> bt_stiffness_corr(0.,num_bond_types_*num_bond_types_);
    for(int d=0; d<dimension(); ++d)
      for (std::size_t bt = 0; bt < num_bond_types_; ++bt) {
        bt_stiffness[bt]+=bond_type_windings_[bt][d]*bond_type_windings_[bt][d]; 
        for (std::size_t bt2 = 0; bt2 < num_bond_types_; ++bt2)
          bt_stiffness_corr[bt+num_bond_types_*bt2]+=bond_type_windings_[bt][d]*bond_type_windings_[bt2][d]; 
      }
    bt_stiffness *= sign/(beta*dimension());
    bt_stiffness_corr *= sign/(beta*dimension());
    measurements["Bond Type Stiffness"] << bt_stiffness;
    measurements["Bond Type Stiffness Correlations"] << bt_stiffness_corr;
  }

  measurements["n"] << (double) current_number_of_non_identity*sign;
  measurements["n^2"] << (double) current_number_of_non_identity*current_number_of_non_identity*sign;
  measurements["n^3"] << (double) current_number_of_non_identity*current_number_of_non_identity*current_number_of_non_identity*sign;
  
#ifdef AUTO_TIMES
  measurements["WormSize"] << worm_size;
#endif
}
