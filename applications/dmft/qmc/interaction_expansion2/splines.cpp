/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/


#include "interaction_expansion.hpp"

///as the name says: linearly interpolate between two points. this is VERY
//EXPENSIVE because of the division.

template<class X, class Y> inline Y linear_interpolate(const X x0, const X x1, const Y y0, const Y y1, const X x)
{
  return y0 + (x-x0)/(x1-x0)*(y1-y0);
}

///Compute the Green's function G0 (the BARE) green's function between two points
double InteractionExpansion::green0_spline(const creator &cdagger, const annihilator &c) const
{
  //for the M matrix we need the bare temporal Green's function. 
  //For the dressed frequency GF we need the bare frequency GF.
  //here we get the bare temporal GF, interpolated, as needed for the M matrix.
  //we will receive this as an input into our solver later.
  spin_t flavor;
  if((flavor=cdagger.flavor()) != c.flavor()){
    return 0; //the delta in spin space.
  }
  itime_t delta_t=cdagger.t()-c.t();
  site_t site1 = cdagger.s();
  site_t site2 = c.s();
  return green0_spline(delta_t, flavor, site1, site2);  
}


///Compute the bare green's function for a given flavor, site, and imaginary time.
double InteractionExpansion::green0_spline(const itime_t delta_t, const spin_t flavor, const site_t site1, const site_t site2) const
{
  if(delta_t*delta_t < almost_zero){
    return bare_green_itime(0,site1, site2, flavor);
  } 
  else if(delta_t>0){
    int time_index_1 = (int)(delta_t*n_tau*temperature);
    int time_index_2 = time_index_1+1;
    return linear_interpolate((double)time_index_1*beta*n_tau_inv, (double)time_index_2*beta*n_tau_inv,
                              bare_green_itime(time_index_1,site1, site2, flavor),
                              bare_green_itime(time_index_2,site1, site2, flavor),delta_t);
  } 
  else{
    int time_index_1 = (int)(delta_t*n_tau*temperature+n_tau);
    int time_index_2 = time_index_1+1;
    return -linear_interpolate((double)time_index_1*beta*n_tau_inv, (double)time_index_2*beta*n_tau_inv,
                               bare_green_itime(time_index_1,site1,site2,flavor),
                               bare_green_itime(time_index_2,site1,site2,flavor),delta_t+beta);
  }
}


double InteractionExpansion::green0_spline(const itime_t delta_t, const spin_t flavor) const
{
  if(delta_t*delta_t < almost_zero){
    return bare_green_itime(0, flavor);
  } 
  else if(delta_t>0){
    int time_index_1 = (int)(delta_t*n_tau*temperature);
    int time_index_2 = time_index_1+1;
    return linear_interpolate((double)time_index_1*beta*n_tau_inv, (double)time_index_2*beta*n_tau_inv,
                              bare_green_itime(time_index_1,flavor),
                              bare_green_itime(time_index_2,flavor),delta_t);
  } 
  else{
    int time_index_1 = (int)(delta_t*n_tau*temperature+n_tau);
    int time_index_2 = time_index_1+1;
    return -linear_interpolate((double)time_index_1*beta*n_tau_inv, (double)time_index_2*beta*n_tau_inv, //time ordering
                               bare_green_itime(time_index_1,flavor),
                               bare_green_itime(time_index_2,flavor),delta_t+beta);
  }
}

