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

/* $Id: SSE.Measurements.cpp,v 1.33 2006/08/08 09:17:27 troyer Exp $ */

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
 
  create_common_observables();
  measurements << SimpleRealVectorObservable("logg");
  measurements << SimpleRealVectorObservable("Histogram");
  measurements << SimpleRealVectorObservable("HistoUp");
  measurements << SimpleRealObservable("Total Measurements");
  measurements << RealObservable("Time Up");
  measurements << RealObservable("Time Down");
  measurements << RealObservable("RealTime Up");
  measurements << RealObservable("RealTime Down");

  // Add the physical observables you want to measure here...  
 
}
   
void SSE::finish_measurements() {
    // In case you need to clean something up...
}

void SSE::do_measurements()
{
// We could measure something if we wanted to...
}
