/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2005 by Fabien Alet <alet@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
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
