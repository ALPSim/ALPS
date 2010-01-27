/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2005 by Fabien Alet <alet@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
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
