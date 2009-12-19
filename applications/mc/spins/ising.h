/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1999-2003 by Matthias Troyer <troyer@comp-phys.org>,
*                            Fabian Stoeckli <fabstoec@student.ethz.ch>
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

#ifndef ALPS_APPLICATIONS_MC_SPIN_H
#define ALPS_APPLICATIONS_MC_SPIN_H

#include "tinyvec.h"
#include "matrices.h"
#include <iostream>

class IsingMoment {
  public:
    // initialize all moments up (true)
    typedef double magnetization_type;
    IsingMoment() : state_(true) {}
    // there is only one possible update type: a spin flip
    enum update_type {flip}; 
    template <class RNG> update_type random_update(RNG&) {return flip;}
    // the update flips(inverts) the spin
    void update(update_type) {state_=!state_;}
    double project(IsingMoment::update_type) {return 1.;}
    // equality comparison is used by energy calculation
    bool operator==(IsingMoment m) const {return state_==m.state_;}
    bool operator!=(IsingMoment m) const {return state_!=m.state_;}
  friend alps::ODump& operator << (alps::ODump& dump, const IsingMoment& m) 
  { return dump << m.state_;}

  friend alps::IDump& operator >> (alps::IDump& dump, IsingMoment& m) 
  { return dump >> m.state_;}
  magnetization_type magnetization() const { return state_ ? 1. : -1.;}
  double mag_h(TinyVector<double,1> h_) const
  { return (h_[0] == 0. ? 0. : (state_ ? 1. : -1. )); }

  static const int dim = 1;
  private:
          // moments stored as a boolean: true=up, false=down
    bool state_; 

  friend inline std::ostream& operator<<(std::ostream& os, IsingMoment const& m)
  {
    return os << (m.state_ ? 1 : -1);
  }


};

inline double bond_energy(IsingMoment m1, IsingMoment m2,
                         MIdMatrix<double,1>& J)
{ return m1==m2 ? -(J.getElement_nc(0,0)) : (J.getElement_nc(0,0)); }

inline double bond_energy_change(IsingMoment& m1, IsingMoment& m2, 
                         MIdMatrix<double,1> J,
                         IsingMoment::update_type)
{ return m1==m2 ? 2.*(J.getElement_nc(0,0)) : -2.*(J.getElement_nc(0,0)); }

inline double site_energy(IsingMoment m, TinyVector<double,1>& h)
{ return m.magnetization()==1. ? -h[0]: h[0]; }

inline double site_energy_change(IsingMoment m, TinyVector<double,1>& h, 
                         IsingMoment::update_type) 
{ return m.magnetization()==1. ? 2.*h[0] : -2.*h[0]; }

inline double onsite_energy(IsingMoment, MIdMatrix<double,1>& D) {
  return D.getElement_nc(0,0);
}

inline double onsite_energy_change(IsingMoment, MIdMatrix<double,1>&,
                         IsingMoment::update_type ) {
// the onsite_energy does not depend on the spin orientation (in this model!!)
// thus, for a given D, there is no change in the onsite_energy, when the 
// spins are flipped.
  return 0.0;
}


#endif
