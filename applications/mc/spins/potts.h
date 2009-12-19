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

#ifndef ALPS_APPLICATIONS_MC_SPIN_POTTS_H
#define ALPS_APPLICATIONS_MC_SPIN_POTTS_H

#include <boost/random/uniform_smallint.hpp>
#include "matrices.h"
#include "tinyvec.h"
#include <iostream>

template <unsigned int Q>
class PottsMoment {
  public:
    static const int dim = 2;
  
    typedef double magnetization_type;
    // all states initialized to 'color' 0
        PottsMoment() : state_(0) {}
        // updates are flips to a new color
    typedef unsigned int update_type;
    template <class RNG> update_type random_update(RNG& rng) 
    { // enforce a different color than the current one
      update_type new_color=int((Q-1)*rng());
      return (new_color>=state_ ? new_color+1 : new_color); 
    }
    void update(update_type new_color) {state_=new_color;}
    double project(typename PottsMoment<Q>::update_type) {return 0.;}
    // comparisons needed for energy calculation
    bool operator==(PottsMoment<Q> m) const {return state_==m.state_;}
    bool operator==(update_type new_color) const {return state_==new_color;}
  friend alps::ODump& operator << (alps::ODump& dump, const PottsMoment<Q>& m) 
  { return dump << m.state_;}

  friend std::ostream& operator << (std::ostream& out, const PottsMoment<Q>& m)
  { return out << m.state_; }

  friend alps::IDump& operator >> (alps::IDump& dump, PottsMoment<Q>& m) 
  { return dump >> m.state_;}

  magnetization_type magnetization() const { return state_==0 ? 1. : 0.;}

  double mag_h(TinyVector<double, 2> h) const {
    double length = h.get_length2();
    return (length==0 ? 0. : (state_==0 ? 1. : 0.));
  }
  
  private:
          // the moment can take any 'color' from 0 ... Q-1
        unsigned int state_;        
};

/**
 * Computes the change in the bond energy between m1 and m2 when the state
 * m1 is changed to new_color.
 * 
 * \param J the couplings matrix.
 */
template <unsigned int Q>
inline double bond_energy_change(PottsMoment<Q>& m1, PottsMoment<Q>& m2,
              MIdMatrix<double,2> J, 
              typename PottsMoment<Q>::update_type& new_color) {
  return (m1==m2 ? (J.getElement_nc(0,0)) : 
               ( m2 == new_color ? -(J.getElement_nc(0,0)) : 0.));
};

/**
 * Computes the bond energy between m1 and m2.
 *
 * \param J the couplings matrix
 */
template <unsigned int Q>
inline double bond_energy(PottsMoment<Q> m1, PottsMoment<Q> m2,
              MIdMatrix<double,2>& J) {
  return m1==m2 ? -(J.getElement_nc(0,0)) : 0.;
};

/**
 * Computes the change in the site energy, when the state m is changed to new_color.
 *
 * \param h the external magnetic field
 */
template <unsigned int Q>
inline double site_energy_change(PottsMoment<Q> m, TinyVector<double,2>& h_,
              typename PottsMoment<Q>::update_type new_color) {
  return m==0 ? h_[1] : (new_color==0 ? -h_[1] : 0.);
};

/**
 * Returns the magnetic energy of the current configuration in the external
 * magnetic field.
 * \param h the external magnetic field.
 */
template <unsigned int Q>
inline double site_energy(PottsMoment<Q> m, TinyVector<double,2>& h_) {
  return m==0 ? -h_[1] : 0.;
};

// guess that there is no onsite_energy term for the potts-modell
/**
 * Dummy function for the onsite_energy of the potts model.
 * always return 0.
 */
template<unsigned int Q>
inline double onsite_energy(PottsMoment<Q>, MIdMatrix<double,2>) {
  return 0.0;
}

/**
 * Dummy function for the onsite_energy_change, returns always 0
 */
template<unsigned int Q>
double onsite_energy_change(PottsMoment<Q>, MIdMatrix<double,2>,
                         typename PottsMoment<Q>::update_type) {
  return 0.0;
}

#endif
