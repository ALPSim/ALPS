/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1999-2003 by Matthias Troyer <troyer@comp-phys.org>
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

/* $Id$ */

#ifndef ALPS_APPLICATIONS_MC_SPIN_CONNECT_H
#define ALPS_APPLICATIONS_MC_SPIN_CONNECT_H

#include "ising.h"
#include "potts.h"
#include "on.h"
#include <boost/random/uniform_real.hpp>

// spinfactors not needed to decide whether two nodes are connected ... ??
template <class Moment, class RNG>
class Connector {
  public:
    // store inverrse temperature and refernce to random number generator 
    Connector(double beta, RNG& rng) 
      : beta_(beta), random_(rng) {}
    
    // calculate connection probability and decide randomly if less than one
    bool operator()(const Moment& m1, const Moment& m2, 
                    const typename Moment::update_type& update, 
                    const MIdMatrix<double, Moment::dim> J) 
    {
      double delta_E = bond_energy_change(m1,m2,J,update);
      return (delta_E> 0 && random_() < (1.-exp(-beta_*delta_E)));
    }

    bool operator()(const Moment& m1, const Moment& m2,
                    const typename Moment::update_type& update,
                    const MIdMatrix<double,Moment::dim> J,
                    const MIdMatrix<double,Moment::dim> D)
    {
      double delta_E = bond_energy_change(m1,m2,J,update);
      delta_E += onsite_energy_change(m1,D,update);
      return (delta_E>0 && random_() <(1.-exp(-beta_*delta_E)));
    }
    
  private:
    double beta_;
    RNG& random_;
};

// specialized implementation for Ising model
template <class RNG> 
class Connector< IsingMoment, RNG>
{
  public:
        // precalculate probability in constructor
    Connector(double beta, RNG& rng) 
      : beta_(beta),probability_(1.-exp(-2.*beta)), random_(rng) {}
    
    bool operator()(const IsingMoment& m1, const IsingMoment& m2, 
                    const IsingMoment::update_type&) 
    { return (m1==m2 && random_() < probability_);}
    
    bool operator()(const IsingMoment& m1, const IsingMoment& m2, 
                    const IsingMoment::update_type&, 
                    MIdMatrix<double,1> Jmat) 
    { 
      double J = Jmat.det();
      if (J>0.)
        return (m1== m2 && random_() <  ( J == 1. ? probability_ : (1.-exp(-2.*beta_*J))));
      else
        return (m1 != m2 && random_() < ( J == -1. ? probability_ : (1.-exp(+2.*beta_*J))));
    }

    bool operator()(const IsingMoment& m1, const IsingMoment& m2,
                    const IsingMoment::update_type& t,
                    MIdMatrix<double,1> Jmat, MIdMatrix<double,1>) {
      /* does the onsite-term have any influence ?? I don't think so!! */
      return operator()(m1,m2,t,Jmat);
    }     
     
  private:
    double beta_;
    double probability_;
    RNG& random_;
};

// specialized implementation for Potts models
template <unsigned int Q, class RNG> 
class Connector<PottsMoment<Q>,RNG>
{
  public:
        // precalculate probability in constructor
    Connector(double beta, RNG& rng) 
     : beta_(beta),probability_(1.-exp(-beta)), random_(rng) {}
    
    bool operator()(const PottsMoment<Q>& m1, const PottsMoment<Q>& m2, 
                    const typename PottsMoment<Q>::update_type&) 
    { return (m1==m2 && random_() < probability_);}
    
    bool operator()(const PottsMoment<Q>& m1, const PottsMoment<Q>& m2, 
                    const typename PottsMoment<Q>::update_type&, 
                    const MIdMatrix<double,2> Jmat) { 
    double J = Jmat.getElement_nc(0,0);
    return (m1==m2 && random_() < ( J == 1. ? probability_ : (1.-exp(-beta_*J))));
    }

    bool operator()(const PottsMoment<Q>& m1, const PottsMoment<Q>& m2, 
                    const typename PottsMoment<Q>::update_type& update, 
                    const MIdMatrix<double,2> Jmat,
                    const MIdMatrix<double,2>) { 
      /* what influence does the D-term have in this model?? */
      return operator()(m1,m2,update,Jmat);
    }

  private:
    double beta_;
    double probability_;
    RNG& random_;
};

#endif
