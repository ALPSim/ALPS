/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1999-2003 by Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_APPLICATIONS_MC_SPIN_XY_H_
#define ALPS_APPLICATIONS_MC_SPIN_XY_H_
#include <cmath>

#define   TWOPI 2.*M_PI

class XYMoment {
public:
  typedef double update_type;

  XYMoment() : state_(0) { }
  template <class RNG> update_type random_update(RNG& rng) {return boost::uniform_real<>(0,TWOPI)(rng);}
/*
  void prepare(update_type dir) {
#ifdef CSE
    scos_ = 2.*std::cos(state_-dir);
#endif
  }
  */
  void update(update_type dir) {state_ = mod2pi(mod2pi(2*dir-state_)-M_PI);}
  friend double energy_change(XYMoment s1,XYMoment s2,update_type dir)
  { 
#ifdef CSE  
  return scos_*std::cos(s2.state_-dir);}
#else
  return 2.*(std::cos(s1.state_-dir)*std::cos(s2.state_-dir));}
#endif
  friend alps::ODump& operator << (alps::ODump& dump, const XYMoment& m) 
  { return dump << m.state_;}
  friend alps::IDump& operator >> (alps::IDump& dump, XYMoment& m) 
  { return dump >> m.state_;}
private:
  double state_;
  double mod2pi(double x) { return (x<0 ? x+TWOPI : (x>= TWOPI ? x-TWOPI : x)); }
#ifdef CSE
  static double scos_; 
#endif
};

#ifdef CSE
double XYMoment::scos_;
#endif

#endif
