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
