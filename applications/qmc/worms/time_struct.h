/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2004 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>
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

#ifndef TIME_STRUCT_H___
#define TIME_STRUCT_H___

#include <alps/osiris.h>
#include <cmath>

class time_struct {
public:
  time_struct(double t=0.) : t_(checkfull(t)) {}
  const time_struct& operator=(double t) { t_=checkfull(t); return *this;}
  operator double() const { return t_;}
  time_struct operator-(double t) const { return time_struct(checkone(t_-t),0);}
  time_struct operator+(double t) const { return time_struct(checkone(t_+t),0);}
  double operator-(time_struct t) const { return t>=t_ ? t_-t+1. : t_-t;}

private:
  time_struct(double t,int) : t_(t) {}
  double checkfull(double t) const { return t>0. ? std::fmod(t,1.) : std::fmod(t,1.)+1.;}
  double checkone(double t) const {return t<=0. ? t+1. : (t>=1. ? t-1. : t);}
  double t_;
};

inline alps::ODump& operator << (alps::ODump& dump, const time_struct& t) {
  return dump << double(t);
}

inline alps::IDump& operator >> (alps::IDump& dump, time_struct& t) {
  double td;
  dump >> td;
  t=td;
  return dump;
}

#endif

