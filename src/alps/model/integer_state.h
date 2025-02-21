/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2004 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
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

#ifndef ALPS_MODEL_INTEGER_STATE_H
#define ALPS_MODEL_INTEGER_STATE_H

#include <boost/integer/static_log2.hpp>

namespace alps {

template <class I, int N=1> class integer_state;

template <class I, int N>
class integer_state {
public:
  BOOST_STATIC_CONSTANT(int, bits = boost::static_log2<N>::value+1);
  BOOST_STATIC_CONSTANT(int, mask = (1<<bits)-1);
  typedef I representation_type;
  
  class reference {
  public:
    reference(I& s, int i) : state_(s), shift_(i*bits) {}
    operator int() const { return (state_ >> shift_) & mask;}
    template <class T>
    reference& operator=(T x)
    {
      state_ &= ~(mask<<shift_);
      state_ |= ((mask & x)<<shift_);
      return *this;
    }
  private:
    I& state_;
    std::size_t shift_;
  };
  
  integer_state(representation_type x=0) : state_(x) {}
  
  template <class J>
  integer_state(const std::vector<J>& x) : state_(0)
  { 
    for (int i=0;i<x.size();++i)  
      state_ |=(x[i]<<(i*bits));
  }
  int operator[](int i) const { return (state_>>i)&mask;}
  reference operator[](int i) { return reference(state_,i);}
  operator representation_type() const { return state_;}
  representation_type state() const { return state_;}
private:
  representation_type state_;
};

template <class I>
class integer_state<I,1> {
public:
  typedef I representation_type;
  
  class reference {
  public:
    reference(I& s, int i) : state_(s), mask_(1<<i) {}
    operator int() const { return (state_&mask_ ? 1 : 0);}
    template <class T>
    reference& operator=(T x)
    {
      if (x)
        state_|=mask_;
      else
        state_&=~mask_;
      return *this;
    }
  private:
    I& state_;
    I mask_;
  };
  
  integer_state(representation_type x=0) : state_(x) {}
  
  template <class J>
  integer_state(const std::vector<J>& x) : state_(0)
  { 
    for (int i=0;i<x.size();++i)  
      if(x[i])
        state_ |=(1<<i);
  }
  int operator[](int i) const { return (state_>>i)&1;}
  reference operator[](int i) { return reference(state_,i);}
  operator representation_type() const { return state_;}
  representation_type state() const { return state_;}
private:
  representation_type state_;
};

template <class I, int N>
bool operator == (integer_state<I,N> x, integer_state<I,N> y)
{ return x.state() == y.state(); }

template <class I, int N>
bool operator < (integer_state<I,N> x, integer_state<I,N> y)
{ return x.state() < y.state(); }

} // namespace alps

#endif
