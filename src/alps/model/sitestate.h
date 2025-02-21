/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2004 by Matthias Troyer <troyer@comp-phys.org>
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

#ifndef ALPS_MODEL_SITESTATE_H
#define ALPS_MODEL_SITESTATE_H

#include <alps/model/sitebasisdescriptor.h>
#include <vector>
#include <iostream>

namespace alps {

template <class I>
class site_state : public std::vector<half_integer<I> > {
public:
  typedef half_integer<I> quantumnumber_type;
  typedef typename std::vector<half_integer<I> >::const_iterator const_iterator;
  site_state() {}
  site_state(const std::vector<half_integer<I> >& x) : std::vector<half_integer<I> >(x)  {}
};


template <class I>
class single_qn_site_state {
public:
  typedef half_integer<I> representation_type;
  typedef half_integer<I> quantumnumber_type;
  typedef std::size_t size_type;
  
  single_qn_site_state() {}
  single_qn_site_state(representation_type x) : state_(x)  {}
  template <class J>
  single_qn_site_state(const std::vector<half_integer<J> >& x) { assert(x.size()==1); state_=x[0];}
  operator representation_type() const { return state_;}
  representation_type state() const { return state_;}
  representation_type& state() { return state_;}
private:
  representation_type state_;
};


template <class I>
bool operator < (const single_qn_site_state<I>& x,  const single_qn_site_state<I>& y)
{
  return x.state() < y.state();
}

template <class I>
bool operator > (const single_qn_site_state<I>& x,  const single_qn_site_state<I>& y)
{
  return x.state() > y.state();
}

template <class I>
bool operator == (const single_qn_site_state<I>& x,  const single_qn_site_state<I>& y)
{
  return x.state() == y.state();
}

template <class I>
bool operator <= (const single_qn_site_state<I>& x,  const single_qn_site_state<I>& y)
{
  return x.state() <= y.state();
}

template <class I>
bool operator >= (const single_qn_site_state<I>& x,  const single_qn_site_state<I>& y)
{
  return x.state() >= y.state();
}


template <class I>
half_integer<I> get_quantumnumber(const site_state<I>& s, typename site_state<I>::size_type i)
{
  if (i>=s.size())
    boost::throw_exception(std::logic_error("Called get_quantumnumber with illegal index"));
  return s[i];
}

template <class I>
half_integer<I> get_quantumnumber(const single_qn_site_state<I>& s, std::size_t i)
{
  if (i!=0)
    boost::throw_exception(std::logic_error("Called get_quantumnumber with illegal index"));
  return s.state();
}

template <class I>
half_integer<I>& get_quantumnumber(site_state<I>& s, typename site_state<I>::size_type i)
{
  if (i>=s.size())
    boost::throw_exception(std::logic_error("Called get_quantumnumber with illegal index"));
  return s[i];
}

template <class I>
half_integer<I>& get_quantumnumber(single_qn_site_state<I>& s, std::size_t i)
{
  if (i!=0)
    boost::throw_exception(std::logic_error("Called get_quantumnumber with illegal index"));
  return s.state();
}

template <class I>
std::size_t get_quantumnumber_index(const std::string& n, const SiteBasisDescriptor<I>& b)
{
  for (std::size_t i=0;i<b.size();++i) {
    if (b[i].name()==n)
      return i;
  }
  return b.size();
}

template <class S, class I>
typename S::quantumnumber_type get_quantumnumber(const S& s, const std::string& n, const SiteBasisDescriptor<I>& b)
{
  return get_quantumnumber(s,get_quantumnumber_index(n,b));
}

template <class I, class S>
bool is_fermionic(const SiteBasisDescriptor<I>& b, const S& s)
{
  bool f=false;
  for (std::size_t i=0;i<b.size();++i)
    if (b[i].fermionic() && is_odd(get_quantumnumber(s,i)))
      f=!f;
  return f;
}

}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace alps {
#endif

template <class I>
std::ostream& operator<<(std::ostream& out, const alps::site_state<I>& s)
{
  out << "|";
  for (typename alps::site_state<I>::const_iterator it=s.begin();it!=s.end();++it)
    out << *it << " ";
  out << ">";
  return out;        
}

template <class I>
std::ostream& operator<<(std::ostream& out, const alps::single_qn_site_state<I>& s)
{
  out << "|" << s.state() << ">";
  return out;        
}


#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // namespace alps
#endif

#endif
