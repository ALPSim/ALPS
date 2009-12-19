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

#ifndef ALPS_APPLICATIONS_WKINK_H___
#define ALPS_APPLICATIONS_WKINK_H___

#define ICC_CYCLIC_ITERATOR

#include "cyclic_iterator.h"
#include "time_struct.h"
#include <boost/graph/graph_traits.hpp>
//#include <boost/pool/pool_alloc.hpp>

template <class Container>
inline cyclic_iterator<Container> adjacent(cyclic_iterator<Container> first, time_struct t1)
{
  Container& c = first.container();
  if (c.empty())
    return first;
#if defined ( USE_SET )
  typedef typename Container::iterator iterator;
  iterator last = c.begin()+c.size()-1;
  first = c.upper_bound(t1);
  if (first==c.end() || first==c.begin())
    first=last;
  else
    --first;
#elif defined( USE_VECTOR )
  typedef typename Container::iterator iterator;
  iterator last = c.begin()+c.size()-1;
  first = std::upper_bound(c.begin(),c.end(),t1);
  if (first==c.end() || first==c.begin())
    first=last;
  else
    --first;
#else
  if (first.valid()) {
    while(first->time() <= t1) {
      double t=first->time();
      ++first;
      if(t>=first->time()) break; // reached end
    }
    --first; // went one too far
  }
#endif
  return first;
}   // adjacent

//- kinks ---------------------------------------------------------------

template<class state_type>
class Kink {
public:
  typedef Kink<state_type> kink_type;

  Kink() {}
  Kink(double t, state_type s) : t_(t), state_(s) {}

  state_type state() const { return state_; }
  void set_state(state_type s) { state_=s; }

  time_struct time() const { return t_; }
  void set_time(time_struct t) { t_=t; }

  template <class C>
  cyclic_iterator<C> adjacent (cyclic_iterator<C> first) const { return ::adjacent(first,time());}
  
  void set_id(unsigned int id) { id_ =id;}
  bool is_linked(Kink const& other) { return id() == other.id();}
  unsigned int id() const { return id_;}

private:
  time_struct t_;    // time of kink
  state_type state_; // "number of particles" after kink
  unsigned int id_;
};   // Kink

template<class state_type>
inline bool operator<(const Kink<state_type>& x, const Kink<state_type>& y) {
  return x.time() < y.time();
}

template<class state_type>
inline bool operator==(const Kink<state_type>& x, const Kink<state_type>& y) {
  return x.time() == y.time();
}

template<class state_type>
inline bool operator==(const Kink<state_type>& x, double t) {
  return x.time() == t;
}

template<class state_type>
inline bool operator==(double t, const Kink<state_type>& x) {
  return x.time() == t;
}

template<class state_type>
inline bool operator!=(const Kink<state_type>& x, double t) {
  return x.time() != t;
}

template<class state_type>
inline bool operator!=(double t, const Kink<state_type>& x) {
  return x.time() != t;
}

template<class state_type>
inline alps::ODump& operator <<(alps::ODump& dump, const Kink<state_type>& k) {
  return dump << k.time() << k.state() << k.id();
}

template<class state_type> 
inline alps::IDump& operator >>(alps::IDump& dump, Kink<state_type>& k) {
  double t;
  state_type s;
  unsigned id;
  dump >> t >> s >> id;
  k.set_time(t);
  k.set_state(s);
  k.set_id(id);
  return dump;
}

template <class state_type>
inline std::ostream& operator<<(std::ostream& os, const Kink<state_type>& k)
{
  return os << "[ " << k.time() << " : " << long(k.state()) << " ]";
}

//- worm head ---------------------------------------------------------------

template <class C, class G>
class Wormhead
{
public:
  typedef G graph_type;
  typedef C container_type;
  typedef ::cyclic_iterator<container_type> cyclic_iterator;
  typedef typename cyclic_iterator::base_iterator iterator;
  
  typedef typename boost::graph_traits<graph_type>::vertex_descriptor vertex_descriptor;

  Wormhead(const graph_type& g, std::vector<container_type>& c) : graph(g), kinks_(c), valid_(false) {};
    
  vertex_descriptor site() const { return s_;}
  time_struct time() const { return time_;}
   
  void set_site(vertex_descriptor s) { 
    s_=s; 
    iterator_=cyclic_iterator(kinks_[s]); 
    valid_=false;
  }
  
  void set_time(time_struct t) { time_=t;}
  
  cyclic_iterator kink() { if (!valid_) revalidate(); return iterator_;}

  const Wormhead& operator=(iterator it) { 
    iterator_=it; 
    set_time(it->time());
    valid_=true; 
    return *this;
  }
  
  void invalidate(int s)
  {  
    if (s==site())
      valid_=false;
  }
  
private:
  void revalidate() {
    container_type& c=kinks_[site()];
#if defined ( USE_SET )
    iterator k= c.find(time_);
    if (k==c.end()) {
#elif defined ( USE_VECTOR )
    iterator k= std::lower_bound(c.begin(),c.end(),time_);
    if (k->time()!=time_) {
#else
    iterator k= std::find(c.begin(),c.end(),time_);
    if (k==c.end()) {
#endif
      for (k=c.begin(); k!=c.end();++k)
        std::cerr << "Have at: " << k->time() << "\n";
      std::cerr << site() << " " << time() << "\n";
      boost::throw_exception(std::logic_error("Reached end!"));
    }
    assert(k!=c.end());
    iterator_=k;
    valid_=true;
  }

  const graph_type& graph;
  std::vector<container_type>& kinks_;
  time_struct time_;
  vertex_descriptor s_;
  cyclic_iterator iterator_;
  bool valid_;
};

#endif

