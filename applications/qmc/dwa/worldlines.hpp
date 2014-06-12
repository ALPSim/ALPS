/*****************************************************************************
*
* ALPS Project Applications: Directed Worm Algorithm 
*
* Copyright (C) 2013 by Matthias Troyer  <troyer@phys.ethz.ch> ,
*                       Lode Pollet      <pollet@phys.ethz.ch> ,
*                       Ping Nang Ma     <pingnang@phys.ethz.ch> 
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
*
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

#ifndef WORLDLINES_HPP
#define WORLDLINES_HPP

#include <iostream>
#include <utility>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <functional>
#include <numeric>
#include <limits>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <alps/hdf5.hpp>

// ==================================================
// kink class
// ==================================================

class kink
{
public:
  kink() {}
  kink(unsigned int siteindicator_, double time_ =0., unsigned short state_ =0) : _siteindicator(siteindicator_), _time(time_), _state(state_)  {}
  kink(std::istream & in) { in  >> _siteindicator >> _time >> _state; }

  std::string representation() const  { return "kink object: site indicator = " + boost::lexical_cast<std::string>(siteindicator()) + " , time = " + boost::lexical_cast<std::string>(time()) + " , state = " + boost::lexical_cast<std::string>(state()); } 

  unsigned int   siteindicator() const  {  return _siteindicator;   }
  double         time()          const  {  return _time;   }
  unsigned short state()         const  {  return _state;  }

  void  set_siteindicator (unsigned int  siteindicator_) { _siteindicator = siteindicator_; }
  void  set_time          (double  time_)          { _time  = time_;  }
  void  set_state         (unsigned short state_)         { _state = state_; }
  void  state_increment()  { ++_state; }
  void  state_decrement()  { --_state; }

  kink& operator+=(double deltatime_)  { _time += deltatime_; return *this; }
  kink& operator-=(double deltatime_)  { _time -= deltatime_; return *this; }
  kink& operator= (double time_)       { _time = time_;       return *this; }
  kink& operator++()                   { ++_state;            return *this; }
  kink& operator--()                   { --_state;            return *this; }
  bool  operator< (double const & time_) const  {  return (_time <  time_);  }
  bool  operator> (double const & time_) const  {  return (_time >  time_);  }
  bool  operator<=(double const & time_) const  {  return (_time <= time_);  }
  bool  operator>=(double const & time_) const  {  return (_time >= time_);  }

  friend void swap_state (kink & obj1_, kink & obj2_)  { std::swap(obj1_._state, obj2_._state); }
  friend std::ostream & operator<<(std::ostream & out, kink const & obj)  { out << "\t" << obj._siteindicator << "\t" << obj._time << "\t" << obj._state;  return out; }

  void save(alps::hdf5::archive & ar) const;
  void load(alps::hdf5::archive & ar);

private:
  unsigned int   _siteindicator;  // partner site for vertex ; site for wormhead, wormtail, dummy
  double         _time;           // time periodic in [0 , 1) 
  unsigned short _state;          // state at time+epsilon
};

void kink::save(alps::hdf5::archive & ar) const
{
  ar << alps::make_pvp("siteindicator", _siteindicator)
     << alps::make_pvp("time"         , _time)
     << alps::make_pvp("state"        , _state)
     ;
}

void kink::load(alps::hdf5::archive & ar)
{
  ar >> alps::make_pvp("siteindicator", _siteindicator)
     >> alps::make_pvp("time"         , _time)
     >> alps::make_pvp("state"        , _state)
     ;
} 

// ==================================================
// worldlines class
// ==================================================

class worldlines
{
public:
  typedef std::vector<kink> line;
  typedef std::vector<line> lines;
  typedef std::pair<lines::iterator, line::iterator>   location_type;

  worldlines() {}
  worldlines(unsigned int num_sites_)
  {
    _worldlines.resize(num_sites_);
    for (unsigned int site=0; site < num_sites_; ++site)
      _worldlines[site].push_back(kink(site));
  }

  worldlines open_worldlines (kink const & kink_) const;

  std::vector<std::vector<unsigned int> >   worldlines_siteindicator() const;
  std::vector<std::vector<double> >         worldlines_time()          const;
  std::vector<std::vector<unsigned short> > worldlines_state()         const;

  unsigned int  num_sites()                   const  { return _worldlines.size(); }
  unsigned int  num_kinks(unsigned int site_) const  { return _worldlines[site_].size(); }

  std::vector<unsigned short> states() const;

  unsigned int site_state(unsigned int site_) const  { return _worldlines[site_][0].state(); }

  int net_number_of_directed_hops(unsigned int site_, unsigned int partnersite_) const;

  location_type  location(unsigned int site_, double time_)  
  { return std::make_pair(_worldlines.begin()+site_, std::lower_bound(_worldlines[site_].begin(), _worldlines[site_].end(), time_)); }

  bool  location_is_kink_unoccupied(location_type const & location_, double time_) const
  { return (time_ == 0. ? false : location_.second == location_.first->end() ? true : location_.second->time() == time_ ? false : true); }

  unsigned short state_before (location_type const & location_) const  { return (location_.second-1)->state(); }
  unsigned short state        (location_type const & location_) const  { return (location_.second == location_.first->end() ? location_.first->begin()->state() : location_.second->state()); }  

  void output(std::ostream & out, unsigned int i);

  friend std::ostream &  operator<<(std::ostream & out, worldlines const & obj_);

  std::string representation() const  { std::ostringstream oss; oss << *this; return oss.str(); }

  void save_old1(alps::hdf5::archive & ar) const;
  void load_old1(alps::hdf5::archive & ar);

  void save(alps::hdf5::archive & ar) const;
  void load(alps::hdf5::archive & ar);

  void save(std::string const & filename) const;
  void load(std::string const & filename);

  bool is_valid(unsigned short Nmax);

private:
  lines _worldlines;
};


// ==================================================
// worldlines member functions
// ==================================================

worldlines worldlines::open_worldlines(kink const & kink_) const
{
  worldlines new_worldlines;
  new_worldlines._worldlines = _worldlines;
  location_type kink_insertion_location = new_worldlines.location(kink_.siteindicator(), kink_.time());
  kink_insertion_location.first->insert(kink_insertion_location.second, kink_);
  return new_worldlines;
}

std::vector<std::vector<unsigned int> > worldlines::worldlines_siteindicator() const
{
  std::vector<std::vector<unsigned int> > _vec;
  _vec.reserve(num_sites());
  for (unsigned int site_=0; site_<num_sites(); ++site_)
  {
    std::vector<unsigned int> _vec_inner;
    _vec_inner.reserve(num_kinks(site_));
    for (unsigned int vertex_=0; vertex_<num_kinks(site_); ++vertex_)
      _vec_inner.push_back(_worldlines[site_][vertex_].siteindicator());
    _vec.push_back(_vec_inner);
  }
  return _vec;
}

std::vector<std::vector<double> > worldlines::worldlines_time() const
{
  std::vector<std::vector<double> > _vec;
  _vec.reserve(num_sites());
  for (unsigned int site_=0; site_<num_sites(); ++site_)
  { 
    std::vector<double> _vec_inner;
    _vec_inner.reserve(num_kinks(site_));
    for (unsigned int vertex_=0; vertex_<num_kinks(site_); ++vertex_)
      _vec_inner.push_back(_worldlines[site_][vertex_].time());
    _vec.push_back(_vec_inner);
  }
  return _vec;
}

std::vector<std::vector<unsigned short> > worldlines::worldlines_state() const
{
  std::vector<std::vector<unsigned short> > _vec;
  _vec.reserve(num_sites());
  for (unsigned int site_=0; site_<num_sites(); ++site_)
  { 
    std::vector<unsigned short> _vec_inner;
    _vec_inner.reserve(num_kinks(site_));
    for (unsigned int vertex_=0; vertex_<num_kinks(site_); ++vertex_)
      _vec_inner.push_back(_worldlines[site_][vertex_].state());
    _vec.push_back(_vec_inner);
  }
  return _vec;
}

std::vector<unsigned short> worldlines::states() const
{
  std::vector<unsigned short> _states;
  _states.reserve(_worldlines.size());
  for (lines::const_iterator it=_worldlines.begin(); it!=_worldlines.end(); ++it)
    _states.push_back(it->begin()->state());
  return _states;
}

int worldlines::net_number_of_directed_hops(unsigned int site_, unsigned int partnersite_) const
{
  int net_number = 0;
  line::const_iterator it = _worldlines[site_].begin();
  ++it;
  for (; it != _worldlines[site_].end(); ++it)
    if (it->siteindicator() == partnersite_)
      (it->state() > (it-1)->state()) ? ++net_number : --net_number;
  return net_number;
}

void worldlines::output(std::ostream & out, unsigned int i) 
{
  out << "\nSite : " << i << "\n";
  std::copy(_worldlines[i].begin(), _worldlines[i].end(), std::ostream_iterator<kink>(out,"\n"));
}

std::ostream & operator<<(std::ostream & out, worldlines const & obj_)
  {
    out << "\n==================================================\n";
    out << "\nworldlines : num_sites() = " << obj_.num_sites()
        << "\n\n"
        << "\nWorldline Details"
        << "\n================="
        << "\n\n";
    for (unsigned int site=0; site < obj_.num_sites(); ++site)  {
      out << "site = " << site << "\t( num_kinks(site) = " << obj_.num_kinks(site) << " )\n"
          << "-----------------------------------------------------\n";
      std::copy((obj_._worldlines[site]).begin(), (obj_._worldlines[site]).end(), std::ostream_iterator<kink>(out,"\n"));
      out << "\n";
    }
    out << "\n==================================================\n";
    return out;
  }

void worldlines::save_old1(alps::hdf5::archive & ar) const
{
  ar << alps::make_pvp("/simulation/worldlines/num_sites"  , num_sites())
     << alps::make_pvp("/simulation/worldlines/worldlines" , _worldlines);
}

void worldlines::load_old1(alps::hdf5::archive & ar)
{
  unsigned int _archive_num_sites;
  ar >> alps::make_pvp("/simulation/worldlines/num_sites"  , _archive_num_sites); 

  if (num_sites() != _archive_num_sites)
    boost::throw_exception(std::runtime_error("Error in loading worldline object. Reason: wrong data structure."));

  ar >> alps::make_pvp("/simulation/worldlines/worldlines" , _worldlines);
}

void worldlines::save(alps::hdf5::archive & ar) const
{
  // number of sites
  ar << alps::make_pvp("/simulation/worldlines/num_sites", num_sites());

  // number of kinks
  std::vector<unsigned int> _local_num_kinks;
  _local_num_kinks.reserve(num_sites());
  for (unsigned int i=0; i<num_sites(); ++i) 
    _local_num_kinks.push_back(num_kinks(i));
  unsigned int _num_kinks = std::accumulate(_local_num_kinks.begin(), _local_num_kinks.end(), 0);

  ar << alps::make_pvp("/simulation/worldlines/num_kinks", _num_kinks);
  ar << alps::make_pvp("/simulation/worldlines/local_num_kinks", _local_num_kinks);

  // siteindicator, time, state
  std::vector<unsigned int>    _siteindicator;
  std::vector<double>          _time;
  std::vector<unsigned short>  _state;

  _siteindicator.reserve(_num_kinks);
  _time.reserve(_num_kinks);
  _state.reserve(_num_kinks);
  for (unsigned int i=0; i<num_sites(); ++i)
  for (unsigned int j=0; j<num_kinks(i); ++j) 
  {
    _siteindicator.push_back(_worldlines[i][j].siteindicator());
    _time.push_back(_worldlines[i][j].time());
    _state.push_back(_worldlines[i][j].state());
  }
  ar << alps::make_pvp("/simulation/worldlines/siteindicator", _siteindicator);
  ar << alps::make_pvp("/simulation/worldlines/time", _time);
  ar << alps::make_pvp("/simulation/worldlines/state", _state);
}

void worldlines::load(alps::hdf5::archive & ar)
{
  _worldlines.clear();

  unsigned int _archive_num_sites;
  ar >> alps::make_pvp("/simulation/worldlines/num_sites" , _archive_num_sites);

  _worldlines.resize(_archive_num_sites);

  std::vector<unsigned int> _local_num_kinks;
  ar >> alps::make_pvp("/simulation/worldlines/local_num_kinks" , _local_num_kinks);

  std::vector<unsigned int>    _siteindicator;
  std::vector<double>          _time;
  std::vector<unsigned short>  _state;
  ar >> alps::make_pvp("/simulation/worldlines/siteindicator", _siteindicator);
  ar >> alps::make_pvp("/simulation/worldlines/time", _time);
  ar >> alps::make_pvp("/simulation/worldlines/state", _state);

  for (unsigned int i=0, idx=0; i<num_sites(); ++i) {
    _worldlines[i].reserve(2*_local_num_kinks[i]);
    for (unsigned int j=0; j<_local_num_kinks[i]; ++j, ++idx) 
      _worldlines[i].push_back(kink(_siteindicator[idx],_time[idx],_state[idx]));
  }
}

void worldlines::save(std::string const & filename) const
{
  alps::hdf5::archive ar(filename.c_str(), "w");
  save(ar);
}

void worldlines::load(std::string const & filename) 
{
  alps::hdf5::archive ar(filename.c_str());
  load(ar);
}

bool worldlines::is_valid(unsigned short Nmax)  
{
  bool valid = true;

  // testing vertex state
  for (unsigned int i=0; i<_worldlines.size(); ++i)
  {
    if (_worldlines[i][0].state() > Nmax)
      valid = false;

    for (unsigned int j=1; j < _worldlines[i].size(); ++j)
    {
      if (_worldlines[i][j].state() > Nmax)
        valid = false;
      short this_state_increment = _worldlines[i][j].state() - _worldlines[i][j-1].state();
      if (!(this_state_increment == 1 || this_state_increment == -1))
        valid = false;
    }

    if (!valid) {
      std::cout << "\nError: testing vertex state fails...\n";
      std::cout << "site " << i << " : ";
      for (unsigned int j=0; j < _worldlines[i].size(); ++j)
        std::cout << _worldlines[i][j].state() << "  ";
      std::cout << "\n";
      return false;
    }
  } 

  // testing vertex time
  for (unsigned int i=0; i<_worldlines.size(); ++i)
  {
    if (_worldlines[i][0].time() != 0.)
      valid = false;
    
    for (unsigned int j=1; j < _worldlines[i].size(); ++j) {
      if (_worldlines[i][j].time() < 0. || _worldlines[i][j].time() > 1.)  // time must be bounded within (0.,1.)
        valid = false;
      if (_worldlines[i][j].time() <= _worldlines[i][j-1].time())   // time must be increasingly ordered...
        valid = false;
    }

    if (!valid) {
      std::cout << "\nError: testing vertex time fails...\n";
      std::cout << "site " << i << " : ";
      for (unsigned int j=0; j < _worldlines[i].size(); ++j)
        std::cout << _worldlines[i][j].time() << "  ";
      std::cout << "\n";
      return false;
    }
  } 

  // testing vertex siteindicator
  for (unsigned int i=0; i<_worldlines.size(); ++i)
  {
    for (unsigned int j=0; j < _worldlines[i].size(); ++j)
      if (_worldlines[i][j].siteindicator() >= _worldlines.size())
        valid = false;

    if (!valid)  {
      std::cout << "\nError: testing vertex siteindicator fails...\n";
      std::cout << "site " << i << " : ";
      for (unsigned int j=0; j < _worldlines[i].size(); ++j)
        std::cout << _worldlines[i][j].siteindicator() << "  ";
      std::cout << "\n";
      return false;
    }
  }

  // testing vertex pairing
  for (unsigned int i=0; i<_worldlines.size(); ++i)
  {
    if (_worldlines[i].size() > 1)
    {
      for (unsigned int j=1; j < _worldlines[i].size(); ++j) {
        kink linkedto(*(location(_worldlines[i][j].siteindicator(), _worldlines[i][j].time()).second)); 

        if (linkedto.time() != _worldlines[i][j].time())
          valid = false;

        if (linkedto.siteindicator() != i)
          valid = false;
      }
    }

    if (!valid) {
      std::cout << "\nError: testing vertex paring fails...\n";    
      std::cout << "site " << i << "\n";
      for (unsigned int j=1; j < _worldlines[i].size(); ++j) {
        std::cout << "kink : " << _worldlines[i][j] << " , linkedto : " << *(location(_worldlines[i][j].siteindicator(), _worldlines[i][j].time()).second) << "\n";
      }
      return false;
    }
  }

  return valid;
}


// ==================================================
// wormpair class
// ==================================================

class wormpair 
{
public:
  typedef worldlines::lines::iterator linesiterator;
  typedef worldlines::line::iterator  lineiterator;
  typedef worldlines::location_type location_type;

  wormpair() {}
  wormpair(location_type location_, kink const & kink_, bool forward_, bool creation_);
 
  unsigned int   site()           const { return _wormhead.siteindicator(); }
  double         time()           const { return _wormhead.time(); }
  double         reference_time() const { return _wormtail.time(); }
  unsigned short state_before()   const { return _creation ? _wormhead.state()-1 : _wormhead.state()+1; }
  unsigned short state()          const { return _wormhead.state(); }

  bool forward()    const  { return _forward;  }
  bool creation()   const  { return _creation; }
  bool increasing() const  { return wormpair::increasing(_forward, _creation); }

  bool neighbor2wormtail() const  { return _neighbor2wormtail; }

  unsigned int   next_partnersite() const { return _next->siteindicator(); }
  double         next_time()        const { return _next->time();  }
  unsigned short next_state()       const { return _next->state(); }

  unsigned int   wormtail_site()  const { return _wormtail.siteindicator(); }
  double         wormtail_time()  const { return _wormtail.time();  }
  unsigned short wormtail_state() const { return _wormtail.state(); }

  kink wormhead() const  { return _wormhead; }
  kink wormtail() const  { return _wormtail; }

  bool wormhead_is_same_type_as_next() const { return _forward ? ((_creation && (state()+1 == next_state())) || (!_creation && (state() == next_state()+1))) : ((_creation && (state_before() == (_next-1)->state()+1)) || (!_creation && (state_before()+1 == (_next-1)->state()))); }

  unsigned short wormpair_state() const { return _wormpair_state; }

  lineiterator forwardlineiterator()  const { return (wormhead_touches_end()   ? _location.first->begin()+1 : _location.second);   }
  lineiterator backwardlineiterator() const { return (wormhead_touches_begin() ? _location.first->end()-1   : _location.second-1); }
  lineiterator nextlineiterator()     const { return (_forward ? forwardlineiterator() : backwardlineiterator()); }

  bool wormhead_touches_begin()    const  { return _location.second == _location.first->begin()+1; }
  bool wormhead_touches_end()      const  { return _location.second == _location.first->end(); }
  bool wormhead_reaches_wormtail() const  { return _next->siteindicator() == site(); } 

  double time2wormtail() const  { return (_forward ? wormpair::modone(wormtail_time() - time()) : wormpair::modone(time() - wormtail_time())); }
  double time2next()     const  { return (_forward ? wormpair::modone(next_time() - time()) : wormpair::modone(time() - next_time())); }

  inline void wormhead_turns_around()  { _forward = !_forward; _next = nextlineiterator(); }
  inline void wormhead_moves_to_new_time(double time_, bool winding_over_time_=false);
  inline void wormhead_inserts_vertex_and_jumps_to_new_site(location_type const & targetlocation_);
  inline void wormhead_deletes_vertex_and_jumps_to_new_site(location_type const & sourcelocation_);
  inline void wormhead_relinks_vertex_and_jumps_to_new_site(location_type const & sourcelocation_, location_type const & targetlocation_);
  inline void wormhead_crosses_vertex();
  inline void wormhead_annihilates_wormtail()  { _location.first->erase(_next); }

  void set_neighbor2wormtail(bool neighbor2wormtail_)  { _neighbor2wormtail = neighbor2wormtail_; }

  static double modone (double value_)
  {
    while(value_<0. || value_>=1.)
      value_ += (value_>=1. ? -std::floor(value_) : (value_<0. ? -std::ceil(value_)+1 : 0.));
    return value_;
  }

  static bool increasing (bool forward_, bool creation_)  { return ((creation_ && !forward_) || (!creation_ && forward_)); }

  friend std::ostream &  operator<<(std::ostream & out, wormpair const & obj_);

  std::string representation() const  { std::ostringstream oss; oss << *this; return oss.str(); }

  bool is_valid(unsigned short Nmax);

private:
  unsigned short _wormpair_state;

  kink _wormtail;
  kink _wormhead;
  bool _forward;
  bool _creation;

  bool _neighbor2wormtail;

  location_type _location;
  lineiterator  _next;
};

// ==================================================
// wormpair member functions
// ==================================================

wormpair::wormpair(location_type location_, kink const & kink_, bool forward_, bool creation_)
  : _wormtail          (kink_)
  , _wormhead          (kink_)
  , _forward           (forward_)
  , _creation          (creation_)
  , _neighbor2wormtail (false)
  , _location          (location_)
{
  _wormhead += _forward ? std::numeric_limits<double>::epsilon() : -std::numeric_limits<double>::epsilon();

  unsigned short _newstate = (increasing() ? kink_.state()+1 : kink_.state()-1);
  _forward ? _wormtail.set_state(_newstate) : _wormhead.set_state(_newstate);

  _wormpair_state = (_creation ? state() : wormtail_state());

  _location.second = _location.first->insert(_location.second, _wormtail);
  if (_forward)
    ++_location.second;
  _next = nextlineiterator();
}

inline void wormpair::wormhead_moves_to_new_time(double time_, bool winding_over_time_)
{
  if (winding_over_time_) 
  {
    increasing() ? _location.first->begin()->state_increment() : _location.first->begin()->state_decrement();
    _location.second = _forward ? _location.first->begin()+1 : _location.first->end();
  }
  _wormhead.set_time(time_);
}

inline void wormpair::wormhead_inserts_vertex_and_jumps_to_new_site(location_type const & targetlocation_)
{
  // insert source vertex
  const kink _source(targetlocation_.first->begin()->siteindicator(), time(), state());
  _location.first->insert(_location.second, _source);

  // insert target vertex
  const unsigned short _targetstate    = (targetlocation_.second-1)->state();
  const unsigned short _targetnewstate = (increasing() ? _targetstate+1 : _targetstate-1);  

  if (_forward)
  {
    const kink _vertex(site(), time(), _targetnewstate);
    _wormhead = kink(_source.siteindicator(), time()+std::numeric_limits<double>::epsilon(), _targetstate);
    _location = targetlocation_;
    _location.second = _location.first->insert(_location.second, _vertex);
    ++_location.second;
  }
  else
  {
    const kink _vertex(site(), time(), _targetstate);
    _wormhead = kink(_source.siteindicator(), time()-std::numeric_limits<double>::epsilon(), _targetnewstate);
    _location = targetlocation_;
    _location.second = _location.first->insert(_location.second, _vertex);
  }

  // set next iterator
  _next = nextlineiterator();
}

inline void wormpair::wormhead_deletes_vertex_and_jumps_to_new_site(location_type const & sourcelocation_)
{
  // deletes target vertex
  _location.first->erase(_next);

  // deletes source vertex
  _wormhead = *sourcelocation_.second;
  _wormhead.set_siteindicator(sourcelocation_.first->begin()->siteindicator());
  _location = sourcelocation_;
  _location.second = _location.first->erase(_location.second);

  // set next iterator
  _next = nextlineiterator();
}

inline void wormpair::wormhead_relinks_vertex_and_jumps_to_new_site(location_type const & sourcelocation_, location_type const & targetlocation_)
{
  // deletes old target vertex
  _location.first->erase(_next);

  // resets partnersite on source and changes direction 
  sourcelocation_.second->set_siteindicator(targetlocation_.first->begin()->siteindicator());
  _forward = !_forward;

  // inserts new target vertex
  const unsigned short _targetstate    = (targetlocation_.second-1)->state();
  const unsigned short _targetnewstate = (increasing() ? _targetstate+1 : _targetstate-1);    

  if (_forward)
  {
    const kink _vertex(sourcelocation_.first->begin()->siteindicator(), sourcelocation_.second->time(), _targetnewstate);
    _wormhead = kink(sourcelocation_.second->siteindicator(), sourcelocation_.second->time()+std::numeric_limits<double>::epsilon(), _targetstate);
    _location = targetlocation_;
    _location.second = _location.first->insert(_location.second, _vertex);
    ++_location.second;
  }
  else
  {
    const kink _vertex(sourcelocation_.first->begin()->siteindicator(), sourcelocation_.second->time(), _targetstate);
    _wormhead = kink(sourcelocation_.second->siteindicator(), sourcelocation_.second->time()-std::numeric_limits<double>::epsilon(), _targetnewstate);
    _location = targetlocation_;
    _location.second = _location.first->insert(_location.second, _vertex);
  }

  // set next iterator
  _next = nextlineiterator();
}

inline void wormpair::wormhead_crosses_vertex()
{
  swap_state(_wormhead,*_next);
  if (_forward)
  {
    _wormhead.set_time(next_time()+std::numeric_limits<double>::epsilon());
    ++_location.second;
  }
  else
  {
    _wormhead.set_time(next_time()-std::numeric_limits<double>::epsilon());
    --_location.second;
  }
  _next = nextlineiterator();
}

std::ostream & operator<<(std::ostream & out, wormpair const & obj_)
{
  out << "\n----------------------------------------------------";
  out << "\nWormtail:\t" << obj_._wormtail << "\t" << (obj_._creation ? "annihilation" : "creation")  
      << "\n"
      << "\nVertex before : \t" << *(obj_.backwardlineiterator()) << "\t" << (obj_.backwardlineiterator()->siteindicator() == obj_._wormhead.siteindicator() ? " -- wormtail " : " -- vertex ")
      << "\nWormhead      : \t" << obj_._wormhead << "\t" << (obj_._forward?"forward":"backward") << "\t" << (obj_._creation ?"creation":"annihilation") 
      << "\nVertex after  : \t" << *(obj_.forwardlineiterator()) << "\t" << (obj_.forwardlineiterator()->siteindicator() == obj_._wormhead.siteindicator() ? " -- wormtail " : " -- vertex ")
      << "\n"
      << "\nKink before : \t" << *(obj_._location.second-1)
      << "\nWormhead    : \t" << obj_._wormhead << "\t" << (obj_._forward?"forward":"backward") << "\t" << (obj_._creation ?"creation":"annihilation")
      << "\nKink after  : \t" << (obj_.wormhead_touches_end() ? *(obj_._location.first->begin()) : *(obj_._location.second))
      << "\n"
      << "\nNext : \t" << *(obj_._next) << "\t" << (obj_.wormhead_reaches_wormtail() ? " -- wormtail " : " -- vertex ")
      << "\n----------------------------------------------------"
      << "\n";
  return out;
}

bool wormpair::is_valid(unsigned short Nmax)
{
  bool valid = true;

  // check state
  if (_wormhead.state() > Nmax)
    valid = false;
  if ((_location.second-1)->state() > Nmax)
    valid = false;
  if (wormhead_touches_end())
  {
    if (_location.first->begin()->state() > Nmax)
      valid = false;
  }
  else
  {
    if (_location.second->state() > Nmax)
      valid = false;
  }

  // check time
  if (_wormhead.time() <= 0. || _wormhead.time() >= 1.)
    valid = false;
  if (_wormhead.time() <= (_location.second-1)->time())
    valid = false;
  if (!wormhead_touches_end())
    if (_wormhead.time() >= _location.second->time())
      valid = false;

  return valid;
}

#endif
