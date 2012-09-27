/*****************************************************************************
*
* ALPS Project Applications: Directed Worm Algorithm 
*
* Copyright (C) 2012 by Lode Pollet      <pollet@phys.ethz.ch>, 
*                       Ping Nang Ma     <pingnang@phys.ethz.ch>,
*                       Matthias Troyer  <troyer@phys.ethz.ch>  
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

#ifndef ALPS_APPLICATIONS_KINK_HPP
#define ALPS_APPLICATIONS_KINK_HPP

#include <cassert>
#include <iostream>

#include <alps/osiris/comm.h>

namespace alps {
namespace applications {

inline bool increasing (bool const forward_, bool const creation_)  { return ((creation_ && !forward_) || (!creation_ && forward_)); }

template <class Site, class Time, class State>
class kink 
{
public:
  typedef kink<Site,Time,State> Kink;

  bool  initialized () const  {  return _initialized;  }
  bool  permanent   () const  {  return _permanent; }
  bool  identity    () const  {  return _identity;  }
  bool  nonidentity () const  {  return !_identity; }
  bool  worm        () const  {  return _worm;      }
  bool  wormhead    () const  {  return _worm && !_permanent;  }
  bool  wormtail    () const  {  return _worm && _permanent;   }
  bool  vertex      () const  {  return !_worm;     }
  bool  creation    () const  {  return _creation;  }
  bool  annihilation() const  {  return !_creation; } 
  bool  forward     () const  {  return _forward;   }
  bool  backward    () const  {  return !_forward;  }
  bool  increasing  () const  {  return alps::applications::increasing(_forward, _creation);  }
  bool  decreasing  () const  {  return !alps::applications::increasing(_forward, _creation); }
  State state       () const  {  return _state;  }
  State state_after () const  {  State _newstate(_state); return (_identity ? _newstate : (_creation ? ++_newstate : --_newstate));  }
  Time  time        () const  {  return _time;   }
  Site  site        () const  {  return _site;   }
  Site  partnersite () const  {  return _partnersite;  }

  Kink generate_vertexpartner(State partnerstate_) const  { return Kink(false, !_creation, false, partnerstate_, _time, _partnersite, _site, true); }
  Kink generate_wormhead     (bool  forward_)      const  { return Kink(true, !_creation, forward_, state_after(), _time, _site, _site, false); }
    
  bool operator< (Time const & time_) const  {  return (time() <  time_);  }
  bool operator> (Time const & time_) const  {  return (time() >  time_);  }
  bool operator<=(Time const & time_) const  {  return (time() <= time_);  }
  bool operator>=(Time const & time_) const  {  return (time() >= time_);  }

  kink()                  { _initialized = false; }
  kink(std::istream & in) { in  >> _initialized >> _permanent >> _identity >> _worm >> _creation >> _forward >> _state >> _time >> _site >> _partnersite; }
  kink(Site site_);
  kink(bool worm_, bool creation_, bool forward_, State state_, Time time_, Site  site_, Site partnersite_, bool permanent_);

  void  vertex2wormhead    (bool forward_)                    {  _worm = true;   _permanent = false;  _forward = forward_;  }
  void  wormhead2vertex    (Site partnersite_)                {  _worm = false;  _permanent = true;   _partnersite = partnersite_;  }

  void  reset_permanent    (bool  permanent_)    {  _permanent   = permanent_;    }  
  void  reset_creation     (bool  creation_)     {  _creation    = creation_;     }
  void  reset_forward      (bool  forward_)      {  _forward     = forward_;      }
  void  reset_state        (State state_)        {  _state       = state_;        }
  void  reset_time         (Time  time_)         {  _time        = time_;         }
  void  reset_site         (Site  site_)         {  _site        = site_;         }
  void  reset_partnersite  (Site  partnersite_)  {  _partnersite = partnersite_;  }

  void  time_increment     (Time  deltatime_)    {  _time += deltatime_;  }
  void  time_decrement     (Time  deltatime_)    {  _time -= deltatime_;  }

  Kink&  turn()         {  _forward = !_forward;  return *this;  } 
  Kink&  operator++()   {  ++_state;              return *this;  }
  Kink&  operator--()   {  --_state;              return *this;  }

  template <class SSite, class TTime, class SState>
  friend inline void swap_state (kink<SSite,TTime,SState> & object1_, kink<SSite,TTime,SState> & object2_)  {  std::swap(object1_._state, object2_._state);  }

  void save (std::ostream & out) const  {  out << " " << _initialized << " " << _permanent << " " << _identity << " " << _worm << " " << _creation << " " << _forward << " " << _state << " " << _time << " " << _site << " " << _partnersite;  }

  template <class SSite, class TTime, class SState>
  friend std::ostream &  operator<< (std::ostream & out, kink<SSite,TTime,SState> const & object_);

private:
  bool  _initialized;
  bool  _permanent;      // are the iterators (including its partner) permanently set ? 
  bool  _identity;       // designed for diagonal measurement convenience (not physical) ; not physical / physical ?
  bool  _worm;           // worm / vertex ?
  bool  _creation;       // creation / annihilation ?
  bool  _forward;        // forward / backward ?
  State _state;          // state at time-epsilon
  Time  _time;           // time periodic in [0 , 1) 
  Site  _site;           // site
  Site  _partnersite;    // site of partner
};

template <class Site, class Time, class State>
kink<Site,Time,State>
  ::kink (Site site_)
  : _initialized   (true)
  , _permanent     (true)
  , _identity      (true)
  , _worm          (false)
  , _creation      (false)
  , _forward       (false)
  , _state         (State())
  , _time          (0.)
  , _site          (site_)
  , _partnersite   (site_)
  {}
  
template <class Site, class Time, class State>
kink<Site,Time,State>
  ::kink
    ( bool  worm_
    , bool  creation_
    , bool  forward_
    , State state_  
    , Time  time_
    , Site  site_
    , Site  partnersite_
    , bool  permanent_
    )
  : _initialized   (true)
  , _permanent     (permanent_)
  , _identity      (false)
  , _worm          (worm_)
  , _creation      (creation_)
  , _forward       (forward_)
  , _state         (state_)
  , _time          (time_)
  , _site          (site_)
  , _partnersite   (partnersite_)
  {}

template <class Site, class Time, class State>
std::ostream &  operator<< (std::ostream & out, kink<Site,Time,State> const & object_)  
{
  if (!object_._initialized)   return out; 
  if (object_._identity)
  {
    out << "\nIdentity"
        << "\nState = " << object_.state() 
        << "\nTime : " << object_._time
        << "\nSite : " << object_._site
        << "\n";
    return out;
  }
  else if (!object_._worm)
  {
    assert(object_._permanent);
    out << "\nVertex"
        << (object_._creation ? "\t-- creation" : "\t--annihilation")
        << "\nState = " << object_._state << " -> " << object_.state_after()
        << "\nTime : " << object_._time
        << "\nSite : " << object_._site << "\t( Partner Site : " << object_._partnersite << " )"
        << "\n";
    return out;
  }
  else if (object_._permanent)
  {
    out << "\nWormtail" 
        << (object_._creation ? "\t-- creation" : "\t--annihilation")
        << "\nState = " << object_._state << " -> " << object_.state_after()
        << "\nTime : " << object_._time
        << "\nSite : " << object_._site
        << "\n";
    return out;
  }
  else 
  {
    out << "\nWormhead" 
        << (object_._creation ? "\t-- creation" : "\t--annihilation")
        << (object_._forward  ? "\t-- forward"  : "\t--backward")
        << "\nState = " << object_._state << " -> " << object_.state_after()
        << "\nTime : " << object_._time
        << "\nSite : " << object_._site
        << "\n";
    return out;
  }
}

} // end namespace applications
} // end namespace alps

#endif
