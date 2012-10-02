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

#ifndef ALPS_APPLICATIONS_WORLDLINE_HPP
#define ALPS_APPLICATIONS_WORLDLINE_HPP

#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>

#include <boost/lexical_cast.hpp>

#include <alps/osiris/comm.h>

#include "numeric.hpp"
#include "kink.hpp"


namespace alps {
namespace applications {

class Worldline
{
public:
  typedef unsigned int          Site;
  typedef double                Time;
  typedef unsigned short        State; 
  typedef kink<Site,Time,State> Kink;
  typedef std::vector<Kink>     Line;
  typedef std::vector<Line>     Lines;
  typedef Line::iterator         LineIterator;
  typedef Line::const_iterator   LineConstIterator;
  typedef Lines::iterator        LinesIterator;
  typedef Lines::const_iterator  LinesConstIterator;     

  Worldline() : _initialized(false)  {}
  Worldline (Site const total_number_of_sites_, Site const line_capacity_);

  bool  closed() const  { return _closed;  } 
  bool  open()   const  { return !_closed; }

  inline std::vector<State> states();
  inline std::vector<int>   vertices();

  Kink  wormhead() const  { return *_lineit; }

  bool  forward()     const  { return _lineit->forward();  }  
  bool  creation()    const  { return _lineit->creation(); }
  Site  site()        const  { return _lineit->site(); }
  Time  time()        const  { return _lineit->time(); }
  State state()       const  { return _lineit->state();       }
  State state_after() const  { return _lineit->state_after(); }

  LinesIterator  linesit()         {  return  _linesit;  }
  LineIterator   lineit()          {  return  _lineit;   }
  LineIterator   forwardlineit()   {  return  ++LineIterator(_lineit); }
  LineIterator   backwardlineit()  {  return  --LineIterator(_lineit); }

  LinesIterator  linesit (Site const site_)                          { return (_worldline.begin())+site_; }
  LineIterator   lineit  (LinesIterator linesit_, Time const time_)  { return std::lower_bound(linesit_->begin(), linesit_->end(), time_); }
  LineIterator   lineit  (Site const site_,       Time const time_)  { return lineit(linesit(site_), time_); } 

  State state (Site const site_, Time const time_)  { return (--lineit(site_, time_))->state_after(); }

  inline bool  unoccupied (Site const site_, Time const time_);
  void         wormpair_insertion (std::pair<Kink,Kink> wormpair_);
  void         wormpair_removal   ();

  bool  wormhead_touches_begin()  {  return (--LineIterator(_lineit) == _linesit->begin());  }
  bool  wormhead_touches_end()    {  return (++LineIterator(_lineit) == _linesit->end());    }

  void  wormhead_changes_its_direction                ()                  {  _lineit->turn();             }
  void  wormhead_moves_to_new_time_without_winding    (Time const time_)  {  _lineit->reset_time(time_);  }
  void  wormhead_moves_to_new_time_with_winding       (Time const time_);
  void  wormhead_jumps_to_new_site_and_inserts_vertex (LineIterator targetlineitprev_); 
  void  wormhead_deletes_vertex_and_jumps_to_new_site (LineIterator sourcelineit_);
  void  wormhead_relinks_vertex_and_jumps_to_new_site (LineIterator sourcelineit_, LineIterator targetlineitprev_);
  void  wormhead_crosses_vertex                       ();

  void  save (std::ostream & out) const;
  void  load (std::istream & in);

  void  save (alps::ODump & dump) const;
  void  load (alps::IDump & dump);

  friend std::ostream &  operator<< (std::ostream & out, Worldline const & object_);

private:
  bool  _initialized;
  bool  _closed;

  std::vector<Line>  _worldline;

  LinesIterator  _linesit;
  LineIterator   _lineit;
};

Worldline
  ::Worldline (Site const total_number_of_sites_, Site const line_capacity_)
    : _initialized       (true)
    , _closed            (true)
  {
    _worldline = Lines(total_number_of_sites_, Line(line_capacity_));
    for (Site _site=0; _site < total_number_of_sites_; ++_site)  {
      _worldline[_site].clear();
      _worldline[_site].push_back(Kink(_site));
    }
  }

inline std::vector<Worldline::State>
  Worldline
    ::states() 
    {
      std::vector<State> _states;
      _states.reserve(_worldline.size());
      for (LinesIterator it=_worldline.begin(); it!=_worldline.end(); ++it)
        _states.push_back(it->begin()->state());
      return _states;
    }

inline std::vector<int>
  Worldline
    ::vertices()
    {
      std::vector<int> _vertices;
      _vertices.reserve(_worldline.size());
      for (LinesIterator it=_worldline.begin(); it!=_worldline.end(); ++it)
        _vertices.push_back(it->size() -1);
      return _vertices;
    }

inline bool
  Worldline
    ::unoccupied (Site const site_, Time const time_)
    {
      if (time_ == 0.)
        return false;
      _linesit = linesit(site_);
      _lineit  = lineit(_linesit, time_);
      if (_lineit == _linesit->end())   
        return true;
      if (_lineit->time() == time_)
        return false;
      return true;
    }

void
  Worldline  
    ::wormpair_insertion (std::pair<Kink, Kink> wormpair_)  
    {
      std::vector<Kink> _wormpair = 
        wormpair_.first.forward() ? numeric::make_vector(wormpair_.second, wormpair_.first) : numeric::make_vector(wormpair_.first, wormpair_.second);

      _linesit = linesit(wormpair_.first.site());
      _linesit->insert(lineit(_linesit,wormpair_.first.time()), _wormpair.begin(), _wormpair.end());
      _lineit  = std::find_if(_linesit->begin(), _linesit->end(), std::mem_fun_ref(&Kink::wormhead));  

      if (forward())
        _lineit->time_increment(1e-8);
      else
        _lineit->time_decrement(1e-8);

      _closed = false;
    }

void
  Worldline
    ::wormpair_removal()
    {
      _lineit->forward() ? _linesit->erase(  LineIterator(_lineit), LineIterator(_lineit)+2) 
                         : _linesit->erase(--LineIterator(_lineit), LineIterator(_lineit)+1);  
      _closed = true;
    }

void
  Worldline
    ::wormhead_moves_to_new_time_with_winding (Time const time_)
    {
      if (forward())
      {
        std::rotate(++(_linesit->begin()), _lineit, _linesit->end());
        _lineit = ++(_linesit->begin());
      }
      else
      {
        std::rotate(++(_linesit->begin()), ++LineIterator(_lineit), _linesit->end());
        _lineit = --(_linesit->end());
      }
      _linesit->begin()->reset_state((++(_linesit->begin()))->state());
      _lineit->reset_time(time_);
    }

void
  Worldline
    ::wormhead_jumps_to_new_site_and_inserts_vertex (LineIterator targetlineitprev_)
    {
      _linesit = linesit(targetlineitprev_->site());
      _lineit->wormhead2vertex(targetlineitprev_->site());

      if (forward())
      {
        Kink  _vertexto = _lineit->generate_vertexpartner(targetlineitprev_->state_after());
        Kink  _wormhead = _vertexto.generate_wormhead(true);
        std::vector<Kink> _kinkpair = numeric::make_vector (_vertexto, _wormhead);
        _linesit->insert(++LineIterator(targetlineitprev_), _kinkpair.begin(), _kinkpair.end());
      }
      else
      {
        State _newstate = LineIterator(targetlineitprev_)->state_after();
        Kink  _vertexto = _lineit->generate_vertexpartner(creation() ? ++_newstate : --_newstate);
        Kink  _wormhead = _vertexto.generate_wormhead(false);
        std::vector<Kink> _kinkpair = numeric::make_vector (_wormhead, _vertexto);
        _linesit->insert(++LineIterator(targetlineitprev_), _kinkpair.begin(), _kinkpair.end());
      }
      _lineit = std::find_if(_linesit->begin(), _linesit->end(), std::mem_fun_ref(&Kink::wormhead));
    }

void
  Worldline
    ::wormhead_deletes_vertex_and_jumps_to_new_site(LineIterator sourcelineit_)
    {
      sourcelineit_->vertex2wormhead(forward());
      forward() ? _linesit->erase(  LineIterator(_lineit), ++(++LineIterator(_lineit)))
                : _linesit->erase(--LineIterator(_lineit),    ++LineIterator(_lineit) );
      _linesit = linesit(sourcelineit_->site());
      _lineit  = sourcelineit_;
    }

void
  Worldline
    ::wormhead_relinks_vertex_and_jumps_to_new_site (LineIterator sourcelineit_, LineIterator targetlineitprev_)
    {
      wormhead_deletes_vertex_and_jumps_to_new_site(sourcelineit_);
      wormhead_changes_its_direction();
      wormhead_jumps_to_new_site_and_inserts_vertex(targetlineitprev_);
    } 

void
  Worldline
    ::wormhead_crosses_vertex ()
    {
      if (forward())
      {
        _lineit->reset_time((++LineIterator(_lineit))->time());
        std::rotate(LineIterator(_lineit), ++LineIterator(_lineit), ++(++LineIterator(_lineit)));
        swap_state(*LineIterator(_lineit), *(++LineIterator(_lineit)));
        ++_lineit;
      }
      else
      {
        _lineit->reset_time((--LineIterator(_lineit))->time());
        std::rotate(--LineIterator(_lineit), LineIterator(_lineit), ++LineIterator(_lineit));
        swap_state(*LineIterator(_lineit), *(--LineIterator(_lineit)));
        --_lineit;
      }
    }

void
  Worldline
    ::save (std::ostream & out) const
  {
    for (LinesConstIterator _linesit = _worldline.begin(); _linesit != _worldline.end(); ++_linesit)
    for (LineConstIterator  _lineit  = _linesit->begin();  _lineit  != _linesit->end();  ++_lineit)
      _lineit->save(out);
  }


void
  Worldline
    ::load (std::istream & in)
  {
    std::for_each(_worldline.begin(), _worldline.end(), std::mem_fun_ref(&Line::clear));
    while (in.good())
    {
      Kink _kink(in);
      LinesIterator _linesit = linesit(_kink.site());  
      LineIterator  _lineit  = lineit(_linesit, _kink.time());
      _linesit->insert(_lineit, _kink);
    }
  }

/*
void
  Worldline
    ::save (alps::ODump & dump) const
  {
    for (LinesConstIterator _linesit=_worldline.begin(); _linesit!=_worldline.end(); ++_linesit)
    {
      dump << _linesit->size();
      for (LineConstIterator _lineit=_linesit->begin(); _lineit!=_linesit->end(); ++_lineit)
        _lineit->save(dump);
    }
  }

void
  Worldline
    ::load (alps::IDump & dump)
  {
    for (LinesIterator _linesit = _worldline.begin(); _linesit != _worldline.end(); ++_linesit)
    {
      _linesit->clear();
      int _line_size;  
      dump >> _line_size;

      for (int line_index_=0; line_index_ < _line_size; ++line_index_)
      {
        Kink _kink(dump);
        LineIterator _lineit  = lineit(_linesit, _kink.time());
        _linesit->insert(_lineit, _kink);
      }
    }
  }
*/

std::ostream &
  operator<< (std::ostream & out, Worldline const & object_)
  {
    if (!object_._initialized)   return out;

    out << "\n" << (object_._closed ? "Closed" : "Open") << " worldline has " << object_._worldline.size() << " number of sites."
        << "\n\n"
        << "\nWorldline Details"
        << "\n================="
        << "\n\n";
    for (Worldline::Site _site=0; _site < object_._worldline.size(); ++_site)  {
      out << "\nSite Nr : " << _site << "\t( Size: " << object_._worldline[_site].size() << " , Capacity: " << object_._worldline[_site].capacity() << " )"
               << "\n-----------------------------------------------------";
      std::copy((object_._worldline[_site]).begin(), (object_._worldline[_site]).end(), std::ostream_iterator<Worldline::Kink>(out,""));
    }
    return out;
  }



} // ending namespace applications
} // ending namespace alps

#endif
