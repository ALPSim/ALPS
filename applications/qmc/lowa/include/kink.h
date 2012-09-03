/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Applications
*
* Copyright (C) 2006-2010 by Lode Pollet <lpollet@physics.harvard.edu>,
*                            Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>
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

/* $Id: kink.h 3520 2010-03-03 16:55:00Z tamama $ */


/*
 *
 * 1) Code modification      -- done
 * 2) Replacing raw pointers -- not done yet
 *
 */

#ifndef ALPS_APPLICATION_LOWA_KINK_H
#define ALPS_APPLICATION_LOWA_KINK_H

#include <iostream>
#include <list>
#include <algorithm>

#include "./typedef.h"


namespace alps {
namespace applications {
namespace lowa {


template<class S, class T>
class kink
{
public:
  // typedef
  typedef          S                               site_type;
  typedef          T                               time_type;   
  typedef          fock_basis_type                 basis_type;     
  typedef          kink_assoc_size_type            size_type;
  typedef          kink_assoc_index_type           index_type;
  typedef typename std::list<kink<S,T> >::iterator iterator_type;


  // set functions
  inline void set_name  (int        const name)   { _name   = name;   }
  inline void set_time  (time_type  const time)   { _time   = time;   }
  inline void set_before(basis_type const before) { _before = before; }
  inline void set_after (basis_type const after)  { _after  = after;  } 
  inline void set_from  (site_type  const from)   { _from   = from;   }
  inline void set_to    (site_type  const to)     { _to     = to;     }
  inline void set_assoc (index_type const index, iterator_type const elem)
  {  _assoc[index] = elem;  }


  // get functions
  inline int        name()   const  { return _name;  }
  inline time_type  time()   const  { return _time;  }
  inline basis_type before() const  { return _before;}
  inline basis_type after()  const  { return _after; }
  inline site_type  from()   const  { return _from;  }
  inline site_type  to()     const  { return _to;    }
  inline size_type  size()   const  { return _size;  }
  inline iterator_type assoc(index_type index) const 
  { return _assoc[index];}


private:
  // member objects
  bool _is_initialized;

  int            _name;    // numerical name of an element; used to save/load associations of elements to/from file
  time_type      _time;    // time of the interaction
  basis_type     _before;  // occupation on the site before the interaction
  basis_type     _after;   // occupation on the site after the interaction
  site_type      _from;    // site from which the particle hops
  site_type      _to;      // site to which the particle hops
  size_type      _size;    // number of associations  

  iterator_type* _assoc;   // associations; i.e. iterator to the kink on nb sites that are equal or just greater in time

public:
  // destructor
  ~kink()  
  {  if (_size != 0) { delete [] _assoc; }  }
    

  // constructors
  kink() {}

  kink(size_type const size)
    : _size(size) 
    , _assoc(new iterator_type [size])
  {}

  kink(size_type const size, basis_type const before, basis_type const after, site_type const from, site_type const to, time_type const time) 
    : _time(time)
    , _before(before)
    , _after(after)
    , _from(from)
    , _to(to)
    , _size(size)
    , _assoc(new iterator_type [size])
  {}

  kink(kink<site_type,time_type> const & mykink)
    : _name(-1)          // *** check its significance...
    , _time(mykink._time)
    , _before(mykink._before)
    , _after(mykink._after)
    , _from(mykink._from)
    , _to(mykink._to)
    , _size(mykink._size)
    , _assoc(new iterator_type [mykink._size])
  {
    std::copy(&mykink._assoc[0],&mykink._assoc[_size],&_assoc[0]);
  }


  void init(size_type size)  {  _size = size;  _assoc = new iterator_type [size]; }


  // assignment operator
  kink& operator=(kink<site_type,time_type> const & mykink)
  {
    if (this == &mykink)   {  return *this;  }

    _time   = mykink._time;
    _before = mykink._before;
    _after  = mykink._after;
    _from   = mykink._from;
    _to     = mykink._to;
    _size   = mykink._size;

    std::copy(&mykink._assoc[0],&mykink._assoc[_size],&_assoc[0]);

    return *this;
  }


  // i/o operator 
  friend std::ostream& operator<<(std::ostream& out, kink<site_type,time_type> const & mykink) {
    out << mykink._time   << "\t" 
        << mykink._from   << "\t"
        << mykink._to     << "\t"
        << mykink._before << "\t" 
        << mykink._after  << "\t"
        << mykink._size   << "\t";
    return (out);
  }

  friend std::istream& operator>>(std::istream& in, kink<site_type,time_type>& mykink) {
    in >> mykink._time 
       >> mykink._from 
       >> mykink._to 
       >> mykink._before 
       >> mykink._after
       >> mykink._size;
    return in;
  }


  // comparison operator
  friend bool operator==(kink<site_type,time_type> const & lhs, kink<site_type,time_type> const & rhs) {
    bool is_equal;
    ((lhs._from == rhs._from) && (lhs.to() == rhs.to()) && (lhs.time() == rhs.time()) && (lhs.before() == rhs.before()) && (lhs.after() == rhs.after()) ) ? is_equal=true : is_equal=false;
    return is_equal;
  }    


  // print    *** I don't care about this, Tama.
  void print() {
    std::cout << "\n" << this << "\t" << _time << "\t" << _from << "\t" << _to << "\t" << _before << "\t" << _after;
    for (site_type j = 0; j < _size; j++) std::cout << "\t" << &(*_assoc[j]) << "\t" << _assoc[j]->time();
  }

};


} // ending namespace lowa
} // ending applications
} // ending alps

#endif
