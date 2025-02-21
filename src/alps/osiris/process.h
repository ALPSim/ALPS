/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2005 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

#ifndef OSIRIS_PROCESS_H
#define OSIRIS_PROCESS_H

#include <alps/config.h>
#include <alps/osiris/std/vector.h>
#include <alps/osiris/dump.h>

#include <cstdlib>
#include <string>

namespace alps {

 
/** a process descriptor. */
    
class Process
{
public:
  
  // CONSTRUCTORS
  
  explicit Process(int); // constructor for process on unknown host
  Process() : tid(-1) {} // invalid Process
  
  // MEMBER FUNCTIONS
  
  void load(IDump&); // load from a dump
  void save(ODump&) const; // save into a dump

  bool valid() const; // is this a valid Process descriptor ?
  bool local() const; // is this the current Process?
  
  inline operator int () const {return tid;}
  
  bool operator==(const Process& p) const
  { return (tid==p.tid);}

  bool operator!=(const Process& p)  const
  { return (tid!=p.tid);}

  /// sorting criterion for two processes?
  bool operator<(const Process& p) const
  { return (tid<p.tid);}

private:
  int tid; // the unique Process id
};

typedef std::vector<Process> ProcessList;
}


#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace alps {
#endif

inline std::ostream& operator<<(std::ostream& out, const alps::Process& p)
{
  out << int(p);
  return out;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace alps
#endif

#endif // OSIRIS_PROCESS_H
