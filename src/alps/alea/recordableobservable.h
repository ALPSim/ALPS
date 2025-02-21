/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2004 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Beat Ammon <ammon@ginnan.issp.u-tokyo.ac.jp>,
*                            Andreas Laeuchli <laeuchli@itp.phys.ethz.ch>,
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

#ifndef ALPS_ALEA_RECORDABLEOBSERVABLE_H
#define ALPS_ALEA_RECORDABLEOBSERVABLE_H

#include <alps/config.h>
#include <boost/throw_exception.hpp>
#include <stdexcept>

namespace alps {

//=======================================================================
// RecordableObservable
//
// an observable that can store new measurements
//-----------------------------------------------------------------------

template <class T=double, class SIGN=double>
class RecordableObservable
{
public:
  typedef T value_type;
  typedef SIGN sign_type;

  /// just a default constructor  
  RecordableObservable() {}
  virtual ~RecordableObservable() {}

  /// add another measurement to the observable
  virtual void operator<<(const value_type& x) =0;
  /// add another measurement to the observable
  virtual void add(const value_type& x) { operator<<(x);}
  /// add an explcitly signed measurement to the observable
  virtual void add(const value_type& x, sign_type s) { 
    if (s==1)
      add(x);
    else
      boost::throw_exception(std::logic_error("Called add of unsigned dobservable with a sign that is not 1"));
  }
 
   };
}

#endif // ALPS_ALEA_SIMPLEOBSERVABLE_H
