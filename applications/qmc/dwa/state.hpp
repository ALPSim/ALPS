/*****************************************************************************
*
* ALPS Project Applications: Directed Worm Algorithm for Boson Hubbard Model 
*
* Copyright (C) 2011 by Lode Pollet      <pollet@itp.phys.ethz.ch>    (Inventor), 
*                       Ping Nang Ma     <pingnang@itp.phys.ethz.ch>  (Coder)   ,
*                       Matthias Troyer  <troyer@itp.phys.ethz.ch>    (Advisor) 
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

#ifndef ALPS_APPLICATIONS_STATE_HPP
#define ALPS_APPLICATIONS_STATE_HPP

#include <cassert>
#include <iostream>

namespace alps {
namespace applications {

template <class Number>
class State
{
public:
  typedef Number         integer_type;

  Number state() const  {  return _state;  }

  bool operator==(State const & rhs)  const  {  return _state == rhs._state;  }
  bool operator< (integer_type rhs)         const  {  return _state <  rhs;         }
  bool operator> (integer_type rhs)         const  {  return _state >  rhs;         }

  State (integer_type state_=0)   : _state(state_)          {}

  State & operator++()  {  ++_state;  return *this;  }
  State & operator--()  {  --_state;  return *this;  }

  template <class NNumber>
  friend std::ostream& operator<<(std::ostream & out, State<NNumber> const & state_)  {  out << state_._state;  return out;  }

  template <class NNumber>
  friend std::istream& operator>>(std::istream & in, State<NNumber> & state_)   {  in >> state_._state;  return in;  }

private:
	integer_type  _state;
};

}  // end namespace applications
}  // end namespace alps

#endif
