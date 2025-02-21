/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2007 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_OPERATOR_H
#define LOOPER_OPERATOR_H

#include "time.h"
#include "type.h"

#include <boost/call_traits.hpp>
#include <boost/throw_exception.hpp>
#include <stdexcept>

// workaround for SuSE 11.4, which defines macro TIME in pyconfig.h
#ifdef TIME
# undef TIME
#endif

namespace looper {

struct local_operator_type {
  BOOST_STATIC_CONSTANT(int, diagonal    = 0 /* 00 */);
  BOOST_STATIC_CONSTANT(int, offdiagonal = 1 /* 01 */);
  BOOST_STATIC_CONSTANT(int, identity    = 2 /* 10 */);
};


template<class QMC, class LOC_G, class TIME = imaginary_time<> >
class local_operator;

template<class LOC_G, class TIME>
class local_operator<sse, LOC_G, TIME> {
public:
  typedef LOC_G                              local_graph_t;
  typedef typename local_graph_t::location_t location_t;

  local_operator() : type_(local_operator_type::identity), loc_() {}
  local_operator(int type, const location_t& loc) : type_(type), loc_(loc) {}
  local_operator(const local_graph_t& g)
    : type_(local_operator_type::diagonal | (g.type() << 2)), loc_(g.loc()) {}

  void flip() { type_ ^= 1; }
  int type() const { return type_ & 3; }
  bool is_diagonal() const { return type() == local_operator_type::diagonal; }
  bool is_offdiagonal() const { return type() == local_operator_type::offdiagonal; }
  bool is_identity() const { return type() == local_operator_type::identity; }

  void clear_graph() { type_ &= 3; }
  local_graph_t graph() const { return local_graph_t(graph_type(), loc_); }
  void assign_graph(const local_graph_t& g) {
    assert(is_identity() || loc_ == g.loc());
    // if type==identity then type will be set to diagonal, otherwize unchanged
    type_ = (type_ & 1) | (g.type() << 2);
    loc_ = g.loc();
  }
  int graph_type() const { return type_ >> 2; }
  bool is_frozen_site_graph() const {
    return local_graph_t::site_graph_t::is_frozen(graph_type());
  }
  bool is_frozen_bond_graph() const {
    return local_graph_t::bond_graph_t::is_frozen(graph_type());
  }


  int loop_0() const { return loop0; }
  int loop_1() const { return loop1; }
  int loop_l() const
  { return local_graph_t::loop_l(graph_type(), loop0, loop1); }
  int loop_u() const
  { return local_graph_t::loop_u(graph_type(), loop0, loop1); }
  int loop_l0() const
  { return local_graph_t::loop_l0(graph_type(), loop0, loop1); }
  int loop_l1() const
  { return local_graph_t::loop_l1(graph_type(), loop0, loop1); }
  int loop_u0() const
  { return local_graph_t::loop_u0(graph_type(), loop0, loop1); }
  int loop_u1() const
  { return local_graph_t::loop_u1(graph_type(), loop0, loop1); }

  const location_t& loc() const { return loc_; }
  int pos() const { return loc_.pos(); }
  int source() const { return loc_.source(); }
  int target() const { return loc_.target(); }
  bool is_site() const { return loc_.is_site(); }
  bool is_bond() const { return loc_.is_bond(); }

  void save(alps::ODump& dump) const { dump << type_ << loc_; }
  void load(alps::IDump& dump) { dump >> type_ >> loc_; }

private:
  int type_;
  location_t loc_;

public:
  int loop0, loop1; // no checkpoint
};

template<class LOC_G, class TIME>
class local_operator<path_integral, LOC_G, TIME> : public local_operator<sse, LOC_G, TIME> {
public:
  typedef local_operator<sse, LOC_G, TIME>   super_type;
  typedef typename super_type::local_graph_t local_graph_t;
  typedef typename super_type::location_t    location_t;
  typedef TIME                               time_t;

  local_operator() : super_type(), time_(0) {}
  local_operator(int type, const location_t& loc); // not defined
  local_operator(int type, const location_t& loc, double time)
    : super_type(type, loc), time_(time) {}
  local_operator(const local_graph_t& g); // not defined
  local_operator(const local_graph_t& g, double time)
    : super_type(g), time_(time) {}

  bool is_identity() const; // not defined
  typename boost::call_traits<time_t>::param_type time() const { return time_; }

  void save(alps::ODump& dump) const { super_type::save(dump); dump << time_; }
  void load(alps::IDump& dump) { super_type::load(dump); dump >> time_; }

private:
  time_t time_;
};

template<class QMC, class LOC_G, class TIME>
int type(const local_operator<QMC, LOC_G, TIME>& op) { return op.type(); }
template<class QMC, class LOC_G, class TIME>
bool is_diagonal(const local_operator<QMC, LOC_G, TIME>& op) { return op.is_diagonal(); }
template<class QMC, class LOC_G, class TIME>
bool is_offdiagonal(const local_operator<QMC, LOC_G, TIME>& op) { return op.is_offdiagonal(); }
template<class LOC_G, class TIME>
bool is_identity(const local_operator<sse, LOC_G, TIME>& op) { return op.is_identity(); }

template<class QMC, class LOC_G, class TIME>
int pos(const local_operator<QMC, LOC_G, TIME>& op) { return op.pos(); }
template<class QMC, class LOC_G, class TIME>
int source(const local_operator<QMC, LOC_G, TIME>& op) { return op.source(); }
template<class QMC, class LOC_G, class TIME>
int target(const local_operator<QMC, LOC_G, TIME>& op) { return op.target(); }
template<class QMC, class LOC_G, class TIME>
bool is_site(const local_operator<QMC, LOC_G, TIME>& op) { return op.is_site(); }
template<class QMC, class LOC_G, class TIME>
bool is_bond(const local_operator<QMC, LOC_G, TIME>& op) { return op.is_bond(); }

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

template<class QMC, class LOC_G, class TIME>
alps::ODump& operator<<(alps::ODump& dump,
                        const looper::local_operator<QMC, LOC_G, TIME>& op)
{ op.save(dump); return dump; }

template<class QMC, class LOC_G, class TIME>
alps::IDump& operator>>(alps::IDump& dump,
                        looper::local_operator<QMC, LOC_G, TIME>& op)
{ op.load(dump); return dump; }

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

#endif // LOOPER_TYPE_H
