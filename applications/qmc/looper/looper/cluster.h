/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_CLUSTER_H
#define LOOPER_CLUSTER_H

#include <boost/mpl/bool.hpp>

namespace looper {

struct cluster_info {
  cluster_info() : weight(0), sign(0) {}
  double weight;
  int sign;

  template<typename F, typename FIELD>
  struct field_accumulator {
    // typedef F cluster_fragment_t;
    template<typename D0, typename D1, typename D2>
    field_accumulator(D0 const&, D1 const&, D2 const&) {}
    void begin_s(int, double, int, int) const {}
    void end_s(int, double, int, int) const {}
  };
  template<typename F>
  struct field_accumulator<F, boost::mpl::true_> {
    typedef F cluster_fragment_t;
    field_accumulator(std::vector<cluster_info>& cl, std::vector<cluster_fragment_t> const& fr,
      std::vector<double> const& fd) : clusters(cl), fragments(fr), field(fd) {}
    void begin_s(int p, double t, int s, int c) {
      clusters[fragments[p].id()].weight -= - field[s] * (0.5-c) * t;
    }
    void end_s(int p, double t, int s, int c) {
      clusters[fragments[p].id()].weight += - field[s] * (0.5-c) * t;
    }
    std::vector<cluster_info>& clusters;
    std::vector<cluster_fragment_t> const& fragments;
    std::vector<double> const& field;
  };

  template<typename F, typename SIGN, typename IMPROVE>
  struct sign_accumulator {
    typedef F cluster_fragment_t;
    template<typename D0, typename D1, typename D2, typename D3>
    sign_accumulator(D0 const&, D1 const&, D2 const&, D3 const&) {}
    void bond_sign(int, int) const {}
    void site_sign(int, int) const {}
  };
  template<typename F>
  struct sign_accumulator<F, boost::mpl::true_, boost::mpl::true_> {
    typedef F cluster_fragment_t;
    sign_accumulator(std::vector<cluster_info>& cl, std::vector<cluster_fragment_t> const& fr,
      std::vector<int> const& bs, std::vector<int> const& ss) :
      clusters(cl), fragments(fr), bond_s(bs), site_s(ss) {}
    void bond_sign(int p, int b) { clusters[fragments[p].id()].sign += bond_s[b]; }
    void site_sign(int p, int s) { clusters[fragments[p].id()].sign += site_s[s]; }
    std::vector<cluster_info>& clusters;
    std::vector<cluster_fragment_t> const& fragments;
    std::vector<int> const& bond_s;
    std::vector<int> const& site_s;
  };

  template<typename F, typename FIELD, typename SIGN, typename IMPROVE>
  struct accumulator :
    private field_accumulator<F, FIELD>, private sign_accumulator<F, SIGN, IMPROVE> {
    typedef F cluster_fragment_t;
    typedef field_accumulator<F, FIELD> field_accumulator_t;
    typedef sign_accumulator<F, SIGN, IMPROVE> sign_accumulator_t;

    template<typename CL, typename FR, typename FD, typename BS, typename SS>
    accumulator(CL& cl, FR const& fr, FD const& fd, BS const& bs, SS const& ss) :
      field_accumulator_t(cl, fr, fd), sign_accumulator_t(cl, fr, bs, ss) {}
    void begin_s(int p, double t, int s, int c) {
      field_accumulator_t::begin_s(p, t, s, c);
      sign_accumulator_t::site_sign(p, s);
    }
    void begin_b(int p0, int p1, double t, int b, int s0, int s1, int c0, int c1) {
      field_accumulator_t::begin_s(p0, t, s0, c0);
      field_accumulator_t::begin_s(p1, t, s1, c1);
      sign_accumulator_t::bond_sign(p0, b);
    }
    void end_s(int p, double t, int s, int c) {
      field_accumulator_t::end_s(p, t, s, c);
      sign_accumulator_t::site_sign(p, s);
    }
    void end_b(int p0, int p1, double t, int b, int s0, int s1, int c0, int c1) {
      field_accumulator_t::end_s(p0, t, s0, c0);
      field_accumulator_t::end_s(p1, t, s1, c1);
      sign_accumulator_t::bond_sign(p0, b);
    }
    void start_bottom(int p, double t, int s, int c) { field_accumulator_t::begin_s(p, t, s, c); }
    void start(int p, double t, int s, int c) { field_accumulator_t::begin_s(p, t, s, c); }
    void stop(int p, double t, int s, int c) { field_accumulator_t::end_s(p, t, s, c); }
    void stop_top(int p, double t, int s, int c) { field_accumulator_t::end_s(p, t, s, c); }
  };
};

} // end namespace looper

#endif // LOOPER_CLUSTER_H
