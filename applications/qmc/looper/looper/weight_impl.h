/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_WEIGHT_IMPL_H
#define LOOPER_WEIGHT_IMPL_H

#include "weight.h"
#include "crop.h"
#include "model_parameter.h"

#include <alps/numeric/is_equal.hpp>
#include <alps/numeric/is_zero.hpp>
#include <alps/parameter.h>
#include <boost/array.hpp>
#include <algorithm> // for std::min std::max
#include <cmath>     // for std::abs

namespace looper {

//
// default graph weights for path-integral and SSE loop algorithms
//

struct site_weight_helper {

  BOOST_STATIC_CONSTANT(int, num_graphs = 1);
  int sign;
  double offset;
  boost::array<double, num_graphs> v;

  // loop equations:
  //
  //   - offset + v0 = 0
  //   - offset + v0 = 0
  //              v0 = |Hx|/2

  // standard solution:
  //
  //   v0 = |Hx|/2
  //   offset = |Hx|/2

  site_weight_helper() : sign(1), offset(0) { v.assign(0); }
  site_weight_helper(const site_parameter& p, alps::Parameters const& params) { init(p, params); }

  void init(const site_parameter& p, alps::Parameters const& /* params */) {
    sign = (p.hx >= 0 ? 1 : -1);
    v[0] = std::abs(p.hx) / 2;
    offset = v[0];
  }

  double weight() const { return v[0]; }
  bool has_weight() const { return alps::numeric::is_positive<1>(weight()); }

  void check(const site_parameter& p) const {
    if (!alps::numeric::is_zero<1>(-offset+v[0]) ||
        !alps::numeric::is_zero<1>(-offset+v[0]) ||
        !alps::numeric::is_equal<1>(v[0], sign * p.hx/2))
      boost::throw_exception(std::logic_error("site_parameter::check 1"));
    site_parameter pp(p.s, p.c, 2 * v[0] * sign, 0, 0);
    if (pp != p) {
      std::cerr << p.c << ' ' << p.hx << ' ' << pp.hz << std::endl;
      std::cerr << pp.c << ' ' << pp.hx << ' ' << pp.hz << std::endl;
      boost::throw_exception(std::logic_error("site_parameter::check 2"));
    }
  }
};

struct xxz_bond_weight_helper {

  BOOST_STATIC_CONSTANT(int, num_graphs = 4);
  int sign;
  double offset;
  boost::array<double, num_graphs> v;

  // loop equations:
  //
  //   - offset + v1 + v3 = - Jz/4
  //   - offset + v0 + v2 = + Jz/4
  //              v0 + v1 =   |Jxy|/2

  // standard solution:
  //
  // i) Jz <= -|Jxy|  (ferro-Ising)
  //      v1 = |Jxy|/2
  //      v3 = -(|Jxy| + Jz)/2
  //      v0 = v2 = 0
  //      offset = -Jz/4
  // ii) -|Jxy| <= Jz <= |Jxy|  (XY)
  //      v0 = (|Jxy| + Jz)/4
  //      v1 = (|Jxy| - Jz)/4
  //      v2 = v3 = 0
  //      offset = |Jxy|/4
  // iii) Jz >= |Jxy|  (antiferro-Ising)
  //      v0 = |Jxy|/2
  //      v2 = -(|Jxy| - Jz)/2
  //      v1 = v3 = 0
  //      offset = Jz/4

  // "ergodic" solutions (with additional parameter 0 < a < 1)
  //
  // i) Jz <= -|Jxy|  (ferro-Ising)
  //      same as the standard solution
  // ii) -|Jxy| <= Jz <= |Jxy|  (XY)
  //   ii-1) |Jxy| - Jz >= 2 a |Jxy|
  //      same as the standard solution
  //   ii-2) |Jxy| - Jz <= 2 a |Jxy|
  //      v0 = (1-a) |Jxy| / 2
  //      v1 = a |Jxy| / 2
  //      v2 = -((1-2a) |Jxy| - Jz)/2
  //      v3 = 0
  //      offset = Jz/4 + a |Jxy| / 2
  // iii) Jz >= |Jxy|  (antiferro-Ising)
  //      same as ii-2)

  xxz_bond_weight_helper() : sign(1), offset(0) { v.assign(0); }
  xxz_bond_weight_helper(const site_parameter& p, alps::Parameters const& params) {
    init(p, params);
  }
  xxz_bond_weight_helper(const bond_parameter_xxz& p, alps::Parameters const& params) {
    init(p, params);
  }
  xxz_bond_weight_helper(const bond_parameter_xyz& p, alps::Parameters const& params) {
    init(p, params);
  }

  void init(const bond_parameter_xxz& p, alps::Parameters const& params) {
    sign = (p.jxy <= 0 ? 1 : -1);
    double jxy = std::abs(p.jxy);
    double jz = p.jz;
    double a = crop_01(params.value_or_default("FORCE_SCATTER", 0.));
    if (alps::numeric::is_nonzero<1>(jxy + std::abs(jz))) {
      if (jxy - jz > 2 * a * jxy) {
        // standard solutions
        v[0] = crop_0(std::min(jxy/2, (jxy + jz)/4));
        v[1] = crop_0(std::min(jxy/2, (jxy - jz)/4));
        v[2] = crop_0(-(jxy - jz)/2.0);
        v[3] = crop_0(-(jxy + jz)/2);
      } else {
        // "ergodic" solutions
        v[0] = (1-a)*jxy/2;
        v[1] = a*jxy/2;
        v[2] = -((1-2*a)*jxy-jz)/2;
        v[3] = 0;
      }
    } else {
      v[0] = v[1] = v[2] = v[3] = 0;
    }
    offset = weight()/2;
  }
  void init(const bond_parameter_xyz& p, alps::Parameters const& params) {
    if (!alps::numeric::is_equal<>(p.jx, p.jy))
      boost::throw_exception(std::runtime_error("not an XXZ model"));
    init(bond_parameter_xxz(p.c, p.jx, p.jz), params);
  }
  void init(const site_parameter& p, alps::Parameters const& params) {
    bond_parameter_xyz bp(0, 0, 0, 2 * p.d);
    init(bp, params);
  }

  double weight() const { return v[0] + v[1] + v[2] + v[3]; }
  bool has_weight() const { return alps::numeric::is_positive<1>(weight()); }

  void check(const bond_parameter_xxz& p) const {
    if (!alps::numeric::is_equal<1>(-offset+v[0]+v[2],  p.jz/4) ||
        !alps::numeric::is_equal<1>(-offset+v[1]+v[3], -p.jz/4) ||
        !alps::numeric::is_equal<1>(v[0]+v[1], -sign * p.jxy/2))
      boost::throw_exception(std::logic_error("bond_parameter_xxz::check 1"));
    bond_parameter_xxz pp(p.c, -2 * (v[0] + v[1]) * sign, 2 * (v[0] - v[1] + v[2] - v[3]));
    if (pp != p) {
      std::cerr << p.c << ' ' << p.jxy << ' ' << p.jz << std::endl;
      std::cerr << pp.c << ' ' << pp.jxy << ' ' << pp.jz << std::endl;
      boost::throw_exception(std::logic_error("bond_parameter_xxz::check 2"));
    }
  }
};

struct xyz_bond_weight_helper {

  BOOST_STATIC_CONSTANT(int, num_graphs = 6);
  int sign;
  double offset;
  boost::array<double, num_graphs> v;

  // loop equations:
  //
  //   - offset + v1 + v3 + v5 = - Jz/4
  //   - offset + v0 + v2 + v4 = + Jz/4
  //              v0 + v1      = |Jx+Jy|/4 = Jp/2
  //              v4 + v5      = |Jx-Jy|/4 = Jm/2

  // standard solution with parameter scattering ratio: r (0 <= r <= 1)
  //
  // i) Jz <= -(Jp + (1-2r)Jm)  (ferro-Ising)
  //      v1 = Jp/2
  //      v3 = -(Jp + Jz + (1-2r)Jm)/2
  //      v4 = r Jm/2
  //      v5 = (1-r) Jm/2
  //      v0 = v2 = 0
  //      offset = (2r Jm-Jz)/4
  // ii) -(Jp + (1-2r)Jm) <= Jz <= Jp - (1-2r)Jm  (XY)
  //      v0 = (Jp + Jz + (1-2r)Jm)/4
  //      v1 = (Jp - Jz - (1-2r)Jm)/4
  //      v4 = r Jm/2
  //      v5 = (1-r) Jm/2
  //      v2 = v3 = 0
  //      offset = (Jp+Jm)/4
  // iii) Jz >= Jp - (1-2r)Jm  (antiferro-Ising)
  //      v0 = Jp/2
  //      v2 = -(Jp - Jz - (1-2r)Jm)/2
  //      v4 = r Jm/2
  //      v5 = (1-r) Jm/2
  //      v1 = v3 = 0
  //      offset = (2(1-r)Jm+Jz)/4

  // "ergodic" solutions (with additional parameter 0 < a < 1)
  //
  // i) Jz <= -(Jp + (1-2r)Jm) (ferro-Ising)
  //      same as the standard solution
  // ii) -(Jp + (1-2r)Jm) <= Jz <= Jp - (1-2r)Jm (XY)
  //   ii-1) Jp - Jz - (1-2r)Jm >= 2 a Jp
  //   ii-2) Jp - Jz - (1-2r)Jm <= 2 a Jp
  //      v0 = (1-a) Jp / 2
  //      v1 = a Jp / 2
  //      v2 = -((1-2a) Jp - Jz - (1-2r)Jm)/2
  //      v4 = r Jm/2
  //      v5 = (1-r) Jm/2
  //      v3 = 0
  //      offset = (2a Jp + Jz + 2(1-r)Jm)/4
  // iii) Jz >= Jp - (1-2r)Jm  (antiferro-Ising)
  //      same as ii-2)

  xyz_bond_weight_helper() : sign(1), offset(0) { v.assign(0); }
  xyz_bond_weight_helper(const site_parameter& p, alps::Parameters const& params) {
    init(p, params);
  }
  xyz_bond_weight_helper(const bond_parameter_xxz& p, alps::Parameters const& params) {
    init(p, params);
  }
  xyz_bond_weight_helper(const bond_parameter_xyz& p, alps::Parameters const& params) {
    init(p, params);
  }

  void init(const bond_parameter_xyz& p, alps::Parameters const& params) {
    if ((p.jx + p.jy) * (p.jx - p.jy) < 0)
      boost::throw_exception(std::runtime_error("(Jx+Jy) and (Jx-Jy) may not have opposite signs"));
    sign = ((p.jx + p.jy) <= 0 ? 1 : -1);
    double jp = std::abs((p.jx + p.jy)) / 2;
    double jm = std::abs((p.jx - p.jy)) / 2;
    double jz = p.jz;
    double a = crop_01(params.value_or_default("FORCE_SCATTER", 0.));
    double r = crop_01(params.value_or_default("SCATTERING_RATIO", 0.5));
    if (alps::numeric::is_nonzero<1>(jp + std::abs(jz))) {
      if (jp - jz - (1-2*r) * jm > 2 * a * jp) {
        // standard solutions
        v[0] = crop_0(std::min(jp/2, (jp + jz + (1-2*r) * jm)/4));
        v[1] = crop_0(std::min(jp/2, (jp - jz - (1-2*r) * jm)/4));
        v[2] = crop_0(-(jp - jz - (1-2*r) * jm)/2);
        v[3] = crop_0(-(jp + jz + (1-2*r) * jm)/2);
      } else {
        // "ergodic" solutions
        v[0] = (1 - a) * jp / 2;
        v[1] = a * jp / 2;
        v[2] = -((1 - 2*a) * jp - jz - (1 - 2*r) * jm)/2;
        v[3] = 0;
      }
    } else {
      v[0] = v[1] = v[2] = v[3] = 0;
    }
    if (alps::numeric::is_nonzero<1>(jm)) {
      v[4] = r * jm/2;
      v[5] = (1-r) * jm/2;
    } else {
      v[4] = v[5] = 0;
    }
    offset = weight()/2;
  }
  void init(const bond_parameter_xxz& p, alps::Parameters const& params) {
    init(bond_parameter_xyz(p.c, p.jxy, p.jxy, p.jz), params);
  }
  void init(const site_parameter& p, alps::Parameters const& params) {
    bond_parameter_xyz bp(0, 0, 0, 2 * p.d);
    init(bp, params);
  }

  double weight() const { return v[0] + v[1] + v[2] + v[3] + v[4] + v[5]; }
  bool has_weight() const { return alps::numeric::is_positive<1>(weight()); }

  void check(const bond_parameter_xyz& p) const {
    if (!alps::numeric::is_equal<1>(-offset+v[0]+v[2]+v[4],  p.jz/4) ||
        !alps::numeric::is_equal<1>(-offset+v[1]+v[3]+v[5], -p.jz/4) ||
        !alps::numeric::is_equal<1>(v[0]+v[1], -sign * (p.jx+p.jy)/4) ||
        !alps::numeric::is_equal<1>(v[4]+v[5], -sign * (p.jx-p.jy)/4)) {
      std::cerr << -offset+v[0]+v[2]+v[4] << '\t' << p.jz/4 << '\t'
                << -offset+v[1]+v[3]+v[5] << '\t' << -p.jz/4 << '\t'
                << v[0]+v[1] << '\t' << -sign * (p.jx+p.jy)/4 << '\t'
                << v[4]+v[5] << '\t' << -sign * (p.jx-p.jy)/4 << std::endl;
      boost::throw_exception(std::logic_error("bond_parameter_xyz::check 1"));
    }
    bond_parameter_xyz pp(p.c, -2 * (v[0] + v[1] + v[4] + v[5]) * sign,
      -2 * (v[0] + v[1] - v[4] - v[5]) * sign, 2 * (v[0] - v[1] + v[2] - v[3] + v[4] - v[5]));
    if (pp != p) {
      std::cerr << p.c << ' ' << p.jx << ' ' << p.jy << ' ' << p.jz << std::endl;
      std::cerr << pp.c << ' ' << pp.jx << ' ' << pp.jy << ' ' << pp.jz << std::endl;
      boost::throw_exception(std::logic_error("bond_parameter_xyz::check 2"));
    }
  }
};

template<typename SITE_WEIGHT_HELPER, typename BOND_WEIGHT_HELPER>
class weight_helper {
public:
  typedef SITE_WEIGHT_HELPER site_t;
  typedef BOND_WEIGHT_HELPER bond_t;

  typedef typename std::pair<int, site_t> site_weight_t;
  typedef typename std::pair<int, bond_t> bond_weight_t;
  typedef typename std::vector<std::pair<int, site_t> >::const_iterator site_weight_iterator;
  typedef typename std::vector<std::pair<int, bond_t> >::const_iterator bond_weight_iterator;

  BOOST_STATIC_CONSTANT(int, num_site_graphs = site_t::num_graphs);
  BOOST_STATIC_CONSTANT(int, num_bond_graphs = bond_t::num_graphs);

  template<class M, class LAT>
  weight_helper(M const& m, const LAT& lat, alps::Parameters const& params) {
    init(m, lat, params);
  }

  template<class M, class LAT>
  void init(M const& m, const LAT& lat, alps::Parameters const& params) {
    weight_ = 0;
    site_weights_.clear();
    bond_weights_.clear();
    BOOST_FOREACH(typename real_site_descriptor<LAT>::type rs, sites(lat.rg())) {
      site_t sw(m.site(rs, lat.rg()), params);
      BOOST_FOREACH(typename virtual_site_descriptor<LAT>::type vs, sites(lat, rs)) {
        site_weights_.push_back(std::make_pair(get(site_index_t(), lat.vg(), vs), sw));
        weight_ += sw.weight();
      }
    }
    BOOST_FOREACH(typename real_bond_descriptor<LAT>::type rb, bonds(lat.rg())) {
      bond_t bw(m.bond(rb, lat.rg()), params);
      BOOST_FOREACH(typename virtual_bond_descriptor<LAT>::type vb, bonds(lat, rb)) {
        bond_weights_.push_back(std::make_pair(get(bond_index_t(), lat.vg(), vb), bw));
        weight_ += bw.weight();
      }
    }
    if (m.has_d_term()) {
      BOOST_FOREACH(typename real_site_descriptor<LAT>::type rs, sites(lat.rg())) {
        bond_t bw(m.site(rs, lat.rg()), params);
        BOOST_FOREACH(typename virtual_bond_descriptor<LAT>::type vb, bonds(lat, rs)) {
          bond_weights_.push_back(std::make_pair(get(bond_index_t(), lat.vg(), vb), bw));
          weight_ += bw.weight();
        }
      }
    }

    site_offset_ = 0;
    bond_offset_ = 0;
    BOOST_FOREACH(site_weight_t const& w, site_weights_) site_offset_ += w.second.offset;
    BOOST_FOREACH(bond_weight_t const& w, bond_weights_) bond_offset_ += w.second.offset;
    if (m.has_d_term()) {
      BOOST_FOREACH(typename real_site_descriptor<LAT>::type rs, sites(lat.rg()))
        site_offset_ -= 0.25 * m.site(rs, lat.rg()).s.get_twice() * m.site(rs, lat.rg()).d;
    }
  }

  std::pair<site_weight_iterator, site_weight_iterator>
  site_weights() const { return std::make_pair(site_weights_.begin(), site_weights_.end()); }
  std::pair<bond_weight_iterator, bond_weight_iterator>
  bond_weights() const { return std::make_pair(bond_weights_.begin(), bond_weights_.end()); }

  double weight() const { return weight_; }
  double rho() const { return weight(); }
  double site_energy_offset() const { return site_offset_; }
  double bond_energy_offset() const { return bond_offset_; }

private:
  double weight_;
  double site_offset_;
  double bond_offset_;
  std::vector<site_weight_t> site_weights_;
  std::vector<bond_weight_t> bond_weights_;
};

} // namespace looper

#endif // LOOPER_WEIGHT_IMPL_H
