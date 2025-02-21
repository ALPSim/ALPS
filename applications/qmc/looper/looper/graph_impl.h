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

#ifndef LOOPER_GRAPH_IMPL_H
#define LOOPER_GRAPH_IMPL_H

#include "graph.h"
#include "location_impl.h"
#include "random_choice.h"

#include <alps/numeric/is_nonzero.hpp>
#include <boost/array.hpp>
#include <boost/throw_exception.hpp>
#include <boost/tuple/tuple.hpp>
#include <stdexcept>

namespace looper {

//
// site graph type
//
//   g = 0 (vertically disconnected)
//

struct site_graph_type {
  static bool is_valid_gid(int g) { return g == 0; }
  static bool is_compatible(int /* g */, int /* c */) { return true; }
  static bool is_frozen(int /* g */) { return false; }
  template<class T>
  static boost::tuple<int /* curr */, int /* loop0 */, int /* loop1 */>
  reconnect(int /* g */, std::vector<T>& fragments, int curr) {
    int loop1 = add(fragments);
    return boost::make_tuple(loop1, curr, loop1);
  }
};


//
// bond graph types
//
//   g = 0   o.o    g = 1   o o
//                           X
//           o.o            o o
//
//   g = 2   o.o    g = 3   o-o
//           |*|            |X|
//           o.o            o-o
//
//   g = 4   o o    g = 5   o-o
//            *
//           o o            o-o

// optimized bond_graph_type for Ising model

struct ising_bond_graph_type {
  static bool is_valid_gid(int g) { return g == 2 || g == 3; }
  static bool is_compatible(int g, int c0, int c1) { return (g & 1) ^ c0 ^ c1; }
  static bool is_frozen(int /* g */) { return true; }
  template<class T>
  static boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
  reconnect(int g, std::vector<T>& fragments, int curr0, int curr1) {
    int loop = unify(fragments, curr0, curr1);
    return boost::make_tuple(loop, loop, loop, loop);
  }

  static bool has_nontrivial_diagonal_choice_helper() { return false; }
  struct diagonal_choice_helper {
    diagonal_choice_helper() {}
    template<std::size_t N>
    diagonal_choice_helper(boost::array<double, N> const&) {}
  };
  template<typename RNG>
  static int choose_diagonal(RNG& rng, std::vector<diagonal_choice_helper> const&, int,
    int c0, int c1) {
    return (c0 == c1) ? 3 : 2;
  }

  static bool has_nontrivial_offdiagonal_choice_helper() { return false; }
  struct offdiagonal_choice_helper {
    offdiagonal_choice_helper() {}
    template<std::size_t N>
    offdiagonal_choice_helper(boost::array<double, N> const&) {}
  };
  template<typename RNG>
  static int choose_offdiagonal(RNG&, std::vector<offdiagonal_choice_helper> const&, int, int,
    int) {
    boost::throw_exception(std::logic_error("ising_bond_graph_type::choose_offdiagonal"));
    return 0;
  }
};

// optimized bond_graph_type for Heisenberg Antiferromagnet

struct haf_bond_graph_type {
  static bool is_valid_gid(int g) { return g == 0; }
  static bool is_compatible(int /* g */, int c0, int c1) { return c0 != c1; }
  static bool is_frozen(int /* g */) { return false; }
  template<class T>
  static boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
  reconnect(int g, std::vector<T>& fragments, int curr0, int curr1) {
    int loop0 = unify(fragments, curr0, curr1);
    int loop1 = curr0 = curr1 = add(fragments);
    return boost::make_tuple(curr0, curr1, loop0, loop1);
  }

  static bool has_nontrivial_diagonal_choice_helper() { return false; }
  struct diagonal_choice_helper {
    diagonal_choice_helper() {}
    template<std::size_t N>
    diagonal_choice_helper(boost::array<double, N> const&) {}
  };
  template<typename RNG>
  static int choose_diagonal(RNG&, std::vector<diagonal_choice_helper> const&, int, int, int) {
    return 0;
  }

  static bool has_nontrivial_offdiagonal_choice_helper() { return false; }
  struct offdiagonal_choice_helper {
    offdiagonal_choice_helper() {}
    template<std::size_t N>
    offdiagonal_choice_helper(boost::array<double, N> const&) {}
  };
  template<typename RNG>
  static int choose_offdiagonal(RNG&, std::vector<offdiagonal_choice_helper> const&, int, int,
    int) {
    return 0;
  }
};

// optimized bond_graph_type for Heisenberg Ferromagnet

struct hf_bond_graph_type {
  static bool is_valid_gid(int g) { return g == 1; }
  static bool is_compatible(int c0, int c1) { return c0 == c1; }
  static bool is_frozen(int /* g */) { return false; }
  template<class T>
  static boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
  reconnect(int g, std::vector<T>& fragments, int curr0, int curr1) {
    int loop0 = curr0;
    int loop1 = curr1;
    std::swap(curr0, curr1);
    return boost::make_tuple(curr0, curr1, loop0, loop1);
  }

  static bool has_nontrivial_diagonal_choice_helper() { return false; }
  struct diagonal_choice_helper {
    diagonal_choice_helper() {}
    template<std::size_t N>
    diagonal_choice_helper(boost::array<double, N> const&) {}
  };
  template<typename RNG>
  static int choose_diagonal(RNG&, std::vector<diagonal_choice_helper> const&, int, int, int) {
    return 1;
  }

  static bool has_nontrivial_offdiagonal_choice_helper() { return false; }
  struct offdiagonal_choice_helper {
    offdiagonal_choice_helper() {}
    template<std::size_t N>
    offdiagonal_choice_helper(boost::array<double, N> const&) {}
  };
  template<typename RNG>
  static int choose_offdiagonal(RNG&, std::vector<offdiagonal_choice_helper> const&, int, int,
    int) {
    return 1;
  }
};

// bond_graph_type for XXZ interaction

struct xxz_bond_graph_type {
  static bool is_valid_gid(int g) { return g >= 0 && g <= 3; }
  static bool is_compatible(int g, int c0, int c1) { return (g & 1) ^ c0 ^ c1; }
  static bool is_frozen(int g) { return ((g & 2) == 2); }
  template<class T>
  static boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
  reconnect(int g, std::vector<T>& fragments, int curr0, int curr1) {
    int loop0, loop1;
    if ((g & 2) == 2) {
      loop0 = loop1 = curr0 = curr1 = unify(fragments, curr0, curr1);
    } else if (g == 0) {
      loop0 = unify(fragments, curr0, curr1);
      loop1 = curr0 = curr1 = add(fragments);
    } else {
      loop0 = curr0;
      loop1 = curr1;
      std::swap(curr0, curr1);
    }
    return boost::make_tuple(curr0, curr1, loop0, loop1);
  }

  static bool has_nontrivial_diagonal_choice_helper() { return true; }
  struct diagonal_choice_helper {
    diagonal_choice_helper() {}
    template<std::size_t N>
    diagonal_choice_helper(boost::array<double, N> const& v) {
      p[0] = alps::numeric::is_nonzero<1>(v[0] + v[2]) ? v[0] / (v[0] + v[2]) : 1; // for antiparallel
      p[1] = alps::numeric::is_nonzero<1>(v[1] + v[3]) ? v[1] / (v[1] + v[3]) : 1; // for parallel
    }
    boost::array<double, 2> p;
  };
  template<typename RNG>
  static int choose_diagonal(RNG& rng, std::vector<diagonal_choice_helper> const& helpers, int b,
    int c0, int c1) {
    int c = 1 ^ c0 ^ c1; // 0 for antiparallel, 1 for parallel,
    return (rng() < helpers[b].p[c] ? 0 : 2) ^ c;
  }

  static bool has_nontrivial_offdiagonal_choice_helper() { return true; }
  struct offdiagonal_choice_helper {
    offdiagonal_choice_helper() {}
    template<std::size_t N>
    offdiagonal_choice_helper(boost::array<double, N> const& v) {
      p = alps::numeric::is_nonzero<1>(v[0] + v[1]) ? v[0] / (v[0] + v[1]) : 1;
    }
    double p;
  };
  template<typename RNG>
  static int choose_offdiagonal(RNG& rng, std::vector<offdiagonal_choice_helper> const& helpers,
    int b, int /* c0 */, int /* c1 */) {
    return rng() < helpers[b].p ? 0 : 1;
  }
};

// bond_graph_type for XYZ interaction

struct xyz_bond_graph_type {
  static bool is_valid_gid(int g) { return g >= 0 && g <= 5; }
  static bool is_compatible(int g, int c0, int c1) { return (g & 1) ^ c0 ^ c1; }
  static bool is_frozen(int g) { return ((g & 2) == 2); }
  template<typename T>
  static boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
  reconnect(int g, std::vector<T>& fragments, int curr0, int curr1) {
    int loop0, loop1;
    if ((g & 2) == 2) {
      loop0 = loop1 = curr0 = curr1 = unify(fragments, curr0, curr1);
    } else if (g == 0 || g == 5) {
      loop0 = unify(fragments, curr0, curr1);
      loop1 = curr0 = curr1 = add(fragments);
    } else {
      loop0 = curr0;
      loop1 = curr1;
      std::swap(curr0, curr1);
    }
    return boost::make_tuple(curr0, curr1, loop0, loop1);
  }

  static bool has_nontrivial_diagonal_choice_helper() { return true; }
  struct diagonal_choice_helper {
    diagonal_choice_helper() {}
    template<std::size_t N>
    diagonal_choice_helper(boost::array<double, N> const& v) {
      // for antiparallel
      p[0] = alps::numeric::is_nonzero<1>(v[0] + v[2] + v[4]) ? v[0] / (v[0] + v[2] + v[4]) : 1;
      p[2] = alps::numeric::is_nonzero<1>(v[0] + v[2] + v[4]) ? (v[0] + v[2])/ (v[0] + v[2] + v[4]) : 1;
      // for parallel
      p[1] = alps::numeric::is_nonzero<1>(v[1] + v[3] + v[5]) ? v[1] / (v[1] + v[3] + v[5]) : 1;
      p[3] = alps::numeric::is_nonzero<1>(v[1] + v[3] + v[5]) ? (v[1] + v[3]) / (v[1] + v[3] + v[5]) : 1;
    }
    boost::array<double, 4> p;
  };
  template<typename RNG>
  static int choose_diagonal(RNG& rng, std::vector<diagonal_choice_helper> const& helpers, int b,
    int c0, int c1) {
    double r = rng();
    if (c0 != c1) {
      if (r < helpers[b].p[0])
        return 0;
      else if (r < helpers[b].p[2])
        return 2;
      else
        return 4;
    } else {
      if (r < helpers[b].p[1])
        return 1;
      else if (r < helpers[b].p[3])
        return 3;
      else
        return 5;
    }
  }

  static bool has_nontrivial_offdiagonal_choice_helper() { return true; }
  struct offdiagonal_choice_helper {
    offdiagonal_choice_helper() {}
    template<std::size_t N>
    offdiagonal_choice_helper(boost::array<double, N> const& v) {
      p[0] = alps::numeric::is_nonzero<1>(v[0] + v[1]) ? v[0] / (v[0] + v[1]) : 1;
      p[1] = alps::numeric::is_nonzero<1>(v[4] + v[5]) ? v[4] / (v[4] + v[5]) : 1;
    }
    boost::array<double, 2> p;
  };
  template<typename RNG>
  static int choose_offdiagonal(RNG& rng, std::vector<offdiagonal_choice_helper> const& helpers,
    int b, int c0, int c1) {
    double r = rng();
    if (c0 != c1)
      return (r < helpers[b].p[0] ? 0 : 1);
    else
      return (r < helpers[b].p[1] ? 4 : 5);
  }
};

template<typename SITE, typename BOND, typename LOC>
class local_graph {
public:
  typedef SITE site_graph_t;
  typedef BOND bond_graph_t;
  typedef LOC location_t;

  local_graph() {}
  local_graph(int type, const location_t& loc) : type_(type), loc_(loc) {}

  int type() const { return type_; }
  bool is_compatible(int c) const {
#ifndef NDEBUG
    if (!is_site())
      boost::throw_exception(std::logic_error("local_graph<>::is_compatible"));
#endif
    return site_graph_t::is_compatible(type_, c);
  }
  bool is_compatible(int c0, int c1) const {
#ifndef NDEBUG
    if (!is_bond())
      boost::throw_exception(std::logic_error("local_graph<>::is_compatible"));
#endif
    return bond_graph_t::is_compatible(type_, c0, c1);
  }

  const location_t& loc() const { return loc_; }
  int pos() const { return loc_.pos(); }
  int source() const { return loc_.source(); }
  int target() const { return loc_.target(); }
  bool is_bond() const { return loc_.is_bond(); }
  bool is_site() const { return loc_.is_site(); }

  template<class T>
  boost::tuple<int /* curr */, int /* loop0 */, int /* loop1 */>
  reconnect(std::vector<T>& fragments, int curr) const {
#ifndef NDEBUG
    if (!is_site())
      boost::throw_exception(std::logic_error("local_graph<>::reconnect"));
#endif
    return site_graph_t::reconnect(type_, fragments, curr);
  }

  template<class T>
  boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
  reconnect(std::vector<T>& fragments, int curr0, int curr1) const {
#ifndef NDEBUG
    if (!is_bond())
      boost::throw_exception(std::logic_error("local_graph<>::reconnect"));
#endif
    return bond_graph_t::reconnect(type_, fragments, curr0, curr1);
  }

  static int loop_l(int /* t */, int loop0, int /* loop1 */) { return loop0; }
  static int loop_u(int /* t */, int /* loop0 */, int loop1) { return loop1; }
  static int loop_l0(int /* t */, int loop0, int /* loop1 */) { return loop0; }
  static int loop_l1(int t, int loop0, int loop1) { return (t ==0 || t == 5) ? loop0 : loop1; }
  static int loop_u0(int /* t */, int /* loop0 */, int loop1) { return loop1; }
  static int loop_u1(int t, int loop0, int loop1) { return (t ==0 || t == 5) ? loop1 : loop0; }

private:
  int type_;
  location_t loc_;
};

template<typename SITE, typename BOND, typename LOC>
inline int type(const local_graph<SITE, BOND, LOC>& g) { return g.type(); }
template<typename SITE, typename BOND, typename LOC>
inline bool is_compatible(const local_graph<SITE, BOND, LOC>& g, int c) {
  return g.is_compatible(c);
}
template<typename SITE, typename BOND, typename LOC>
inline bool is_compatible(const local_graph<SITE, BOND, LOC>& g, int c0, int c1) {
  return g.is_compatible(c0, c1);
}

template<typename SITE, typename BOND, typename LOC>
inline int pos(const local_graph<SITE, BOND, LOC>& g) { return g.pos(); }
template<typename SITE, typename BOND, typename LOC>
inline bool is_site(const local_graph<SITE, BOND, LOC>& g) { return g.is_site(); }
template<typename SITE, typename BOND, typename LOC>
inline bool is_bond(const local_graph<SITE, BOND, LOC>& g) { return g.is_bond(); }

template<typename T, typename SITE, typename BOND, typename LOC>
boost::tuple<int /* curr */, int /* loop0 */, int /* loop1 */>
reconnect(std::vector<T>& fragments, const local_graph<SITE, BOND, LOC>& g, int curr) {
  return g.reconnect(fragments, curr);
}

template<typename T, typename SITE, typename BOND, typename LOC>
boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
reconnect(std::vector<T>& fragments, const local_graph<SITE, BOND, LOC>& g, int curr0, int curr1) {
  return g.reconnect(fragments, curr0, curr1);
}


//
// graph_chooser
//

template<typename LOCAL_GRAPH> class graph_chooser;

template<typename SITE, typename BOND, typename LOC>
class graph_chooser<local_graph<SITE, BOND, LOC> > {
public:
  typedef SITE site_graph_t;
  typedef BOND bond_graph_t;
  typedef LOC location_t;
  typedef local_graph<site_graph_t, bond_graph_t, location_t> local_graph_t;

  graph_chooser() : dist_graph_() {}
  template<class WEIGHT_TABLE>
  graph_chooser(const WEIGHT_TABLE& wt, bool is_path_integral) : dist_graph_() {
    init(wt, is_path_integral);
  }

  template<class WEIGHT_TABLE>
  void init(const WEIGHT_TABLE& wt, bool is_path_integral) {
    graph_.resize(0);
    diag_.resize(0);
    offdiag_.resize(0);
    std::vector<double> w;
    weight_ = 0;
    BOOST_FOREACH(typename WEIGHT_TABLE::site_weight_t const& sw, wt.site_weights()) {
      for (int g = 0; g < WEIGHT_TABLE::num_site_graphs; ++g) {
        if (alps::numeric::is_nonzero<1>(sw.second.v[g])) {
          if (site_graph_t::is_valid_gid(g)) {
            graph_.push_back(local_graph_t(g, location_t::site_location(sw.first)));
            w.push_back(sw.second.v[g]);
            weight_ += sw.second.v[g];
          } else {
            boost::throw_exception(std::runtime_error("graph_chooser::init"));
          }
        }
      }
    }
    BOOST_FOREACH(typename WEIGHT_TABLE::bond_weight_t const& bw, wt.bond_weights()) {
      for (int g = 0; g < WEIGHT_TABLE::num_bond_graphs; ++g) {
        if (alps::numeric::is_nonzero<1>(bw.second.v[g])) {
          if (bond_graph_t::is_valid_gid(g)) {
            graph_.push_back(local_graph_t(g, location_t::bond_location(bw.first)));
            w.push_back(bw.second.v[g]);
            weight_ += bw.second.v[g];
          } else {
            boost::throw_exception(std::runtime_error("graph_chooser::init"));
          }
        }
      }
      if (!is_path_integral && bond_graph_t::has_nontrivial_diagonal_choice_helper())
        diag_.push_back(typename bond_graph_t::diagonal_choice_helper(bw.second.v));
      if (bond_graph_t::has_nontrivial_offdiagonal_choice_helper())
        offdiag_.push_back(typename bond_graph_t::offdiagonal_choice_helper(bw.second.v));
    }
    if (w.size()) dist_graph_.init(w);
    if (alps::numeric::is_zero<2>(weight_)) weight_ = 1.0e-20;
  }

  template<typename RNG>
  const local_graph_t& graph(RNG& rng) const { return graph_[dist_graph_(rng)]; }
  template<typename RNG>
  local_graph_t diagonal(RNG& /* rng */, const location_t& loc, int /* c */) const {
    return local_graph_t(0, loc);
  }
  template<typename RNG>
  local_graph_t diagonal(RNG& rng, const location_t& loc, int c0, int c1) const {
    return local_graph_t(bond_graph_t::choose_diagonal(rng, diag_, pos(loc), c0, c1), loc);
  }
  template<typename RNG>
  local_graph_t offdiagonal(RNG& rng, const location_t& loc, int c0, int c1) const {
    return local_graph_t(bond_graph_t::choose_offdiagonal(rng, offdiag_, pos(loc), c0, c1), loc);
  }
  double weight() const { return weight_; }

private:
  random_choice_walker_d<> dist_graph_;
  double weight_;
  std::vector<local_graph_t> graph_;
  std::vector<typename bond_graph_t::diagonal_choice_helper> diag_;
  std::vector<typename bond_graph_t::offdiagonal_choice_helper> offdiag_;
};

} // end namespace looper

#endif // LOOPER_GARPH_IMPL_H
