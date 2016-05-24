/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_LATTICE_H_
#define LOOPER_LATTICE_H_

#include "alternating_tensor.h"
#include <alps/lattice.h>
#include <alps/parapack/integer_range.h>
#include <boost/throw_exception.hpp>
#include <stdexcept>
#include <utility>                 // std::pair, std::make_pair
#include <vector>                  // std::vector

//
// extension to boost
//

namespace boost {

template<class T0, class T1, class T2, class T3, class T4, class T5, class T6>
typename boost::graph_traits<
  boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::vertex_descriptor
source(int b, boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> const& g) {
  return boost::source(*(edges(g).first + b), g);
}

template<class T0, class T1, class T2, class T3, class T4, class T5, class T6>
typename boost::graph_traits<
  boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::vertex_descriptor
target(int b, boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> const& g) {
  return boost::target(*(edges(g).first + b), g);
}

template<class T0, class T1, class T2, class T3, class T4, class T5, class T6>
typename boost::graph_traits<
  boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::edge_descriptor
edge(boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> const& g, int b) {
  return *(edges(g).first + b);
}

template<class T0, class T1, class T2, class T3, class T4, class T5, class T6>
typename boost::graph_traits<
  boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::edge_descriptor
bond(boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> const& g, int b) { return edge(g, b); }

} // end namespace boost

namespace looper {

using alps::site_index_t;
using alps::site_type_t;
using alps::coordinate_t;
using alps::parity_t;
struct gauge_t { typedef boost::vertex_property_tag kind; };

using alps::bond_index_t;
using alps::bond_type_t;
using alps::bond_vector_t;
using alps::bond_vector_relative_t;

using alps::graph_name_t;
using alps::dimension_t;

using alps::graph_traits;
using alps::coordinate_type;

struct real_vertex_t { typedef boost::vertex_property_tag kind; };
typedef real_vertex_t real_site_t;
struct real_edge_t { typedef boost::edge_property_tag kind; };
typedef real_edge_t real_bond_t;

template<typename RG>
struct virtual_graph {
private:
  typedef typename graph_traits<RG>::site_descriptor real_site_descriptor;
  typedef typename graph_traits<RG>::bond_descriptor real_bond_descriptor;

public:
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    // site property
    boost::property<real_site_t, real_site_descriptor,
    boost::property<gauge_t, double> >,
    // bond property
    boost::property<bond_index_t, unsigned int,
    boost::property<real_bond_t, real_bond_descriptor> >,
    // graph property
    boost::no_property,
    boost::vecS> type;
};

template<class T0, class T1, class T2, class T3, class T4, class T5, class T6>
int
gauge(boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> const& g,
  typename graph_traits<boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::site_descriptor s) {
  return 2 * get(parity_t(), g, s) - 1;
}

template<class G>
int
gauge(alps::graph_helper<G> const& g, typename graph_traits<G>::site_descriptor s) {
  return 2 * get(parity_t(), g.graph(), s) - 1;
}


//
// helper class: linear_distance_helper
//
// for calculating 'real' linear distance between vertices

template<typename G>
class linear_distance_helper {
public:
  typedef typename alps::graph_helper<G>::vector_type vector_type;
  typedef typename alps::graph_helper<G>::site_descriptor site_descriptor;

  linear_distance_helper(alps::graph_helper<G> const& lattice) : lattice_(lattice), span_(lattice.dimension(), 0) {
    int d = 0;
    BOOST_FOREACH(vector_type const& vec, lattice.basis_vectors()) {
      int ext = lattice.lattice().extent(d);
      for (int i = 0; i < lattice.dimension(); ++i) {
        span_[i] += ext * vec[i];
      }
      ++d;
    }
    for (int i = 0; i < lattice.dimension(); ++i) span_[i] = std::abs(span_[i]);
  }

  double distance(site_descriptor const& s, site_descriptor const& t) const {
    vector_type cs = lattice_.coordinate(s);
    vector_type ct = lattice_.coordinate(t);
    double distance = 0;
    for (int i = 0; i < lattice_.dimension(); ++i) {
      double x = std::abs(ct[i] - cs[i]);
      distance += std::pow(std::min(x, span_[i] - x), 2.);
    }
    return std::sqrt(distance);
  }

private:
  alps::graph_helper<G> const& lattice_;
  vector_type span_;
};


//
// class template virtual_mapping
//
// for describing a mapping from a real site/bond to virtual ones

template<typename RG, typename VG>
class virtual_mapping {
public:
  typedef RG real_graph_type;
  typedef typename graph_traits<real_graph_type>::site_descriptor real_site_descriptor;
  typedef typename graph_traits<real_graph_type>::bond_descriptor real_bond_descriptor;
  typedef typename graph_traits<real_graph_type>::site_iterator real_site_iterator;
  typedef typename graph_traits<real_graph_type>::bond_iterator real_bond_iterator;

  typedef VG virtual_graph_type;
  typedef typename graph_traits<virtual_graph_type>::site_iterator virtual_site_iterator;
  typedef typename graph_traits<virtual_graph_type>::bond_iterator virtual_bond_iterator;

  typedef std::pair<virtual_site_iterator, virtual_site_iterator> virtual_site_range_type;
  typedef std::pair<virtual_bond_iterator, virtual_bond_iterator> virtual_bond_range_type;

  virtual_mapping() :
    site_map_(1, virtual_site_iterator()),
    bond_map_(1, virtual_bond_iterator()),
    s2b_map_(1, virtual_bond_iterator()),
    s2b_offset_(0), max_vv_(0) {}

  virtual_site_range_type
  virtual_vertices(const real_graph_type& rg, const real_site_descriptor& rv) const {
    return std::make_pair(site_map_[get(site_index_t(), rg, rv)],
                          site_map_[get(site_index_t(), rg, rv) + 1]);
  }

  virtual_bond_range_type
  virtual_bonds(const real_graph_type& rg, const real_bond_descriptor& re) const {
    return std::make_pair(bond_map_[get(bond_index_t(), rg, re)],
                          bond_map_[get(bond_index_t(), rg, re) + 1]);
  }

  virtual_bond_range_type
  virtual_bonds(const real_graph_type& rg, const real_site_descriptor& rv) const {
    return std::make_pair(s2b_map_[get(site_index_t(), rg, rv)],
                          s2b_map_[get(site_index_t(), rg, rv) + 1]);
  }

  void add_vertices(const real_graph_type& rg, const real_site_descriptor& rv,
    const virtual_site_iterator& first, const virtual_site_iterator& last) {
    if (get(site_index_t(), rg, rv) != site_map_.size() - 1)
      boost::throw_exception(std::invalid_argument("virtual_mapping<G>::add_vertices()"));
    if (site_map_.size() == 1) {
      site_map_.back() = first;
    } else {
      assert(first == site_map_.back());
    }
    site_map_.push_back(last);
    max_vv_ = std::max(max_vv_, static_cast<int>(last - first));
  }

  void add_bonds(const real_graph_type& rg, const real_bond_descriptor& re,
    const virtual_bond_iterator& first, const virtual_bond_iterator& last) {
    if (get(bond_index_t(), rg, re) != bond_map_.size() - 1)
      boost::throw_exception(std::invalid_argument("virtual_mapping<G>::add_bonds()"));
    if (bond_map_.size() == 1) {
      bond_map_.back() = first;
    } else {
      assert(first == bond_map_.back());
    }
    bond_map_.push_back(last);
  }

  void add_s2bonds(const real_graph_type& rg, const real_site_descriptor& rv,
    const virtual_bond_iterator& first, const virtual_bond_iterator& last) {
    if (get(site_index_t(), rg, rv) != s2b_map_.size() - 1)
      boost::throw_exception(std::invalid_argument("virtual_mapping<G>::add_s2bonds()"));
    if (s2b_map_.size() == 1) {
      s2b_map_.back() = first;
    } else {
      assert(first == s2b_map_.back());
    }
    s2b_map_.push_back(last);
  }

  void clear() {
    site_map_.clear();
    site_map_.push_back(virtual_site_iterator());
    bond_map_.clear();
    bond_map_.push_back(virtual_bond_iterator());
    s2b_map_.clear();
    s2b_map_.push_back(virtual_bond_iterator());
    s2b_offset_ = 0;
    max_vv_ = 0;
  }

  void set_s2bond_type_offset(int t) { s2b_offset_ = t; }
  int s2bond_type_offset() const { return s2b_offset_; }

  int max_virtual_vertices() const { return max_vv_; }

  bool operator==(const virtual_mapping& rhs) const {
    return site_map_ == rhs.site_map_ && bond_map_ == rhs.bond_map_ && s2b_map_ == rhs.s2b_map_;
  }
  bool operator!=(const virtual_mapping& rhs) { return !(*this == rhs); }

  void output(std::ostream& os, const real_graph_type& rg, const virtual_graph_type& vg) const;

private:
  std::vector<virtual_site_iterator> site_map_;
  std::vector<virtual_bond_iterator> bond_map_;
  std::vector<virtual_bond_iterator> s2b_map_;
  int s2b_offset_;
  int max_vv_;
};


template<typename RG>
class lattice_helper {
public:
  typedef RG                                  real_graph_type;
  typedef typename virtual_graph<RG>::type    virtual_graph_type;
  typedef alps::graph_helper<real_graph_type> graph_helper_type;
  typedef virtual_mapping<real_graph_type, virtual_graph_type>
                                              mapping_type;
  typedef alps::integer_range<int>            type_range_type;

  typedef real_graph_type    rg_type;
  typedef virtual_graph_type vg_type;
  typedef mapping_type       mp_type;

  lattice_helper(alps::Parameters const& params) : helper_(params) {
    volume_ = calc_volume(params);
    convert_type(params);
  }
  template<typename M>
  lattice_helper(alps::Parameters const& params, M const& model, bool has_d_term = false) :
    helper_(params) {
    volume_ = calc_volume(params);
    convert_type(params);
    generate_virtual_graph(model, has_d_term);
  }

  double volume() const { return volume_; }

  real_graph_type const& rg() const { return helper_.graph(); }
  virtual_graph_type const& vg() const { return vgraph_; }
  mapping_type const& mp() const { return mapping_; }

  graph_helper_type const& graph_helper() const { return helper_; }

  template<typename M>
  void generate_virtual_graph(M const& model, bool has_d_term = false);

protected:
  void convert_type(alps::Parameters const& p);

  double calc_volume(alps::Parameters const& p) {
    double vol;
    if (helper_.has_lattice(p["LATTICE"])) {
      vol = helper_.volume();
      std::vector<std::vector<double> > vecs;
      BOOST_FOREACH(std::vector<double> const& v, helper_.basis_vectors()) vecs.push_back(v);
      double f = 0;
      switch (vecs.size()) {
      case 1 :
        f = vecs[0][0];
        break;
      case 2 :
        for (int i = 0; i < 2; ++i)
          for (int j = 0; j < 2; ++j)
            f += alternating_tensor(i,j,2) * vecs[0][i] * vecs[1][j];
        break;
      case 3 :
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
              f += alternating_tensor(i,j,k) * vecs[0][i] * vecs[1][j] * vecs[2][k];
        break;
      default :
        std::cout << "WARNING: can not calculate unit lattice volume.  Assume 1.\n";
        f = 1;
      }
      vol *= std::abs(f);
    } else {
      vol = num_sites(rg());
    }
    return vol;
  }

private:
  graph_helper_type helper_;
  double volume_;
  virtual_graph_type vgraph_;
  mapping_type mapping_;
  type_range_type site_type_range_;
  type_range_type bond_type_range_;
};


template<typename L>
struct real_bond_iterator {
  typedef typename graph_traits<typename L::rg_type>::bond_iterator type;
};

template<typename L>
struct real_site_iterator {
  typedef typename graph_traits<typename L::rg_type>::site_iterator type;
};

template<typename L>
struct real_bond_descriptor {
  typedef typename graph_traits<typename L::rg_type>::bond_descriptor type;
};

template<typename L>
struct real_site_descriptor {
  typedef typename graph_traits<typename L::rg_type>::site_descriptor type;
};

template<typename L>
struct virtual_bond_iterator {
  typedef typename graph_traits<typename L::vg_type>::bond_iterator type;
};

template<typename L>
struct virtual_site_iterator {
  typedef typename graph_traits<typename L::vg_type>::site_iterator type;
};

template<typename L>
struct virtual_bond_descriptor {
  typedef typename graph_traits<typename L::vg_type>::bond_descriptor type;
};

template<typename L>
struct virtual_site_descriptor {
  typedef typename graph_traits<typename L::vg_type>::site_descriptor type;
};

template<typename L>
struct momentum_iterator {
  typedef typename L::graph_helper_type::momentum_iterator type;
};


template<typename RG>
std::pair<
  typename virtual_site_iterator<lattice_helper<RG> >::type,
  typename virtual_site_iterator<lattice_helper<RG> >::type>
sites(lattice_helper<RG> const& lat,
  typename real_site_descriptor<lattice_helper<RG> >::type s) {
  return lat.mp().virtual_vertices(lat.rg(), s);
}

template<typename RG>
std::pair<
  typename virtual_bond_iterator<lattice_helper<RG> >::type,
  typename virtual_bond_iterator<lattice_helper<RG> >::type>
bonds(lattice_helper<RG> const& lat,
  typename real_bond_descriptor<lattice_helper<RG> >::type b) {
  return lat.mp().virtual_bonds(lat.rg(), b);
}

template<typename RG>
std::pair<
  typename virtual_bond_iterator<lattice_helper<RG> >::type,
  typename virtual_bond_iterator<lattice_helper<RG> >::type>
bonds(lattice_helper<RG> const& lat,
  typename real_site_descriptor<lattice_helper<RG> >::type s) {
  return lat.mp().virtual_bonds(lat.rg(), s);
}

template<typename RG>
int max_virtual_sites(lattice_helper<RG> const& lat) {
  return lat.mp().max_virtual_vertices();
}

//
// wrappers
//

template<typename RG>
std::vector<std::string>
distance_labels(lattice_helper<RG> const& lat, boost::optional<int> const& origin) {
  if (origin) {
    std::vector<std::string> label;
    for (int s = 0; s < num_sites(lat.rg()); ++s)
      label.push_back(lat.graph_helper().coordinate_string(origin.get()) +
                      " -- " + lat.graph_helper().coordinate_string(s));
    return label;
  } else {
    return lat.graph_helper().distance_labels();
  }
}

template<typename RG>
int distance(lattice_helper<RG> const& lat, int s0, int s1) {
  return lat.graph_helper().distance(s0, s1);
}

template<typename RG>
std::vector<unsigned int> distance_multiplicities(lattice_helper<RG> const& lat) {
  std::vector<unsigned int> mult = lat.graph_helper().distance_multiplicities();
  return lat.graph_helper().distance_multiplicities();
}

template<typename RG>
bool is_bipartite(lattice_helper<RG> const& lat) {
  return lat.graph_helper().is_bipartite();
}

template<typename RG>
std::vector<std::string> momenta_labels(lattice_helper<RG> const& lat) {
  return lat.graph_helper().momenta_labels();
}

template<typename RG>
std::pair<
  typename momentum_iterator<lattice_helper<RG> >::type,
  typename momentum_iterator<lattice_helper<RG> >::type>
momenta(lattice_helper<RG> const& lat) {
  return lat.graph_helper().momenta();
}

template<typename RG>
typename lattice_helper<RG>::graph_helper_type::vector_type
momentum(typename lattice_helper<RG>::graph_helper_type::vector_type const& m, lattice_helper<RG> const& lat) {
    return lat.graph_helper().momentum(m);
}

} // end namespace looper

namespace looper {

///////////////////////////////////////////////////////////////////////////////
//
// Implementations
//

//
// class template virtual_mapping
//

template<typename RG, typename VG>
void virtual_mapping<RG, VG>::output(std::ostream& os, RG const& rg, VG const& vg) const {
  os << "[[vitual_mapping]]\n";
  os << "  number of site groups = " << num_vertices(rg) << '\n';
  os << "  site mapping:\n";
  BOOST_FOREACH(real_site_descriptor v, vertices(rg)) {
    os << "    " << get(site_index_t(), rg, v) << " -> ";
    virtual_site_range_type vr = virtual_vertices(rg, v);
    if (vr.first == vr.second) {
      os << "null\n";
    } else if (vr.first == boost::prior(vr.second)) {
      os << get(site_index_t(), vg, *vr.first) << '\n';
    } else {
      os << '['
         << get(site_index_t(), vg, *vr.first) << ','
         << get(site_index_t(), vg, *boost::prior(vr.second))
         << "]\n";
    }
  }
  os << "  number of bond groups = " << num_bonds(rg) << '\n';
  os << "  bond mapping:\n";
  BOOST_FOREACH(real_bond_descriptor e, bonds(rg)) {
    os << "    " << get(bond_index_t(), rg, e) << " -> ";
    virtual_bond_range_type er = virtual_bonds(rg, e);
    if (er.first == er.second) {
      os << "null\n";
    } else if (er.first == boost::prior(er.second)) {
      os << get(bond_index_t(), vg, *er.first) << '\n';
    } else {
      os << '['
         << get(bond_index_t(), vg, *er.first) << ','
         << get(bond_index_t(), vg, *boost::prior(er.second))
         << "]\n";
    }
  }
  os << "  site2bond mapping:\n";
  BOOST_FOREACH(real_site_descriptor v, vertices(rg)) {
    os << "    " << get(site_index_t(), rg, v) << " -> ";
    virtual_bond_range_type er = virtual_bonds(rg, v);
    if (er.first == er.second) {
      os << "null\n";
    } else if (er.first == boost::prior(er.second)) {
      os << get(bond_index_t(), vg, *er.first) << '\n';
    } else {
      os << '['
         << get(bond_index_t(), vg, *er.first) << ','
         << get(bond_index_t(), vg, *boost::prior(er.second))
         << "]\n";
    }
  }
}


//
// class template lattice_helper
//

template<typename RG>
template<typename M>
void lattice_helper<RG>::generate_virtual_graph(M const& model, bool has_d_term) {
  vgraph_.clear();
  mapping_.clear();

  // setup s2bond_type_offset
  int tmin = 0;
  typename alps::property_map<site_type_t, const rg_type, int>::type
    site_type = alps::get_or_default(site_type_t(), rg(), 0);
  BOOST_FOREACH(typename graph_traits<rg_type>::site_descriptor rv, vertices(rg()))
    tmin = std::min(tmin, int(site_type[rv]));

  int tmax = 0;
  typename alps::property_map<bond_type_t, const rg_type, int>::type
    bond_type = alps::get_or_default(bond_type_t(), rg(), 0);
  BOOST_FOREACH(typename graph_traits<rg_type>::bond_descriptor re, bonds(rg()))
    tmax = std::max(tmax, int(bond_type[re]));
  mapping_.set_s2bond_type_offset(tmax-tmin+1);

  // add vertices to virtual graph
  BOOST_FOREACH(typename graph_traits<rg_type>::site_descriptor rv, vertices(rg())) {
    for (int i = 0; i < model.site(rv, rg()).s.get_twice(); ++i) {
      typename graph_traits<vg_type>::site_descriptor vvd = add_vertex(vgraph_);
      put(real_site_t(), vgraph_, vvd, rv);
      put(gauge_t(), vgraph_, vvd, 2. * get(parity_t(), rg(), rv) - 1);
    }
  }

  // setup site mapping
  typename graph_traits<vg_type>::site_iterator vvi_first = vertices(vgraph_).first;
  typename graph_traits<vg_type>::site_iterator vvi_last = vvi_first;
  BOOST_FOREACH(typename graph_traits<rg_type>::site_descriptor rv, vertices(rg())) {
    vvi_last += model.site(rv, rg()).s.get_twice();
    mapping_.add_vertices(rg(), rv, vvi_first, vvi_last);
    vvi_first = vvi_last;
  }

  // add bonds to virtual graph
  BOOST_FOREACH(typename graph_traits<rg_type>::bond_descriptor re, bonds(rg())) {
    typename graph_traits<rg_type>::site_descriptor rs = source(re, rg());
    typename graph_traits<rg_type>::site_descriptor rt = target(re, rg());
    typename graph_traits<vg_type>::site_iterator vvsi, vvsi_end;
    for (boost::tie(vvsi, vvsi_end) = mapping_.virtual_vertices(rg(), rs);
         vvsi != vvsi_end; ++vvsi) {
      typename graph_traits<vg_type>::site_iterator vvti, vvti_end;
      for (boost::tie(vvti, vvti_end) = mapping_.virtual_vertices(rg(), rt);
           vvti != vvti_end; ++vvti) {
        typename graph_traits<vg_type>::bond_descriptor
          ved = add_edge(*vvsi, *vvti, vgraph_).first;
        put(bond_index_t(), vgraph_, ved, num_bonds(vgraph_) - 1);
        put(real_bond_t(), vgraph_, ved, re);
      }
    }
  }

  // add `in-real-site' bonds to virtual graph
  if (has_d_term) {
    BOOST_FOREACH(typename graph_traits<rg_type>::site_descriptor rv, vertices(rg())) {
      typename graph_traits<vg_type>::site_iterator vvsi, vvsi_end;
      for (boost::tie(vvsi, vvsi_end) =
             mapping_.virtual_vertices(rg(), rv); vvsi != vvsi_end;
           ++vvsi)
        for (typename graph_traits<vg_type>::site_iterator
               vvti = boost::next(vvsi); vvti != vvsi_end; ++vvti) {
          typename graph_traits<vg_type>::bond_descriptor
            ved = add_edge(*vvsi, *vvti, vgraph_).first;
          put(bond_index_t(), vgraph_, ved, num_bonds(vgraph_) - 1);
        }
    }
  }

  // setup bond and s2bond mapping
  typename graph_traits<vg_type>::bond_iterator vei_first = bonds(vgraph_).first;
  typename graph_traits<vg_type>::bond_iterator vei_last = vei_first;
  // BOOST_FOREACH(typename graph_traits<rg_type>::bond_descriptor re, bonds(rg())) {
  typename graph_traits<rg_type>::bond_iterator rei, rei_end;
  for (boost::tie(rei, rei_end) = bonds(rg()); rei != rei_end; ++rei) {
    typename graph_traits<rg_type>::bond_descriptor re = *rei;
    vei_last += model.site(source(re, rg()), rg()).s.get_twice() *
      model.site(target(re, rg()), rg()).s.get_twice();
    mapping_.add_bonds(rg(), re, vei_first, vei_last);
    vei_first = vei_last;
  }
  BOOST_FOREACH(typename graph_traits<rg_type>::site_descriptor rv, vertices(rg())) {
    if (has_d_term)
      vei_last += model.site(rv, rg()).s.get_twice() *
        (model.site(rv, rg()).s.get_twice() - 1) / 2;
    mapping_.add_s2bonds(rg(), rv, vei_first, vei_last);
    vei_first = vei_last;
  }

  // set range of types
  site_type_range_ = 0;
  BOOST_FOREACH(typename graph_traits<rg_type>::site_descriptor rv, vertices(rg()))
    site_type_range_.include(get(site_index_t(), rg(), rv));
  bond_type_range_ = 0;
  BOOST_FOREACH(typename graph_traits<rg_type>::bond_descriptor re, bonds(rg()))
    bond_type_range_.include(get(bond_index_t(), rg(), re));
}

template<typename RG>
void lattice_helper<RG>::convert_type(alps::Parameters const& p) {
  typedef typename graph_traits<RG>::site_descriptor site_descriptor;
  typedef typename graph_traits<RG>::bond_descriptor bond_descriptor;
  if (p.value_or_default("USE_SITE_INDICES_AS_TYPES", false)) {
    BOOST_FOREACH(site_descriptor s, sites(helper_.graph()))
      put(site_type_t(), helper_.graph(), s, s);
  }
  if (p.value_or_default("USE_BOND_INDICES_AS_TYPES", false)) {
    BOOST_FOREACH(bond_descriptor b, bonds(helper_.graph()))
      put(bond_type_t(), helper_.graph(), b, get(bond_index_t(), helper_.graph(), b));
  }
}

} // end namespace looper

#endif // LOOPER_LATTICE_H
