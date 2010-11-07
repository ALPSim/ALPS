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

#ifndef LOOPER_MODEL_PARAMETER_H
#define LOOPER_MODEL_PARAMETER_H

#include "lapack.h"
#include "lattice.h"
#include "matrix.h"
#include <alps/numeric/is_equal.hpp>
#include <alps/numeric/is_positive.hpp>
#include <alps/numeric/round.hpp>
#include <alps/model.h>

// IBM workaround
#ifdef hz
#undef hz
#endif

namespace looper {

//
// site_parameter
//

struct site_parameter {
  typedef alps::half_integer<int> spin_type;

  spin_type s; // size of spin
  double c;    // constant energy offset: + c
  double hx;   // transverse field:       - hx * Sx
  double hz;   // longitudinal field:     - hz * Sz
  double d;    // single-ion anisotropy:  + d * (Sz)^2

  site_parameter(double s_in = 0.5) : s(s_in), c(0), hx(0), hz(0), d(0) {}
  site_parameter(double s_in, double c_in, double hx_in, double hz_in, double d_in) :
    s(s_in), c(c_in), hx(hx_in), hz(hz_in), d(d_in) {}
  template<typename J>
  site_parameter(const alps::half_integer<J>& s_in) : s(s_in), c(0), hx(0), hz(0), d(0) {}
  template<typename J>
  site_parameter(const alps::half_integer<J>& s_in, double c_in, double hx_in, double hz_in,
    double d_in) : s(s_in), c(c_in), hx(hx_in), hz(hz_in), d(d_in) {}
  site_parameter(const boost::multi_array<double, 2>& mat) { do_fit(mat); }

  bool operator==(const site_parameter& rhs) const {
    return (s == rhs.s) && alps::numeric::is_equal<1>(c, rhs.c) && alps::numeric::is_equal<1>(hx, rhs.hx) &&
      alps::numeric::is_equal<1>(hz, rhs.hz) && alps::numeric::is_equal<1>(d, rhs.d);
  }
  bool operator!=(const site_parameter& rhs) const { return !operator==(rhs); }

  bool is_quantal() const { return alps::numeric::is_nonzero<1>(hx); }
  bool has_hz() const { return alps::numeric::is_nonzero<1>(hz); }
  bool has_d() const { return alps::numeric::is_nonzero<1>(d); }

  void output(std::ostream& os) const {
    os << '(' << s << ", " << c << ", " << hx << ", " << hz << ", " << d << ')';
  }

protected:
  void do_fit(const boost::multi_array<double, 2>& mat);
};

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

inline std::ostream& operator<<(std::ostream& os, looper::site_parameter const& p) {
  p.output(os);
  return os;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

namespace looper {

//
// bond_paraemter
//

struct bond_parameter_xxz {
  double c;   // constant energy offset:   + c
  double jxy; // off-diagonal interaction: + jxy/2 * (s1+ s2- + hc)
  double jz;  // diagonal interaction:     + jz s1z s2z

  bond_parameter_xxz() : c(0), jxy(0), jz(0) {}
  bond_parameter_xxz(double c_in, double jxy_in, double jz_in) : c(c_in), jxy(jxy_in), jz(jz_in) {}
  bond_parameter_xxz(const boost::multi_array<double, 4>& mat) { do_fit(mat); }

  bool operator==(const bond_parameter_xxz& rhs) const {
    return alps::numeric::is_equal<1>(c, rhs.c) && alps::numeric::is_equal<1>(jxy, rhs.jxy) &&
      alps::numeric::is_equal<1>(jz, rhs.jz);
  }
  bool operator!=(const bond_parameter_xxz& rhs) const { return !operator==(rhs); }

  bool is_quantal() const { return alps::numeric::is_nonzero<1>(jxy); }

  void print(std::ostream& os) const { os << '(' << c << ", " << jxy << ", " << jz << ')'; }

protected:
  void do_fit(const boost::multi_array<double, 4>& mat);
};

struct bond_parameter_xyz {
  double c;   // constant energy offset:   + c
  double jx;  // off-diagonal interaction: + jx S1x S2x
  double jy;  // off-diagonal interaction: + jy S1y S2y
  double jz;  // diagonal interaction:     + jz S1z S2z

  bond_parameter_xyz() : c(0), jx(0), jy(0), jz(0) {}
  bond_parameter_xyz(double c_in, double jx_in, double jy_in, double jz_in) :
    c(c_in), jx(jx_in), jy(jy_in), jz(jz_in) {}
  bond_parameter_xyz(const boost::multi_array<double, 4>& mat) { do_fit(mat); }

  bool operator==(const bond_parameter_xyz& rhs) const {
    return alps::numeric::is_equal<1>(c, rhs.c) && alps::numeric::is_equal<1>(jx, rhs.jx) &&
      alps::numeric::is_equal<1>(jy, rhs.jy) && alps::numeric::is_equal<1>(jz, rhs.jz);
  }
  bool operator!=(const bond_parameter_xyz& rhs) const { return !operator==(rhs); }

  bool is_quantal() const { return alps::numeric::is_nonzero<1>(jx) || alps::numeric::is_nonzero<1>(jy); }

  void print(std::ostream& os) const {
    os << '(' << c << ", " << jx << ", " << jy << ", " << jz << ')';
  }

protected:
  void do_fit(const boost::multi_array<double, 4>& mat);
};

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

inline std::ostream& operator<<(std::ostream& os, looper::bond_parameter_xxz const& p) {
  p.print(os);
  return os;
}

inline std::ostream& operator<<(std::ostream& os, looper::bond_parameter_xyz const& p) {
  p.print(os);
  return os;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

namespace looper {

template<typename G, typename I>
site_parameter get_site_parameter(alps::Parameters const& param, G const& graph,
  alps::HamiltonianDescriptor<I> const& hd, typename alps::graph_traits<G>::site_descriptor s) {
  int t = get(site_type_t(), graph, s);
  return site_parameter(get_matrix(double(), hd.site_term(t), hd.basis().site_basis(t), param));
}

template<typename G, typename I>
bond_parameter_xxz get_bond_parameter_xxz(alps::Parameters const& param, G const& graph,
  alps::HamiltonianDescriptor<I> const& hd, typename alps::graph_traits<G>::bond_descriptor b) {
  int t = get(bond_type_t(), graph, b);
  int ts = get(site_type_t(), graph, source(b, graph));
  int tt = get(site_type_t(), graph, target(b, graph));
  return bond_parameter_xxz(get_matrix(double(), hd.bond_term(t),
    hd.basis().site_basis(ts), hd.basis().site_basis(tt), param));
}

template<typename G, typename I>
bond_parameter_xyz get_bond_parameter_xyz(alps::Parameters const& param, G const& graph,
  alps::HamiltonianDescriptor<I> const& hd, typename alps::graph_traits<G>::bond_descriptor b) {
  int t = get(bond_type_t(), graph, b);
  int ts = get(site_type_t(), graph, source(b, graph));
  int tt = get(site_type_t(), graph, target(b, graph));
  return bond_parameter_xyz(get_matrix(double(), hd.bond_term(t),
    hd.basis().site_basis(ts), hd.basis().site_basis(tt), param));
}

//
// model_parameter
//

class model_parameter {
public:
  typedef site_parameter::spin_type spin_type;

  model_parameter() :
    sites_(), bonds_(), use_site_indices_(false), use_bond_indices_(false), quantal_(false),
    has_hz_(false), has_d_(false), signed_(false), frustrated_(false),
    site_offset_(0), bond_offset_(0) {
  }
  template<typename G, typename I>
  model_parameter(const alps::Parameters& params, const alps::graph_helper<G>& gh,
    const alps::model_helper<I>& mh) {
    set_parameters(params, gh.graph(), gh.inhomogeneous_sites(), gh.inhomogeneous_bonds(),
      mh.model(), alps::has_sign_problem(mh.model(), gh, params));
  }

  // set_parameters

  template<typename G>
  void set_parameters(G const& g, site_parameter const& sp, bond_parameter_xyz const& bp) {
    set_parameters_impl(g, sp, bp);
    signed_ = check_sign(g);
    frustrated_ = check_frustration(g);
  }
  template<typename G, typename I>
  void set_parameters(const alps::Parameters& params, const G& g, bool inhomogeneous_sites,
    bool inhomogeneous_bond, const alps::HamiltonianDescriptor<I>& hd) {
    set_parameters_impl(params, g, inhomogeneous_sites, inhomogeneous_bond, hd);
    signed_ = check_sign(g);
    frustrated_ = check_frustration(g);
  }
  template<typename G, typename I>
  void set_parameters(const alps::Parameters& params, const G& g, bool inhomogeneous_sites,
    bool inhomogeneous_bond, const alps::HamiltonianDescriptor<I>& hd, bool is_signed)
 {
    set_parameters_impl(params, g, inhomogeneous_sites, inhomogeneous_bond, hd);
    signed_ = is_signed;
    frustrated_ = check_frustration(g);
  }

  // bool uniform_site() const { return sites_.size() == 1; }
  bool inhomogeneous_site() const {
    return use_site_indices_;
  }
  template<class G>
  const site_parameter& site(const typename alps::graph_traits<G>::site_descriptor& v,
    const G& g) const {
    return inhomogeneous_site() ? sites_[v] : sites_map_.find(get(site_type_t(), g, v))->second;
  }

  bool inhomogeneous_bond() const {
    return use_bond_indices_;
  }
  template<class G>
  const bond_parameter_xyz& bond(const typename alps::graph_traits<G>::bond_descriptor& e,
    const G& g) const {
    return inhomogeneous_bond() ? bonds_[get(bond_index_t(), g, e)] :
      bonds_map_.find(boost::make_tuple(get(bond_type_t(), g, e),
        get(site_type_t(), g, source(e, g)), get(site_type_t(), g, target(e, g))))->second;
  }

  bool is_quantal() const { return quantal_; }
  bool has_field() const { return has_hz_; }
  bool has_d_term() const { return has_d_; }
  bool is_signed() const { return signed_; }
  bool is_frustrated() const { return frustrated_; }
  double site_energy_offset() const { return site_offset_; }
  double bond_energy_offset() const { return bond_offset_; }

  template<typename G>
  void output(std::ostream& os, G const& g) const;

protected:
  template<typename G>
  void set_parameters_impl(G const& g, site_parameter const& sp, bond_parameter_xyz const& bp);

  template<typename G, typename I>
  void set_parameters_impl(alps::Parameters params, const G& g, bool inhomogeneous_sites,
    bool inhomogeneous_bond, const alps::HamiltonianDescriptor<I>& hd);

  template<typename G>
  bool check_sign(const G& g) const;

  template<typename G>
  bool check_frustration(const G& g) const;

private:
  std::vector<site_parameter> sites_;
  std::map<int, site_parameter> sites_map_;
  std::vector<bond_parameter_xyz> bonds_;
  std::map<boost::tuple<int, int, int>, bond_parameter_xyz> bonds_map_;
  bool use_site_indices_;
  bool use_bond_indices_;
  bool quantal_;
  bool has_hz_;
  bool has_d_;
  bool signed_;
  bool frustrated_;
  double site_offset_;
  double bond_offset_;
};

//
// implementation of site_parameter
//

inline void site_parameter::do_fit(const boost::multi_array<double, 2>& mat) {
  int dim = mat.shape()[0];
  int m = dim * dim;
  int n = (dim > 2 ? 4 : 3); // D-term is meaningful only for S > 1/2
  s = (double)(dim-1)/2;

  if (!(mat.shape()[0] == mat.shape()[1] && s.get_twice() + 1 == dim))
    boost::throw_exception(std::runtime_error("Error: fitting to site_parameter failed.  "
      "This model is not supported by the current version of looper."));

  site_matrix mat_c(site_parameter(s, 1, 0, 0, 0));
  site_matrix mat_hx(site_parameter(s, 0, 1, 0, 0));
  site_matrix mat_hz(site_parameter(s, 0, 0, 1, 0));
  site_matrix mat_d(site_parameter(s, 0, 0, 0, 1));

  boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> a(m, n);
  boost::numeric::ublas::vector<double> b(m);
  boost::numeric::ublas::vector<double> x(n);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      int k = dim * i + j;
      b(k) = mat[i][j];
      a(k, 0) = mat_c(i,j);
      a(k, 1) = mat_hx(i,j);
      a(k, 2) = mat_hz(i,j);
      if (n == 4) a(k, 3) = mat_d(i,j);
    }
  }

  // call linear least sqaure problem solver
  bool success = alps::numeric::is_zero<1>(solve_llsp(a, b, x));
  if (!success)
    boost::throw_exception(std::runtime_error("Error: fitting to site_parameter failed.  "
      "This model is not supported by the current version of looper."));

  for (int i = 0; i < n; ++i) x(i) = alps::numeric::round<1>(x(i));
  c = x(0);
  hx = x(1);
  hz = x(2);
  d = (n == 4 ? x(3) : 0);
}

//
// implementation of bond_paraemter
//

inline void bond_parameter_xxz::do_fit(const boost::multi_array<double, 4>& mat) {
  int d0 = mat.shape()[0];
  int d1 = mat.shape()[1];
  int dim = d0 * d1;
  int m = dim * dim;
  int n = 3;
  alps::half_integer<int> s0((double)(d0-1)/2);
  alps::half_integer<int> s1((double)(d1-1)/2);

  if (!(mat.shape()[0] == mat.shape()[2] && mat.shape()[1] == mat.shape()[3] &&
        s0.get_twice() + 1 == d0 && s1.get_twice() + 1 == d1))
    boost::throw_exception(std::runtime_error("Error: fitting to bond_parameter_xxz failed.  "
      "This model is not supported by the current version of looper."));

  bond_matrix_xxz mat_c(s0, s1, bond_parameter_xxz(1, 0, 0));
  bond_matrix_xxz mat_jxy(s0, s1, bond_parameter_xxz(0, 1, 0));
  bond_matrix_xxz mat_jz(s0, s1, bond_parameter_xxz(0, 0, 1));

  boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> a(m, n);
  boost::numeric::ublas::vector<double> b(m);
  boost::numeric::ublas::vector<double> x(n);
  for (int i0 = 0; i0 < d0; ++i0)
    for (int i1 = 0; i1 < d1; ++i1)
      for (int j0 = 0; j0 < d0; ++j0)
        for (int j1 = 0; j1 < d1; ++j1) {
          int k = dim * (i0 * d1 + i1) + (j0 * d1 + j1);
          b(k) = mat[i0][i1][j0][j1];
          a(k, 0) = mat_c(i0,i1,j0,j1);
          a(k, 1) = mat_jxy(i0,i1,j0,j1);
          a(k, 2) = mat_jz(i0,i1,j0,j1);
        }

  // call linear least squaure problem solver
  bool success = alps::numeric::is_zero<1>(solve_llsp(a, b, x));
  if (!success)
    boost::throw_exception(std::runtime_error("Error: fitting to bond_parameter_xxz failed.  "
      "This model is not supported by the current version of looper."));

  for (int i = 0; i < n; ++i) x(i) = alps::numeric::round<1>(x(i));
  c = x(0);
  jxy = x(1);
  jz = x(2);
}

inline void bond_parameter_xyz::do_fit(const boost::multi_array<double, 4>& mat) {
  int d0 = mat.shape()[0];
  int d1 = mat.shape()[1];
  int dim = d0 * d1;
  int m = dim * dim;
  int n = 4;
  alps::half_integer<int> s0((double)(d0-1)/2);
  alps::half_integer<int> s1((double)(d1-1)/2);

  if (!(mat.shape()[0] == mat.shape()[2] && mat.shape()[1] == mat.shape()[3] &&
        s0.get_twice() + 1 == d0 && s1.get_twice() + 1 == d1))
    boost::throw_exception(std::runtime_error("Error: fitting to bond_parameter_xyz failed.  "
      "This model is not supported by the current version of looper."));

  bond_matrix_xyz mat_c(s0, s1, bond_parameter_xyz(1, 0, 0, 0));
  bond_matrix_xyz mat_jx(s0, s1, bond_parameter_xyz(0, 1, 0, 0));
  bond_matrix_xyz mat_jy(s0, s1, bond_parameter_xyz(0, 0, 1, 0));
  bond_matrix_xyz mat_jz(s0, s1, bond_parameter_xyz(0, 0, 0, 1));

  boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> a(m, n);
  boost::numeric::ublas::vector<double> b(m);
  boost::numeric::ublas::vector<double> x(n);
  for (int i0 = 0; i0 < d0; ++i0)
    for (int i1 = 0; i1 < d1; ++i1)
      for (int j0 = 0; j0 < d0; ++j0)
        for (int j1 = 0; j1 < d1; ++j1) {
          int k = dim * (i0 * d1 + i1) + (j0 * d1 + j1);
          b(k) = mat[i0][i1][j0][j1];
          a(k, 0) = mat_c(i0,i1,j0,j1);
          a(k, 1) = mat_jx(i0,i1,j0,j1);
          a(k, 2) = mat_jy(i0,i1,j0,j1);
          a(k, 3) = mat_jz(i0,i1,j0,j1);
        }

  // call linear least squaure problem solver
  bool success = alps::numeric::is_zero<1>(solve_llsp(a, b, x));
  if (!success)
    boost::throw_exception(std::runtime_error("Error: fitting to bond_parameter_xyz failed.  "
      "This model is not supported by the current version of looper."));

  for (int i = 0; i < n; ++i) x(i) = alps::numeric::round<1>(x(i));
  c = x(0);
  jx = x(1);
  jy = x(2);
  jz = x(3);
}

//
// implementation of model_parameter
//

template<typename G>
void model_parameter::set_parameters_impl(G const& g, site_parameter const& sp,
  bond_parameter_xyz const& bp) {

  use_site_indices_ = false;
  sites_.clear();
  sites_map_.clear();
  BOOST_FOREACH(typename alps::graph_traits<G>::site_descriptor s, sites(g)) {
    int t = get(site_type_t(), g, s);
    if (sites_map_.find(t) == sites_map_.end())
      sites_map_[t] = sp;
  }

  use_bond_indices_ = false;
  bonds_.clear();
  bonds_map_.clear();
  BOOST_FOREACH(typename alps::graph_traits<G>::bond_descriptor b, bonds(g)) {
    boost::tuple<int, int, int>
      t(get(bond_type_t(), g, b), get(site_type_t(), g, source(b, g)),
      get(site_type_t(), g, target(b, g)));
    if (bonds_map_.find(t) == bonds_map_.end())
      bonds_map_[t] = bp;
  }

  quantal_ = bp.is_quantal();
  has_hz_ = false;
  has_d_ = false;
  signed_ = false;
  frustrated_ = false;
  site_offset_ = 0;
  bond_offset_ = 0;
}

template<typename G, typename I>
void model_parameter::set_parameters_impl(alps::Parameters params, const G& g,
  bool inhomogeneous_sites, bool inhomogeneous_bonds, const alps::HamiltonianDescriptor<I>& hd) {

  quantal_ = false;
  has_hz_ = false;
  has_d_ = false;
  signed_ = false;
  frustrated_ = false;

  alps::basis_states_descriptor<I> basis(hd.basis(), g);
  alps::Disorder::seed(params.value_or_default("DISORDER_SEED",0));

  if (inhomogeneous_sites || inhomogeneous_bonds)
    alps::throw_if_xyz_defined(params, g);

  // site terms

  sites_.clear();
  sites_map_.clear();
  if (inhomogeneous_sites || params.value_or_default("USE_SITE_INDICES_AS_TYPES", false)) {
    use_site_indices_ = true;
    sites_.resize(num_sites(g));
  } else {
    use_site_indices_ = false;
  }

  // generate site matrices and set site parameters
  if (use_site_indices_) {
    alps::Parameters p(params);
    BOOST_FOREACH(typename alps::graph_traits<G>::site_descriptor s, sites(g)) {
      if (inhomogeneous_sites)
        p << alps::coordinate_as_parameter(g, s);
      sites_[s] = get_site_parameter(p, g, hd, s);
      quantal_ |= sites_[s].is_quantal();
      has_hz_ |= sites_[s].has_hz();
      has_d_ |= sites_[s].has_d();
    }
  } else {
    BOOST_FOREACH(typename alps::graph_traits<G>::site_descriptor s, sites(g)) {
      int t = get(site_type_t(), g, s);
      if (sites_map_.find(t) == sites_map_.end()) {
        sites_map_[t] = get_site_parameter(params, g, hd, s);
        quantal_ |= sites_map_[t].is_quantal();
        has_hz_ |= sites_map_[t].has_hz();
        has_d_ |= sites_map_[t].has_d();
      }
    }
  }

  // bond terms

  bonds_.clear();
  bonds_map_.clear();
  if (inhomogeneous_bonds || params.value_or_default("USE_BOND_INDICES_AS_TYPES", false)) {
    use_bond_indices_ = true;
    bonds_.resize(num_bonds(g));
  } else {
    use_bond_indices_ = false;
  }

  // generate bond matrices and set bond parameters
  if (use_bond_indices_) {
    alps::Parameters p(params);
    BOOST_FOREACH(typename alps::graph_traits<G>::bond_descriptor b, bonds(g)) {
      if (inhomogeneous_bonds)
        p << alps::coordinate_as_parameter(g, b);
      int i = get(bond_index_t(), g, b);
      bonds_[i] = get_bond_parameter_xyz(p, g, hd, b);
      quantal_ |= bonds_[i].is_quantal();
    }
  } else {
    BOOST_FOREACH(typename alps::graph_traits<G>::bond_descriptor b, bonds(g)) {
      boost::tuple<int, int, int>
        t(get(bond_type_t(), g, b), get(site_type_t(), g, source(b, g)),
        get(site_type_t(), g, target(b, g)));
      if (bonds_map_.find(t) == bonds_map_.end()) {
        bonds_map_[t] = get_bond_parameter_xyz(params, g, hd, b);
        quantal_ |= bonds_map_[t].is_quantal();
      }
    }
  }

  site_offset_ = 0;
  BOOST_FOREACH(typename alps::graph_traits<G>::site_descriptor s, sites(g))
    site_offset_ += site(s, g).c;
  bond_offset_ = 0;
  BOOST_FOREACH(typename alps::graph_traits<G>::bond_descriptor b, bonds(g))
    bond_offset_ += bond(b, g).c;
  if (has_d_term()) {
    BOOST_FOREACH(typename alps::graph_traits<G>::site_descriptor s, sites(g))
      site_offset_ += 0.5 * site(s, g).s.get_twice() * site(s, g).d;
  }

  if (params.defined("LOOPER_DEBUG[MODEL OUTPUT]")) {
    if (params["LOOPER_DEBUG[MODEL OUTPUT]"] == "cerr")
      output(std::cerr, g);
    else
      output(std::cout, g);
  }

}

template<typename G>
void model_parameter::output(std::ostream& os, G const& g) const {
  os << "site terms: number of sites = " << num_sites(g) << std::endl;
  BOOST_FOREACH(typename alps::graph_helper<G>::site_descriptor s, sites(g))
    os << "  site " << s << ": (S, C, Hx, Hz, D) = " << site(s, g) << std::endl;
  os << "bond terms: number of bonds = " << num_bonds(g) << std::endl;
  BOOST_FOREACH(typename alps::graph_helper<G>::bond_descriptor b, bonds(g))
    os << "  bond " << b << ": (C, Jx, Jy, Jz) = " << bond(b, g) << std::endl;
}

template<typename G>
bool model_parameter::check_sign(const G& g) const {
  std::vector<double> wb, ws;
  BOOST_FOREACH(typename alps::graph_traits<G>::bond_descriptor b, bonds(g)) {
    double w = 1;
    if (alps::numeric::is_zero<1>(bond(b, g).jx + bond(b, g).jy)) {
      if (alps::numeric::is_zero<1>(bond(b, g).jx - bond(b, g).jy))
        w = 0;
      else
        w = (bond(b, g).jx - bond(b, g).jy < 0 ? 1 : -1);
    } else if (bond(b, g).jx + bond(b, g).jy < 0) {
      if (bond(b, g).jx - bond(b, g).jy > 0)
        return true;
      else
        w = 1;
    } else {
      if (bond(b, g).jx - bond(b, g).jy < 0)
        return true;
      else
        w = -1;
    }
    wb.push_back(w);
  }
  BOOST_FOREACH(typename alps::graph_traits<G>::site_descriptor s, sites(g))
    ws.push_back(alps::numeric::is_zero<1>(site(s, g).hx) ? 0 : (site(s, g).hx > 0 ? 1 : -1));
  return alps::is_frustrated(g,
    boost::make_iterator_property_map(wb.begin(), get(bond_index_t(), g)),
    boost::make_iterator_property_map(ws.begin(), get(site_index_t(), g)));
}

template<typename G>
bool model_parameter::check_frustration(const G& g) const {
  std::vector<double> w;
  BOOST_FOREACH(typename alps::graph_traits<G>::bond_descriptor b, bonds(g))
    w.push_back(alps::numeric::is_zero<1>(bond(b, g).jz) ? 0 : (bond(b, g).jz < 0 ? 1 : -1));
  bool frustrated =
    alps::is_frustrated(g, boost::make_iterator_property_map(w.begin(), get(bond_index_t(), g)));
  if (has_d_ && !frustrated) {
    BOOST_FOREACH(typename alps::graph_traits<G>::site_descriptor s, sites(g))
      if (site(s, g).s.get_twice() > 1 && alps::numeric::is_positive<1>(site(s, g).d)) {
        frustrated = true;
        break;
      }
  }
  return frustrated;
}

} // end namespace looper

#endif // LOOPER_MODEL_PARAMETER_H
