/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/


// alps/model/fcidump.hpp: write an ALPS fermion model as an FCIDUMP file.
//
// FCIDUMP (Knowles and Handy, Comput. Phys. Commun. 54, 75 (1989)) is the
// integral format read by quantum chemistry ED, CI and DMRG codes (PySCF,
// Molpro, Block2, Dice, NECI, ...). It holds a spin-independent Hamiltonian
//
//   H = E_core + sum_{pq,a} h_pq c+_{p,a} c_{q,a}
//       + 1/2 sum_{pqrs,ab} (pq|rs) c+_{p,a} c+_{r,b} c_{s,b} c_{q,a}
//
// with spins a, b and one spatial orbital per ALPS site.
//
// The integrals are not derived from the model's term strings. They are read
// off the matrix ALPS itself builds (hamiltonian_matrix, with site, bond and
// n-site terms), in the sectors with at most one up and one down electron:
//
//   (Nup,Ndown) = (0,0):  E_core           = <0|H|0>
//                 (1,0):  h_pq             = <p up|H|q up> - E_core delta_pq
//                 (1,1):  (pq|rs)          = <p up, r down|H|q up, s down> - ...
//
// so the FCIDUMP follows ALPS's operator definitions and fermion signs by
// construction. An ALPS model need not be of this form (spin-dependent
// couplings, three-body n-site terms, ...). export_fcidump then rebuilds every
// sector with up to `check` electrons from the integrals and compares it with
// ALPS's matrix, element by element; any difference is an error and no file
// is written.

#ifndef ALPS_MODEL_FCIDUMP_HPP
#define ALPS_MODEL_FCIDUMP_HPP

#include <alps/lattice.h>
#include <alps/model.h>
#include <alps/model/hamiltonian_matrix.hpp>
#include <alps/model/model_helper.h>
#include <alps/parameter.h>
#include <boost/filesystem/operations.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace alps {
namespace fcidump {


typedef boost::numeric::ublas::mapped_vector_of_mapped_vector<
    double, boost::numeric::ublas::row_major> sparse_matrix;
typedef alps::hamiltonian_matrix<sparse_matrix> hamiltonian;

// A Slater determinant: occupied spin orbitals (modes) in ascending order.
// Mode 2*site is up, 2*site+1 is down: site-major, up before down, the order
// in which ALPS counts fermion signs.
typedef std::vector<int> det;

// Apply ops (rightmost first) to d. Each op is (creation?, mode). Returns
// false if the result is zero; otherwise d is the result and sign its sign
// relative to the ordered product of creators.
inline bool apply(const std::vector<std::pair<bool, int> >& ops, det& d, int& sign)
{
  for (std::size_t k = ops.size(); k-- > 0;) {
    bool create = ops[k].first;
    int m = ops[k].second;
    det::iterator it = std::lower_bound(d.begin(), d.end(), m);
    bool occupied = it != d.end() && *it == m;
    if (occupied == create)
      return false;
    if ((it - d.begin()) % 2)
      sign = -sign;
    if (create)
      d.insert(it, m);
    else
      d.erase(it);
  }
  return true;
}

inline std::string not_fixed()
{
  return "the model's BASIS does not fix Nup and Ndown from Nup_total and "
         "Ndown_total (it needs <CONSTRAINT quantumnumber=\"Nup\" "
         "value=\"Nup_total\"/> and the same for Ndown)";
}

inline std::string not_fermion(std::size_t site)
{
  return "site " + boost::lexical_cast<std::string>(site) +
         " does not use a spin-1/2 fermion basis with quantum numbers "
         "Nup and Ndown (e.g. BASIS \"fermion\"); FCIDUMP needs one "
         "spatial orbital per site";
}

struct sector {
  std::vector<det> dets;               // ALPS basis state i as a determinant
  std::map<det, std::size_t> index;
  std::vector<std::map<std::size_t, double> > rows;  // nonzeros of ALPS's H
};

class extractor {
public:
  explicit extractor(const alps::Parameters& p) : parms_(p), n_(0) {}

  // Refuse, before any basis is enumerated, a site basis without Nup and
  // Ndown, and a BASIS that does not constrain them (every state would be
  // built). Its other constraints would leave sectors incomplete: a
  // constraint of quantum number Q to the parameter Q_total (N_total and
  // Sz_total of "alternative fermion") is lifted by unsetting Q_total, and the
  // values given for N and Sz are kept in totals(). Any other constraint is
  // refused.
  void check_constraints()
  {
    alps::Parameters p(sector_parameters(0, 0));
    alps::model_helper<> m(p);
    alps::model_helper<>::basis_descriptor_type b = m.basis();
    alps::graph_helper<> g(p);
    std::size_t site = 0;
    for (alps::graph_helper<>::site_iterator it = g.sites().first; it != g.sites().second; ++it, ++site) {
      const alps::SiteBasisDescriptor<short>& sb = b.site_basis(g.site_type(*it));
      if (alps::get_quantumnumber_index("Nup", sb) == sb.size() ||
          alps::get_quantumnumber_index("Ndown", sb) == sb.size())
        throw std::runtime_error(not_fermion(site));
    }
    for (std::size_t k = 0; k < b.all_constraints().size(); ++k) {
      const std::string& q = b.all_constraints()[k].first;
      std::string value = boost::lexical_cast<std::string>(b.all_constraints()[k].second);
      if (q == "Nup" || q == "Ndown" || value != q + "_total" || !parms_.defined(value))
        continue;
      if (q == "N" || q == "Sz")
        totals_[q] = alps::evaluate<double>(static_cast<std::string>(parms_[value]),
                                            alps::expression::ParameterEvaluator<double>(parms_));
      parms_.erase(value);
    }
    p = sector_parameters(0, 0);
    b.set_parameters(p);
    bool up = false, down = false;
    for (std::size_t k = 0; k < b.constraints().size(); ++k) {
      const std::string& q = b.constraints()[k].first;
      if (q == "Nup")
        up = true;
      else if (q == "Ndown")
        down = true;
      else {
        std::string value;
        for (std::size_t j = 0; j < b.all_constraints().size(); ++j)
          if (b.all_constraints()[j].first == q)
            value = boost::lexical_cast<std::string>(b.all_constraints()[j].second);
        throw std::runtime_error(
            "the model's BASIS also constrains " + q + " (<CONSTRAINT quantumnumber=\"" +
            q + "\" value=\"" + value + "\"/>) to a fixed value; FCIDUMP export "
            "needs every (Nup,Ndown) sector");
      }
    }
    if (!up || !down)
      throw std::runtime_error(not_fixed());
  }

  // ALPS's Hamiltonian in sector (nup, ndown). The sectors with at most two
  // electrons, used by both extract and verify, are built once; a larger one
  // is valid until the next call.
  const sector& build(int nup, int ndown)
  {
    if (nup + ndown > 2) {
      last_ = sector();
      last_ = make(nup, ndown);
      return last_;
    }
    std::pair<int, int> key(nup, ndown);
    std::map<std::pair<int, int>, sector>::iterator it = small_.find(key);
    if (it == small_.end())
      it = small_.insert(std::make_pair(key, make(nup, ndown))).first;
    return it->second;
  }

  std::size_t sites() const { return n_; }

  // Values the parameters gave N and Sz before check_constraints unset them.
  const std::map<std::string, double>& totals() const { return totals_; }

private:
  // ALPS's Hamiltonian in sector (nup, ndown), basis states as determinants.
  sector make(int nup, int ndown)
  {
    hamiltonian h(sector_parameters(nup, ndown));
    const hamiltonian::basis_states_type& states = h.states_vector();
    const alps::basis_states_descriptor<short>& basis = states.basis();
    n_ = basis.size();

    sector out;
    for (std::size_t i = 0; i < states.size(); ++i) {
      det d;
      for (std::size_t s = 0; s < n_; ++s) {
        const alps::SiteBasisDescriptor<short>& sb = basis.get_site_basis(s);
        std::size_t iu = alps::get_quantumnumber_index("Nup", sb);
        std::size_t id = alps::get_quantumnumber_index("Ndown", sb);
        if (iu == sb.size() || id == sb.size() || basis[s].size() != 4)
          throw std::runtime_error(not_fermion(s));
        int u = alps::to_double(alps::get_quantumnumber(basis[s][states[i][s]], iu));
        int dn = alps::to_double(alps::get_quantumnumber(basis[s][states[i][s]], id));
        if (u) d.push_back(2 * s);
        if (dn) d.push_back(2 * s + 1);
      }
      int nu = 0, nd = 0;
      for (std::size_t k = 0; k < d.size(); ++k) (d[k] % 2 ? nd : nu)++;
      if (nu != nup || nd != ndown)
        throw std::runtime_error(not_fixed());
      out.index[d] = out.dets.size();
      out.dets.push_back(d);
    }
    double expected = choose(n_, nup) * choose(n_, ndown);
    if (out.dets.size() != expected) {
      std::ostringstream o;
      o << "sector (Nup,Ndown) = (" << nup << "," << ndown << ") has " << out.dets.size()
        << " states instead of " << expected << "; the model's BASIS restricts it "
        << "beyond Nup_total and Ndown_total";
      throw std::runtime_error(o.str());
    }

    const sparse_matrix& m = h.matrix();
    out.rows.resize(out.dets.size());
    for (sparse_matrix::const_iterator1 r = m.begin1(); r != m.end1(); ++r)
      for (sparse_matrix::const_iterator2 c = r.begin(); c != r.end(); ++c)
        if (*c != 0.)
          out.rows[c.index1()][c.index2()] = *c;
    return out;
  }

  alps::Parameters sector_parameters(int nup, int ndown) const
  {
    alps::Parameters p(parms_);
    p["Nup_total"] = nup;
    p["Ndown_total"] = ndown;
    p["CONSERVED_QUANTUMNUMBERS"] = "Nup,Ndown";
    return p;
  }

  static double choose(std::size_t n, int k)
  {
    double c = 1.;
    for (int j = 0; j < k; ++j) c = c * (n - j) / (j + 1);
    return c;
  }

  alps::Parameters parms_;
  std::size_t n_;
  std::map<std::string, double> totals_;
  std::map<std::pair<int, int>, sector> small_;
  sector last_;
};

inline double element(const sector& s, std::size_t i, std::size_t j)
{
  std::map<std::size_t, double>::const_iterator it = s.rows[i].find(j);
  return it == s.rows[i].end() ? 0. : it->second;
}

struct integrals {
  std::size_t n;
  double ecore;
  std::vector<double> h;   // n*n
  std::vector<double> g;   // n^4, chemists' (pq|rs) at ((p*n+q)*n+r)*n+s
  double& H(std::size_t p, std::size_t q) { return h[p * n + q]; }
  double H(std::size_t p, std::size_t q) const { return h[p * n + q]; }
  double& G(std::size_t p, std::size_t q, std::size_t r, std::size_t s)
  { return g[((p * n + q) * n + r) * n + s]; }
  double G(std::size_t p, std::size_t q, std::size_t r, std::size_t s) const
  { return g[((p * n + q) * n + r) * n + s]; }
};

inline std::string orbitals(int p, int q) { std::ostringstream o; o << "(" << p + 1 << "," << q + 1 << ")"; return o.str(); }

inline std::string integral(std::size_t p, std::size_t q, std::size_t r, std::size_t s)
{
  std::ostringstream o;
  o << "(" << p + 1 << " " << q + 1 << "|" << r + 1 << " " << s + 1 << ")";
  return o.str();
}

inline void require(double a, double b, double tol, const std::string& what)
{
  if (std::abs(a - b) > tol) {
    std::ostringstream o;
    o << what << ": " << std::setprecision(12) << a << " vs " << b;
    throw std::runtime_error(o.str());
  }
}

inline integrals extract(extractor& x, double tol)
{
  x.check_constraints();
  integrals I;
  const sector& s00 = x.build(0, 0);
  I.n = x.sites();
  const std::size_t n = I.n;
  I.ecore = element(s00, 0, 0);
  I.h.assign(n * n, 0.);
  I.g.assign(n * n * n * n, 0.);

  const sector& s10 = x.build(1, 0);
  const sector& s01 = x.build(0, 1);
  for (std::size_t p = 0; p < n; ++p)
    for (std::size_t q = 0; q < n; ++q) {
      det dp(1, 2 * p), dq(1, 2 * q);
      double up = element(s10, s10.index.at(dp), s10.index.at(dq)) - (p == q ? I.ecore : 0.);
      dp[0] += 1; dq[0] += 1;
      double down = element(s01, s01.index.at(dp), s01.index.at(dq)) - (p == q ? I.ecore : 0.);
      require(up, down, tol, "h" + orbitals(p, q) + " differs between up and down "
              "electrons; FCIDUMP (restricted) needs a spin-independent Hamiltonian");
      I.H(p, q) = up;
    }
  for (std::size_t p = 0; p < n; ++p)
    for (std::size_t q = 0; q < p; ++q)
      require(I.H(p, q), I.H(q, p), tol, "h" + orbitals(p, q) + " is not symmetric");

  // <p up, r down| H |q up, s down> with |p up, r down> = c+_{p up} c+_{r down}|0>
  const sector& s11 = x.build(1, 1);
  std::vector<std::vector<int> > sgn(n, std::vector<int>(n));
  std::vector<std::vector<std::size_t> > idx(n, std::vector<std::size_t>(n));
  for (std::size_t p = 0; p < n; ++p)
    for (std::size_t r = 0; r < n; ++r) {
      std::vector<std::pair<bool, int> > ops;
      ops.push_back(std::make_pair(true, int(2 * p)));
      ops.push_back(std::make_pair(true, int(2 * r + 1)));
      det d;
      int sign = 1;
      apply(ops, d, sign);
      sgn[p][r] = sign;
      idx[p][r] = s11.index.at(d);
    }
  for (std::size_t p = 0; p < n; ++p)
    for (std::size_t q = 0; q < n; ++q)
      for (std::size_t r = 0; r < n; ++r)
        for (std::size_t s = 0; s < n; ++s) {
          double v = sgn[p][r] * sgn[q][s] * element(s11, idx[p][r], idx[q][s]);
          if (p == q && r == s) v -= I.ecore;
          if (r == s) v -= I.H(p, q);
          if (p == q) v -= I.H(r, s);
          I.G(p, q, r, s) = v;
        }

  // FCIDUMP stores (pq|rs) with the 8-fold symmetry of real orbitals.
  for (std::size_t p = 0; p < n; ++p)
    for (std::size_t q = 0; q < n; ++q)
      for (std::size_t r = 0; r < n; ++r)
        for (std::size_t s = 0; s < n; ++s) {
          double v = I.G(p, q, r, s);
          if (std::abs(v - I.G(r, s, p, q)) > tol)
            require(v, I.G(r, s, p, q), tol, integral(p, q, r, s) + " changes when up and "
                    "down are swapped; FCIDUMP needs a spin-independent interaction");
          if (std::abs(v - I.G(q, p, r, s)) > tol || std::abs(v - I.G(p, q, s, r)) > tol)
            require(v, I.G(q, p, r, s) != v ? I.G(q, p, r, s) : I.G(p, q, s, r), tol,
                    integral(p, q, r, s) + " lacks the symmetry (pq|rs) = (qp|rs) = (pq|sr) of real "
                    "orbitals that FCIDUMP assumes");
        }
  return I;
}

// The nonzero integrals, listed by the orbitals they annihilate, so that
// apply_integrals visits only the terms that act on a determinant.
struct terms {
  std::size_t n;
  double ecore;
  std::vector<std::vector<std::pair<int, double> > > one;   // q: (p, h_pq)
  // q*n+s: (p, r, (pq|rs)/2)
  std::vector<std::vector<std::pair<std::pair<int, int>, double> > > two;

  explicit terms(const integrals& I) : n(I.n), ecore(I.ecore), one(n), two(n * n)
  {
    for (std::size_t q = 0; q < n; ++q)
      for (std::size_t p = 0; p < n; ++p)
        if (I.H(p, q) != 0.)
          one[q].push_back(std::make_pair(int(p), I.H(p, q)));
    for (std::size_t p = 0; p < n; ++p)
      for (std::size_t q = 0; q < n; ++q)
        for (std::size_t r = 0; r < n; ++r)
          for (std::size_t s = 0; s < n; ++s)
            if (I.G(p, q, r, s) != 0.)
              two[q * n + s].push_back(std::make_pair(std::make_pair(int(p), int(r)),
                                                      0.5 * I.G(p, q, r, s)));
  }
};

// H|d> from the integrals, as {basis index: amplitude} in sector s.
inline void apply_integrals(const terms& T, const det& d, const sector& s,
                     std::map<std::size_t, double>& out)
{
  if (T.ecore != 0.) out[s.index.at(d)] += T.ecore;
  std::vector<std::pair<bool, int> > ops(2);
  for (std::size_t a = 0; a < d.size(); ++a) {
    int q = d[a] / 2, sp = d[a] % 2;
    for (std::size_t k = 0; k < T.one[q].size(); ++k) {
      ops[0] = std::make_pair(true, 2 * T.one[q][k].first + sp);
      ops[1] = std::make_pair(false, d[a]);
      det e = d; int sign = 1;
      if (apply(ops, e, sign)) out[s.index.at(e)] += sign * T.one[q][k].second;
    }
  }
  // 1/2 (pq|rs) c+_{p,sigma} c+_{r,tau} c_{s,tau} c_{q,sigma}
  ops.resize(4);
  for (std::size_t a = 0; a < d.size(); ++a)
    for (std::size_t b = 0; b < d.size(); ++b) {
      if (a == b) continue;
      int q = d[a] / 2, sig = d[a] % 2, s_ = d[b] / 2, tau = d[b] % 2;
      const std::vector<std::pair<std::pair<int, int>, double> >& t = T.two[q * T.n + s_];
      for (std::size_t k = 0; k < t.size(); ++k) {
        ops[0] = std::make_pair(true, 2 * t[k].first.first + sig);
        ops[1] = std::make_pair(true, 2 * t[k].first.second + tau);
        ops[2] = std::make_pair(false, d[b]);
        ops[3] = std::make_pair(false, d[a]);
        det e = d; int sign = 1;
        if (apply(ops, e, sign)) out[s.index.at(e)] += sign * t[k].second;
      }
    }
}

inline std::string show(const det& d)
{
  std::ostringstream o;
  o << "|";
  for (std::size_t k = 0; k < d.size(); ++k)
    o << (k ? " " : "") << d[k] / 2 + 1 << (d[k] % 2 ? "dn" : "up");
  return o.str() + ">";
}

// Compare ALPS's H with the integrals' H in every sector with at most
// max_electrons electrons. Returns the largest difference.
inline double verify(extractor& x, const integrals& I, int max_electrons, double tol)
{
  double worst = 0.;
  const int n = I.n;
  const terms T(I);
  for (int ne = 0; ne <= max_electrons; ++ne)
    for (int nup = std::max(0, ne - n); nup <= std::min(ne, n); ++nup) {
      const sector& s = x.build(nup, ne - nup);
      for (std::size_t j = 0; j < s.dets.size(); ++j) {
        std::map<std::size_t, double> col;
        apply_integrals(T, s.dets[j], s, col);
        // ALPS stores y[i] += H(i,j) x[j]; compare row j of ALPS's matrix
        // (its H acting from j) with the integrals acting on j. For a real
        // symmetric H the two conventions agree; a non-symmetric ALPS H shows
        // up as a difference here.
        std::map<std::size_t, double> alps = s.rows[j];
        for (std::map<std::size_t, double>::const_iterator it = col.begin(); it != col.end(); ++it)
          alps[it->first] -= it->second;
        for (std::map<std::size_t, double>::const_iterator it = alps.begin(); it != alps.end(); ++it) {
          worst = std::max(worst, std::abs(it->second));
          if (std::abs(it->second) > tol) {
            std::ostringstream o;
            o << "the model is not a spin-independent two-body Hamiltonian: in "
              << "sector (Nup,Ndown) = (" << nup << "," << ne - nup << ") the "
              << "element <" << show(s.dets[j]).substr(1) << "|H" << show(s.dets[it->first])
              << " differs by " << it->second << " from the extracted integrals "
              << "(e.g. a spin-dependent coupling or a three-body term)";
            throw std::runtime_error(o.str());
          }
        }
      }
    }
  return worst;
}

// The integrals go to file.tmp, which replaces file only once all of it is
// written, so a failed write (full disk, size limit) leaves no truncated file.
inline void write(const integrals& I, const std::string& file, int nelec, int ms2, double cutoff)
{
  const std::string tmp = file + ".tmp";
  std::ofstream out(tmp.c_str());
  if (!out)
    throw std::runtime_error("cannot write " + file);
  const std::size_t n = I.n;
  out << " &FCI NORB=" << n << ",NELEC=" << nelec << ",MS2=" << ms2 << ",\n  ORBSYM=";
  for (std::size_t p = 0; p < n; ++p) out << "1,";
  out << "\n  ISYM=1,\n &END\n";
  out << std::scientific << std::setprecision(16);
  // Unique (pq|rs): p >= q, r >= s, pair (p,q) >= pair (r,s).
  for (std::size_t p = 0; p < n; ++p)
    for (std::size_t q = 0; q <= p; ++q)
      for (std::size_t r = 0; r <= p; ++r)
        for (std::size_t s = 0; s <= (r == p ? q : r); ++s) {
          double v = I.G(p, q, r, s);
          if (std::abs(v) > cutoff)
            out << std::setw(24) << v << std::setw(5) << p + 1 << std::setw(5) << q + 1
                << std::setw(5) << r + 1 << std::setw(5) << s + 1 << "\n";
        }
  for (std::size_t p = 0; p < n; ++p)
    for (std::size_t q = 0; q <= p; ++q)
      if (std::abs(I.H(p, q)) > cutoff)
        out << std::setw(24) << I.H(p, q) << std::setw(5) << p + 1 << std::setw(5) << q + 1
            << std::setw(5) << 0 << std::setw(5) << 0 << "\n";
  out << std::setw(24) << I.ecore << std::setw(5) << 0 << std::setw(5) << 0
      << std::setw(5) << 0 << std::setw(5) << 0 << "\n";
  out.close();
  boost::system::error_code ec;
  if (out)
    boost::filesystem::rename(tmp, file, ec);
  if (!out || ec) {
    boost::filesystem::remove(tmp, ec);
    throw std::runtime_error("cannot write " + file);
  }
}

inline int as_int(const alps::Parameters& p, const std::string& key)
{
  double v = alps::evaluate<double>(static_cast<std::string>(p[key]),
                                    alps::expression::ParameterEvaluator<double>(p));
  if (v != std::floor(v))
    throw std::runtime_error(key + " must be an integer");
  return int(v);
}

struct summary {
  std::size_t orbitals;
  int nelec, ms2;
  double ecore;
  int checked;    // sectors with up to this many electrons were compared
  double worst;   // largest |difference| between ALPS's H and the integrals
};

// Write the Hamiltonian described by p (LATTICE/GRAPH, MODEL, couplings,
// libraries) to file. NELEC and MS2 come from Nup_total and Ndown_total, or
// from the values of the BASIS's N and Sz constraints (N_total, Sz_total),
// or half filling if none is set. TOTAL_MOMENTUM is ignored: the file holds
// the real-space Hamiltonian. Throws std::runtime_error, without writing, if
// the model cannot be held by an FCIDUMP.
inline summary export_fcidump(const alps::Parameters& parms, const std::string& file,
                              int check = 4, double tol = 1e-10, double cutoff = 1e-14)
{
  alps::Parameters p(parms);
  if (p.defined("TOTAL_MOMENTUM"))
    p.erase("TOTAL_MOMENTUM");
  check = std::max(check, 2);
  // A diag task may fix one of them and scan the other's sectors.
  const bool up = p.defined("Nup_total"), down = p.defined("Ndown_total");
  if (up != down)
    std::cerr << "Note: only " << (up ? "Nup_total" : "Ndown_total")
              << " is set; NELEC and MS2 do not use it.\n";

  extractor x(p);
  integrals I = extract(x, tol);
  const int n = I.n;

  int nup, ndown;
  const std::map<std::string, double>& totals = x.totals();
  if (up && down) {
    nup = as_int(p, "Nup_total");
    ndown = as_int(p, "Ndown_total");
    if ((totals.count("N") && totals.at("N") != nup + ndown) ||
        (totals.count("Sz") && 2 * totals.at("Sz") != nup - ndown)) {
      std::ostringstream o;
      o << "Nup_total = " << nup << " and Ndown_total = " << ndown << " disagree with ";
      if (totals.count("N")) o << "N_total = " << totals.at("N") << (totals.count("Sz") ? " and " : "");
      if (totals.count("Sz")) o << "Sz_total = " << totals.at("Sz");
      throw std::runtime_error(o.str());
    }
  } else if (!totals.empty()) {
    double ne = totals.count("N") ? totals.at("N") : n;
    double ms2 = totals.count("Sz") ? 2 * totals.at("Sz") : std::fmod(ne, 2.);
    if (ne != std::floor(ne) || ms2 != std::floor(ms2) || std::fmod(ne + ms2, 2.) != 0.) {
      std::ostringstream o;
      o << "N = " << ne << " and Sz = " << ms2 / 2 << " do not give whole numbers of "
        << "up and down electrons";
      throw std::runtime_error(o.str());
    }
    nup = int(ne + ms2) / 2;
    ndown = int(ne - ms2) / 2;
    if (!totals.count("N") || !totals.count("Sz"))
      std::cerr << "Note: " << (totals.count("N") ? "Sz" : "N") << " not set; writing NELEC="
                << nup + ndown << ", MS2=" << nup - ndown << ".\n";
  } else {
    nup = (n + 1) / 2;
    ndown = n / 2;
    if (!up && !down)
      std::cerr << "Note: Nup_total/Ndown_total not set; writing half filling, NELEC="
                << n << ", MS2=" << nup - ndown << ".\n";
  }
  if (nup < 0 || ndown < 0 || nup > n || ndown > n) {
    std::ostringstream o;
    o << nup << " up and " << ndown << " down electrons do not fit in " << n << " orbitals";
    throw std::runtime_error(o.str());
  }

  summary s;
  s.checked = std::min(check, 2 * n);
  s.worst = verify(x, I, s.checked, tol);

  write(I, file, nup + ndown, nup - ndown, cutoff);
  s.orbitals = n;
  s.nelec = nup + ndown;
  s.ms2 = nup - ndown;
  s.ecore = I.ecore;
  return s;
}

}  // namespace fcidump
}  // namespace alps

#endif  // ALPS_MODEL_FCIDUMP_HPP
