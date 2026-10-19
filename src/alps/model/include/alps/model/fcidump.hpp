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
#include <alps/parameter.h>
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

struct sector {
  std::vector<det> dets;               // ALPS basis state i as a determinant
  std::map<det, std::size_t> index;
  std::vector<std::map<std::size_t, double> > rows;  // nonzeros of ALPS's H
};

class extractor {
public:
  explicit extractor(const alps::Parameters& p) : parms_(p), n_(0) {}

  // ALPS's Hamiltonian in sector (nup, ndown), basis states as determinants.
  sector build(int nup, int ndown)
  {
    alps::Parameters p(parms_);
    p["Nup_total"] = nup;
    p["Ndown_total"] = ndown;
    p["CONSERVED_QUANTUMNUMBERS"] = "Nup,Ndown";
    hamiltonian h(p);
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
          throw std::runtime_error(
              "site " + boost::lexical_cast<std::string>(s) +
              " does not use a spin-1/2 fermion basis with quantum numbers "
              "Nup and Ndown (e.g. BASIS \"fermion\"); FCIDUMP needs one "
              "spatial orbital per site");
        int u = alps::to_double(alps::get_quantumnumber(basis[s][states[i][s]], iu));
        int dn = alps::to_double(alps::get_quantumnumber(basis[s][states[i][s]], id));
        if (u) d.push_back(2 * s);
        if (dn) d.push_back(2 * s + 1);
      }
      int nu = 0, nd = 0;
      for (std::size_t k = 0; k < d.size(); ++k) (d[k] % 2 ? nd : nu)++;
      if (nu != nup || nd != ndown)
        throw std::runtime_error(
            "the model's BASIS does not fix Nup and Ndown from Nup_total and "
            "Ndown_total (it needs <CONSTRAINT quantumnumber=\"Nup\" "
            "value=\"Nup_total\"/> and the same for Ndown)");
      out.index[d] = out.dets.size();
      out.dets.push_back(d);
    }

    const sparse_matrix& m = h.matrix();
    out.rows.resize(out.dets.size());
    for (sparse_matrix::const_iterator1 r = m.begin1(); r != m.end1(); ++r)
      for (sparse_matrix::const_iterator2 c = r.begin(); c != r.end(); ++c)
        if (*c != 0.)
          out.rows[c.index1()][c.index2()] = *c;
    return out;
  }

  std::size_t sites() const { return n_; }

private:
  alps::Parameters parms_;
  std::size_t n_;
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
  double& G(std::size_t p, std::size_t q, std::size_t r, std::size_t s)
  { return g[((p * n + q) * n + r) * n + s]; }
};

inline std::string orbitals(int p, int q) { std::ostringstream o; o << "(" << p + 1 << "," << q + 1 << ")"; return o.str(); }

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
  integrals I;
  sector s00 = x.build(0, 0);
  I.n = x.sites();
  const std::size_t n = I.n;
  I.ecore = element(s00, 0, 0);
  I.h.assign(n * n, 0.);
  I.g.assign(n * n * n * n, 0.);

  sector s10 = x.build(1, 0), s01 = x.build(0, 1);
  for (std::size_t p = 0; p < n; ++p)
    for (std::size_t q = 0; q < n; ++q) {
      det dp(1, 2 * p), dq(1, 2 * q);
      double up = element(s10, s10.index[dp], s10.index[dq]) - (p == q ? I.ecore : 0.);
      dp[0] += 1; dq[0] += 1;
      double down = element(s01, s01.index[dp], s01.index[dq]) - (p == q ? I.ecore : 0.);
      require(up, down, tol, "h" + orbitals(p, q) + " differs between up and down "
              "electrons; FCIDUMP (restricted) needs a spin-independent Hamiltonian");
      I.H(p, q) = up;
    }
  for (std::size_t p = 0; p < n; ++p)
    for (std::size_t q = 0; q < p; ++q)
      require(I.H(p, q), I.H(q, p), tol, "h" + orbitals(p, q) + " is not symmetric");

  // <p up, r down| H |q up, s down> with |p up, r down> = c+_{p up} c+_{r down}|0>
  sector s11 = x.build(1, 1);
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
      idx[p][r] = s11.index[d];
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
          std::string e = "(" + boost::lexical_cast<std::string>(p + 1) + " " +
                          boost::lexical_cast<std::string>(q + 1) + "|" +
                          boost::lexical_cast<std::string>(r + 1) + " " +
                          boost::lexical_cast<std::string>(s + 1) + ")";
          double v = I.G(p, q, r, s);
          if (std::abs(v - I.G(r, s, p, q)) > tol)
            require(v, I.G(r, s, p, q), tol, e + " changes when up and down are "
                    "swapped; FCIDUMP needs a spin-independent interaction");
          if (std::abs(v - I.G(q, p, r, s)) > tol || std::abs(v - I.G(p, q, s, r)) > tol)
            require(v, I.G(q, p, r, s) != v ? I.G(q, p, r, s) : I.G(p, q, s, r), tol,
                    e + " lacks the symmetry (pq|rs) = (qp|rs) = (pq|sr) of real "
                    "orbitals that FCIDUMP assumes");
        }
  return I;
}

// H|d> from the integrals, as {basis index: amplitude} in sector s.
inline void apply_integrals(integrals& I, const det& d, const sector& s,
                     std::map<std::size_t, double>& out)
{
  const int n = I.n;
  if (I.ecore != 0.) out[s.index.at(d)] += I.ecore;
  for (std::size_t a = 0; a < d.size(); ++a) {
    int q = d[a] / 2, sp = d[a] % 2;
    for (int p = 0; p < n; ++p) {
      double v = I.H(p, q);
      if (v == 0.) continue;
      std::vector<std::pair<bool, int> > ops;
      ops.push_back(std::make_pair(true, 2 * p + sp));
      ops.push_back(std::make_pair(false, d[a]));
      det e = d; int sign = 1;
      if (apply(ops, e, sign)) out[s.index.at(e)] += sign * v;
    }
  }
  // 1/2 (pq|rs) c+_{p,sigma} c+_{r,tau} c_{s,tau} c_{q,sigma}
  for (std::size_t a = 0; a < d.size(); ++a)
    for (std::size_t b = 0; b < d.size(); ++b) {
      if (a == b) continue;
      int q = d[a] / 2, sig = d[a] % 2, s_ = d[b] / 2, tau = d[b] % 2;
      for (int p = 0; p < n; ++p)
        for (int r = 0; r < n; ++r) {
          double v = 0.5 * I.G(p, q, r, s_);
          if (v == 0.) continue;
          std::vector<std::pair<bool, int> > ops;
          ops.push_back(std::make_pair(true, 2 * p + sig));
          ops.push_back(std::make_pair(true, 2 * r + tau));
          ops.push_back(std::make_pair(false, d[b]));
          ops.push_back(std::make_pair(false, d[a]));
          det e = d; int sign = 1;
          if (apply(ops, e, sign)) out[s.index.at(e)] += sign * v;
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
inline double verify(extractor& x, integrals& I, int max_electrons, double tol)
{
  double worst = 0.;
  const int n = I.n;
  for (int ne = 0; ne <= max_electrons; ++ne)
    for (int nup = std::max(0, ne - n); nup <= std::min(ne, n); ++nup) {
      sector s = x.build(nup, ne - nup);
      for (std::size_t j = 0; j < s.dets.size(); ++j) {
        std::map<std::size_t, double> col;
        apply_integrals(I, s.dets[j], s, col);
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

inline void write(const integrals& I0, const std::string& file, int nelec, int ms2, double cutoff)
{
  integrals I = I0;
  std::ofstream out(file.c_str());
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
// half filling if they are not set. Throws std::runtime_error, without
// writing, if the model cannot be held by an FCIDUMP.
inline summary export_fcidump(const alps::Parameters& p, const std::string& file,
                              int check = 4, double tol = 1e-10, double cutoff = 1e-14)
{
  if (p.defined("TOTAL_MOMENTUM"))
    throw std::runtime_error("TOTAL_MOMENTUM is set; FCIDUMP export needs the real-space basis");
  check = std::max(check, 2);

  extractor x(p);
  integrals I = extract(x, tol);
  const int n = I.n;

  int nup, ndown;
  if (p.defined("Nup_total") && p.defined("Ndown_total")) {
    nup = as_int(p, "Nup_total");
    ndown = as_int(p, "Ndown_total");
  } else {
    nup = (n + 1) / 2;
    ndown = n / 2;
    std::cerr << "Note: Nup_total/Ndown_total not set; writing half filling, NELEC="
              << n << ", MS2=" << nup - ndown << ".\n";
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
