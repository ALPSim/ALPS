// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Written for this tree, not ported: the quantum-chemistry Hamiltonian path
// of the molecular solver (src/molecular.ipp). The method is the normal and
// complementary operator construction of ab initio DMRG (S. R. White and
// R. L. Martin, J. Chem. Phys. 110, 4127 (1999); the operator classes follow
// G. K.-L. Chan and M. Head-Gordon, J. Chem. Phys. 116, 4462 (2002)); the
// code is new and works on dmtk's blocks, bases and products.
//
// Every block keeps, between DMRG steps, its Hamiltonian, the creation and
// annihilation operators of its own spin orbitals, their normal pair
// products and the complementary operators of the spin orbitals outside it
// (pair operators P and Q, three-operator sums S). A step enlarges each block
// by its neighbouring site, builds the enlarged blocks' operators from the
// stored ones and the site's, applies the Hamiltonian as a sum of products of
// one left and one right enlarged-block operator on the condensed
// (B1 B2)(B3 B4) state, and after the truncation rotates the grown block's
// operators into the new basis. The work per step is O(n^2) operators instead
// of the O(n^4) terms the generic path rebuilds from the term list.
//
// Conventions. Spin orbital I = 2 * orbital + spin (spin 0 up, 1 down);
// sites are dmtk sites in chain order, and an operator product written
// A B with A on the earlier and B on the later modes is the second-quantized
// product, the Jordan-Wigner sign of B taken over the earlier modes' ket state
// as in dmtk's products. With t the one-body and w(P,Q,R,S) the coefficient
// of 1/2 a+_P a+_Q a_R a_S, w(P,Q,R,S) = (p s|q r) for spin(P) = spin(S) and
// spin(Q) = spin(R), a block Z with spin orbitals z keeps
//   P[K,L] = sum_{z,z'} w(K,L,z,z') a_z a_z'                 (K < L outside)
//   Q[K,L] = sum_{z,z'} g(K,L,z,z') a+_z a_z',
//            g(A,B,C,D) = w(A,C,D,B) - w(A,C,B,D)             (K <= L outside)
//   S[K]   = 1/2 sum_z t(z,K) a+_z + sum_{P,Q,R} w(P,Q,R,K) a+_P a+_Q a_R
//   A[I,J] = a+_I a+_J (I < J inside), B[I,J] = a+_I a_J (I <= J inside).
// For disjoint mode sets X before Y the terms of H coupling them are
//   sum_{y in Y} (S^X[y] a_y + h.c.) - sum_{x in X} (a_x S^Y[x] + h.c.)
//   + sum_{x<x' in X} (A^X[x,x'] P^Y[x,x'] + h.c.)
//   + sum_{x<=x' in X} (B^X[x,x'] Q^Y[x,x'] + h.c. for x != x'),
// or, with the normal operators on Y, the last two lines with X and Y swapped
// (P^X[y,y'] A^Y[y,y'] and Q^X[y,y'] B^Y[y,y']). The same expression couples
// a block and its new site inside an enlarged block and the two enlarged
// blocks of the superblock; the normal operators sit on the side with fewer
// spin orbitals.
//
// Measurement (measure_products). Products of at most four single-mode
// factors, such as the RDM elements, are read from the final state of
// final_sweep() by a walk of the wavefunction centre from the closing
// configuration to the left end, instead of as dmtk measurement terms (whose
// stored in-block pieces grow as n^4 m^2). At centre site k the state is
// sum C[l,s,r] |l>|s>|r>: l a state of the left block of k sites (the final
// left-to-right half's block, whose a, a+, A and B were filed when it grew),
// s a state of site k, r a state of a right block built here: an
// orthonormal basis of the span of the centre's rows over (site, right
// block), by Householder QR without truncation, so the walk never
// approximates the state. A product is measured at k = its third-lowest site
// (second for two factors), clamped to [1, ls-3]: the left block then holds at
// most two of its factors (a stored normal product; at k = 1 the block is site
// 0 itself), the right block at most one (a single rotated here; at k = ls-3
// the right block is the last two sites, exact). Grouped by part, with the
// permutation sign, a product is the plain tensor product
// (O_L S_L^pl) x (O_D S_D^pr) x O_R on dmtk's product basis, O the part's
// normal product, S a part's parity and pl, pr the parities of the left and
// right factors, and its value is
//   sum_{s,s'} O_D S_D^pr[s,s'] <O_L S_L^pl, C_s O_R C_s'^T>,
// <,> the Frobenius product on the left block.

#ifndef DMTK_MOLECULAR_SYSTEM_H_
#define DMTK_MOLECULAR_SYSTEM_H_

#include "dmtk/system.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cstddef>
#include <deque>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace dmtk
{

// Restricted real integrals: h[p*n+q], eri[((p*n+q)*n+r)*n+s] = (pq|rs).
struct MolecularIntegrals
{
  int n = 0;
  std::vector<double> h;
  std::vector<double> eri;

  double t(int P, int Q) const
    { return (P & 1) == (Q & 1) ? h[size_t(P >> 1) * n + (Q >> 1)] : 0.; }
  double w(int P, int Q, int R, int S) const
    {
      if((P & 1) != (S & 1) || (Q & 1) != (R & 1)) return 0.;
      return eri[((size_t(P >> 1) * n + (S >> 1)) * n + (Q >> 1)) * n + (R >> 1)];
    }
  double g(int A, int B, int C, int D) const
    { return w(A, C, D, B) - w(A, C, B, D); }
};

namespace molecular
{

using Op = BasicOp<double>;
using Sub = SubMatrix<double>;

// A dense operator on one site, <r|O|c> at m[r * 4 + c], in the order of the
// site's charge-sorted basis.
struct Dense
{
  std::array<double, 16> m{};
  double& operator()(int r, int c) { return m[r * 4 + c]; }
  double operator()(int r, int c) const { return m[r * 4 + c]; }
  Dense& add(const Dense& o, double s)
    { for(int k = 0; k < 16; k++) m[k] += s * o.m[k]; return *this; }
};

inline Dense mul(const Dense& a, const Dense& b)
{
  Dense r;
  for(int i = 0; i < 4; i++)
    for(int k = 0; k < 4; k++){
      double x = a(i, k);
      if(x == 0.) continue;
      for(int j = 0; j < 4; j++) r(i, j) += x * b(k, j);
    }
  return r;
}

inline Dense transposed(const Dense& a)
{
  Dense r;
  for(int i = 0; i < 4; i++) for(int j = 0; j < 4; j++) r(i, j) = a(j, i);
  return r;
}

// An operator of charge change dqn with zero blocks on the given subspaces.
// Its fermion flag follows the parity of dqn (BasicOp::fermion()).
inline Op make_op(const PackedBasis& basis, const QN& dqn)
{
  Op op;
  op.dqn = dqn;
  op.set_fermion(true).set_diagonal(false);
  op.resize(basis);
  return op;
}

// A block's elements are its first cols() * rows() entries: dmtk matrices
// keep their storage when assigned a smaller shape.
inline size_t elements(const Sub& b) { return b.cols() * b.rows(); }

inline bool is_zero(const Sub& b)
{
  const double* x = b.array();
  for(size_t k = 0; k < elements(b); k++) if(x[k] != 0.) return false;
  return true;
}

// Drops the blocks that hold only zeros; products skip missing blocks.
inline void prune(Op& op)
{
  if(std::none_of(op.begin(), op.end(), is_zero)) return;
  std::vector<Sub> kept;
  for(const Sub& b : op) if(!is_zero(b)) kept.push_back(b);
  op.clear();
  for(const Sub& b : kept) op.push_back(b);
}

inline void axpy(Op& dst, const Op& src, double s)
{
  if(s == 0.) return;
  for(const Sub& b : src){
    Sub* t = dst.block(b.qn());
    if(!t || t->cols() != b.cols() || t->rows() != b.rows())
      throw std::logic_error("molecular dmtk: operator blocks do not match");
    double* x = t->array();
    const double* y = b.array();
    for(size_t k = 0; k < elements(b); k++) x[k] += s * y[k];
  }
}

inline Op transpose(const Op& src)
{
  Op dst = make_op(src.subspaces(), QN() - src.dqn);
  for(const Sub& b : src){
    Sub* t = dst.block(b.qn() + src.dqn);
    if(!t) throw std::logic_error("molecular dmtk: transposed block missing");
    for(size_t i = 0; i < b.cols(); i++)
      for(size_t j = 0; j < b.rows(); j++) (*t)(j, i) = b(i, j);
  }
  prune(dst);
  return dst;
}

// dst = rho+ src rho, src on the enlarged basis, dst on the kept states.
inline Op rotate(const Op& src, const BMatrix<double>& rho, const PackedBasis& kept)
{
  Op dst = make_op(kept, src.dqn);
  Matrix<double> aux;
  for(const Sub& b : src){
    const Sub* u1 = rho.block(b.qn());
    const Sub* u2 = rho.block(b.qn() + src.dqn);
    if(!u1 || !u2) continue;
    Sub* t = dst.block(b.qn());
    if(!t) continue;
    aux.reshape(u1->cols(), b.rows());
    matrix_matrix_product('N', 'N', b, *u1, aux);
    matrix_matrix_product('C', 'N', *u2, aux, *t, 1., 1.);
  }
  prune(dst);
  return dst;
}

// The product basis of two parts, first part before second in mode order, in
// the reordered layout that dmtk uses for (B1 B2) and (B3 B4) (hint(),
// rotate(), VectorState::condense).
struct Enlarged
{
  Basis basis;
  PackedBasis first, second;
  size_t d2 = 0;
  std::vector<int> pos; // pos[i1 * d2 + i2] = index in basis

  Enlarged() = default;
  Enlarged(const Basis& b1, const Basis& b2): basis(b1, b2)
    {
      basis.reorder();
      first = b1.subspaces();
      second = b2.subspaces();
      d2 = b2.size();
      pos.assign(b1.size() * b2.size(), -1);
      for(size_t k = 0; k < basis.size(); k++){
        State s = basis[k];
        pos[s.i1 * d2 + s.i2] = int(k);
      }
    }
};

struct Entry { int c, r; double v; };

inline void nonzeros(const Sub& b, std::vector<Entry>& out)
{
  out.clear();
  int const c0 = b.col_range().begin(), r0 = b.row_range().begin();
  for(size_t i = 0; i < b.cols(); i++)
    for(size_t j = 0; j < b.rows(); j++){
      double v = b(i, j);
      if(v != 0.) out.push_back({c0 + int(i), r0 + int(j), v});
    }
}

// E += coef (o1 on the first part) (o2 on the second part).
inline void kron(Op& E, const Op& o1, const Op& o2, double coef, const Enlarged& en)
{
  if(coef == 0. || o1.empty() || o2.empty()) return;
  bool const f2 = o2.fermion();
  std::vector<std::vector<Entry> > e2(o2.size());
  for(size_t k = 0; k < o2.size(); k++) nonzeros(*o2[k], e2[k]);
  std::vector<Entry> e1;
  for(const Sub& b1 : o1){
    nonzeros(b1, e1);
    if(e1.empty()) continue;
    double const s = (f2 && b1.qn().fermion_sign() == -1) ? -coef : coef;
    for(size_t k = 0; k < o2.size(); k++){
      if(e2[k].empty()) continue;
      const Sub& b2 = *o2[k];
      Sub* t = E.block(b1.qn() + b2.qn());
      if(!t) throw std::logic_error("molecular dmtk: enlarged block missing");
      int const tc = t->col_range().begin(), tr = t->row_range().begin();
      int const nc = int(t->cols()), nr = int(t->rows());
      for(const Entry& y : e2[k])
        for(const Entry& x : e1){
          int const i = en.pos[x.c * en.d2 + y.c] - tc, j = en.pos[x.r * en.d2 + y.r] - tr;
          if(i < 0 || i >= nc || j < 0 || j >= nr) throw std::logic_error("molecular dmtk: charge change does not match the enlarged operator");
          (*t)(i, j) += s * x.v * y.v;
        }
    }
  }
}

// E += coef (o1 on the first part) x identity.
inline void embed_first(Op& E, const Op& o1, double coef, const Enlarged& en)
{
  if(coef == 0.) return;
  std::vector<Entry> e1;
  for(const Sub& b1 : o1){
    nonzeros(b1, e1);
    if(e1.empty()) continue;
    for(const SubSpace& s2 : en.second){
      Sub* t = E.block(b1.qn() + s2.qn());
      if(!t) throw std::logic_error("molecular dmtk: enlarged block missing");
      int const tc = t->col_range().begin(), tr = t->row_range().begin();
      for(int i2 = s2.begin(); i2 <= s2.end(); i2++)
        for(const Entry& x : e1)
          (*t)(en.pos[x.c * en.d2 + i2] - tc, en.pos[x.r * en.d2 + i2] - tr) += coef * x.v;
    }
  }
}

// E += coef identity x (o2 on the second part).
inline void embed_second(Op& E, const Op& o2, double coef, const Enlarged& en)
{
  if(coef == 0.) return;
  bool const f2 = o2.fermion();
  std::vector<Entry> e2;
  for(const Sub& b2 : o2){
    nonzeros(b2, e2);
    if(e2.empty()) continue;
    for(const SubSpace& s1 : en.first){
      double const s = (f2 && s1.qn().fermion_sign() == -1) ? -coef : coef;
      Sub* t = E.block(s1.qn() + b2.qn());
      if(!t) throw std::logic_error("molecular dmtk: enlarged block missing");
      int const tc = t->col_range().begin(), tr = t->row_range().begin();
      for(int i1 = s1.begin(); i1 <= s1.end(); i1++)
        for(const Entry& y : e2)
          (*t)(en.pos[i1 * en.d2 + y.c] - tc, en.pos[i1 * en.d2 + y.r] - tr) += s * y.v;
    }
  }
}

// The operators a block keeps (see the file comment). Inside spin orbitals
// are [lo, hi); an empty range is [nso, nso). Outside ones are all others.
struct Set
{
  int lo = 0, hi = 0, nso = 0;
  Op H;
  std::vector<Op> a, ad, A, B, P, Q, S;
  std::vector<int> outside, opos;

  static size_t lt(int i, int j) { return size_t(j) * (j - 1) / 2 + i; } // i < j
  static size_t le(int i, int j) { return size_t(j) * (j + 1) / 2 + i; } // i <= j

  void init(int lo_, int hi_, int nso_)
    {
      lo = lo_; hi = hi_; nso = nso_;
      outside.clear();
      opos.assign(nso, -1);
      for(int K = 0; K < nso; K++)
        if(K < lo || K >= hi) { opos[K] = int(outside.size()); outside.push_back(K); }
      int const ni = hi - lo, no = int(outside.size());
      a.assign(ni, Op()); ad.assign(ni, Op());
      A.assign(lt(0, ni), Op()); B.assign(le(0, ni), Op());
      P.assign(lt(0, no), Op()); Q.assign(le(0, no), Op()); S.assign(no, Op());
    }
  bool inside(int I) const { return I >= lo && I < hi; }
  Op& a_(int I) { return a[I - lo]; }
  Op& ad_(int I) { return ad[I - lo]; }
  Op& A_(int I, int J) { return A[lt(I - lo, J - lo)]; }
  Op& B_(int I, int J) { return B[le(I - lo, J - lo)]; }
  Op& P_(int K, int L) { return P[lt(opos[K], opos[L])]; }
  Op& Q_(int K, int L) { return Q[le(opos[K], opos[L])]; }
  Op& S_(int K) { return S[opos[K]]; }
  const Op& a_(int I) const { return a[I - lo]; }
  const Op& ad_(int I) const { return ad[I - lo]; }
  const Op& A_(int I, int J) const { return A[lt(I - lo, J - lo)]; }
  const Op& B_(int I, int J) const { return B[le(I - lo, J - lo)]; }
  const Op& P_(int K, int L) const { return P[lt(opos[K], opos[L])]; }
  const Op& Q_(int K, int L) const { return Q[le(opos[K], opos[L])]; }
  const Op& S_(int K) const { return S[opos[K]]; }

  // Every operator, in a fixed order (file layout and rotation).
  template<class F> void each(F f)
    {
      f(H);
      for(auto* v : {&a, &ad, &A, &B, &P, &Q, &S}) for(Op& o : *v) f(o);
    }

  void write(std::ostream& s)
    {
      s.write(reinterpret_cast<const char*>(&lo), sizeof(int));
      s.write(reinterpret_cast<const char*>(&hi), sizeof(int));
      s.write(reinterpret_cast<const char*>(&nso), sizeof(int));
      each([&](Op& o) { o.write(s); });
    }
  void read(std::istream& s)
    {
      int l, h, n;
      s.read(reinterpret_cast<char*>(&l), sizeof(int));
      s.read(reinterpret_cast<char*>(&h), sizeof(int));
      s.read(reinterpret_cast<char*>(&n), sizeof(int));
      init(l, h, n);
      each([&](Op& o) { o.read(s); });
    }
};

// One factor of a measured product: a+ (dagger) or a of a spin orbital.
struct Factor { int mode; bool dagger; };

// A charge sector of a block basis.
struct Sector { QN qn; int begin, dim; };

inline std::vector<Sector> sectors(const PackedBasis& b)
{
  std::vector<Sector> out;
  for(const SubSpace& s : b) out.push_back({s.qn(), int(s.begin()), int(s.dim())});
  return out;
}

inline int sector_of(const std::vector<Sector>& v, const QN& q)
{
  for(size_t i = 0; i < v.size(); i++) if(v[i].qn == q) return int(i);
  return -1;
}

// The sector holding state i.
inline int containing(const std::vector<Sector>& v, int i)
{
  for(size_t k = 0; k < v.size(); k++) if(i >= v[k].begin && i < v[k].begin + v[k].dim) return int(k);
  return -1;
}

using Mat = Eigen::MatrixXd;

inline Eigen::Map<const Mat> view(const Sub& b)
  { return Eigen::Map<const Mat>(b.array(), Eigen::Index(b.rows()), Eigen::Index(b.cols())); }

// The final state with the centre on one site (see the file comment):
// m[a * d + s] holds C[l,s,r] for l in left sector a, r in right sector
// rs[a * d + s] (-1: no state of that charge).
struct Centre
{
  std::vector<Sector> left, right;
  std::vector<QN> site;  // charge of each site state
  int d = 0;
  std::vector<Mat> m;
  std::vector<int> rs;

  Centre(std::vector<Sector> l, std::vector<Sector> r, const Basis& sb, const QN& target)
    : left(std::move(l)), right(std::move(r)), d(int(sb.size()))
    {
      for(int s = 0; s < d; s++) site.push_back(sb[s].qn());
      m.resize(left.size() * d);
      rs.assign(left.size() * d, -1);
      for(size_t a = 0; a < left.size(); a++)
        for(int s = 0; s < d; s++){
          int const x = sector_of(right, target - left[a].qn - site[s]);
          rs[a * d + s] = x;
          if(x >= 0) m[a * d + s] = Mat::Zero(left[a].dim, right[x].dim);
        }
    }
  Mat& at(int a, int s) { return m[a * d + s]; }
  const Mat& at(int a, int s) const { return m[a * d + s]; }
};

// An operator's blocks with their sector indices (bra, ket).
struct Blocks
{
  struct Ref { const Sub* b; int bra, ket; };
  std::vector<Ref> refs;
  Blocks(const Op& op, const std::vector<Sector>& v)
    {
      for(const Sub& b : op){
        int const ket = containing(v, b.col_range().begin()), bra = containing(v, b.row_range().begin());
        if(ket < 0 || bra < 0) throw std::logic_error("molecular dmtk: operator block outside its basis");
        refs.push_back({&b, bra, ket});
      }
    }
};

} // namespace molecular

class MolecularSystem: public System<double>
{
  public:
    using Op = molecular::Op;
    using Set = molecular::Set;
    using Dense = molecular::Dense;
    using Enlarged = molecular::Enlarged;

    // h holds the sites and lattice only; the Hamiltonian is ints. c_up and
    // c_dn are the orbital site's annihilators on its charge-sorted basis.
    MolecularSystem(const Hami<double>& h, const Lattice& lattice, const char* name,
                    MolecularIntegrals ints, const Op& c_up, const Op& c_dn)
      : System<double>(h, lattice, name), _ints(std::move(ints)), _nso(2 * _ints.n)
      {
        const Op* c[2] = {&c_up, &c_dn};
        for(int s = 0; s < 2; s++){
          _dq_c[s] = c[s]->dqn;
          _dq_cd[s] = QN() - c[s]->dqn;
          for(int r = 0; r < 4; r++)
            for(int col = 0; col < 4; col++) _c[s](r, col) = (*c[s])(col, r);
          _cd[s] = molecular::transposed(_c[s]);
        }
        _sites.resize(lattice.size());
      }

    // Expectations of products of two or four single-mode factors, in
    // operator order, in the final state of final_sweep() (full sweep): read
    // at its closing configuration from the left rotations and left block
    // operators of its left-to-right half (see the file comment).
    std::vector<double> measure_products(const std::vector<std::vector<molecular::Factor> >& products);

  protected:
    MolecularIntegrals _ints;
    int _nso;
    Dense _c[2], _cd[2];
    QN _dq_c[2], _dq_cd[2];
    std::vector<Set> _sites;                // site sets, built on first use
    std::vector<int> _site_built;
    Set _left, _right;                      // enlarged blocks of this step
    std::vector<ProductTerm<double> > _terms1, _terms2;

    int lo_of(int site) const { return site < _ints.n ? 2 * site : _nso; }

    const PackedBasis& site_basis(int site) const
      { return h.get_site(site).basis().subspaces(); }

    QN dq_pair_c(int K, int L) const { return _dq_c[K & 1] + _dq_c[L & 1]; }

    // Dense operators of one orbital site's spin orbitals z, z' (z = lo + s).
    const Dense& c(int z) const { return _c[z & 1]; }
    const Dense& cd(int z) const { return _cd[z & 1]; }

    // An operator on a site from its dense matrix and charge change; throws
    // if a nonzero element does not have that charge change.
    Op site_op(const Dense& m, const QN& dqn, int site) const
      {
        const PackedBasis& basis = site_basis(site);
        Op op = molecular::make_op(basis, dqn);
        std::array<bool, 16> placed{};
        for(molecular::Sub& b : op){
          int const c0 = b.col_range().begin(), r0 = b.row_range().begin();
          for(size_t i = 0; i < b.cols(); i++)
            for(size_t j = 0; j < b.rows(); j++){
              b(i, j) = m(r0 + int(j), c0 + int(i));
              placed[(r0 + j) * 4 + c0 + i] = true;
            }
        }
        for(int k = 0; k < 16; k++)
          if(m.m[k] != 0. && !placed[k]) throw std::logic_error("molecular dmtk: site operator has mixed charge changes");
        molecular::prune(op);
        return op;
      }

    // Dense site parts used by the complementary operators.
    Dense dense_P(int K, int L, int lo) const
      {
        Dense m;
        if(lo >= _nso) return m; // padding site
        for(int z = lo; z < lo + 2; z++)
          for(int zp = lo; zp < lo + 2; zp++){
            double v = _ints.w(K, L, z, zp);
            if(v != 0.) m.add(molecular::mul(c(z), c(zp)), v);
          }
        return m;
      }
    Dense dense_Q(int K, int L, int lo) const
      {
        Dense m;
        if(lo >= _nso) return m; // padding site
        for(int z = lo; z < lo + 2; z++)
          for(int zp = lo; zp < lo + 2; zp++){
            double v = _ints.g(K, L, z, zp);
            if(v != 0.) m.add(molecular::mul(cd(z), c(zp)), v);
          }
        return m;
      }
    Dense dense_S(int K, int lo) const
      {
        Dense m;
        if(lo >= _nso) return m; // padding site
        for(int z = lo; z < lo + 2; z++) m.add(cd(z), .5 * _ints.t(z, K));
        for(int P = lo; P < lo + 2; P++)
          for(int Q = lo; Q < lo + 2; Q++){
            if(P == Q) continue;
            Dense pq = molecular::mul(cd(P), cd(Q));
            for(int R = lo; R < lo + 2; R++){
              double v = _ints.w(P, Q, R, K);
              if(v != 0.) m.add(molecular::mul(pq, c(R)), v);
            }
          }
        return m;
      }
    Dense dense_H(int lo) const
      {
        Dense m;
        if(lo >= _nso) return m; // padding site
        for(int P = lo; P < lo + 2; P++)
          for(int Q = lo; Q < lo + 2; Q++){
            double v = _ints.t(P, Q);
            if(v != 0.) m.add(molecular::mul(cd(P), c(Q)), v);
          }
        for(int P = lo; P < lo + 2; P++)
          for(int Q = lo; Q < lo + 2; Q++){
            if(P == Q) continue;
            Dense pq = molecular::mul(cd(P), cd(Q));
            for(int R = lo; R < lo + 2; R++)
              for(int S = lo; S < lo + 2; S++){
                if(R == S) continue;
                double v = _ints.w(P, Q, R, S);
                if(v != 0.) m.add(molecular::mul(molecular::mul(pq, c(R)), c(S)), .5 * v);
              }
          }
        return m;
      }

    // The operators of a one-site block, complementary over all other modes.
    const Set& site_set(int site)
      {
        if(_site_built.empty()) _site_built.assign(_sites.size(), 0);
        Set& X = _sites[site];
        if(_site_built[site]) return X;
        _site_built[site] = 1;
        int const lo = lo_of(site), hi = site < _ints.n ? lo + 2 : _nso;
        X.init(lo, hi, _nso);
        if(lo == hi){ // padding: no modes, zero Hamiltonian
          X.H = molecular::make_op(site_basis(site), QN());
          return X;
        }
        X.H = site_op(dense_H(lo), QN(), site);
        for(int I = lo; I < hi; I++){
          X.a_(I) = site_op(c(I), _dq_c[I & 1], site);
          X.ad_(I) = site_op(cd(I), _dq_cd[I & 1], site);
          for(int J = I; J < hi; J++){
            if(I < J) X.A_(I, J) = site_op(molecular::mul(cd(I), cd(J)), _dq_cd[I & 1] + _dq_cd[J & 1], site);
            X.B_(I, J) = site_op(molecular::mul(cd(I), c(J)), _dq_cd[I & 1] + _dq_c[J & 1], site);
          }
        }
        for(size_t k = 0; k < X.outside.size(); k++){
          int const K = X.outside[k];
          X.S_(K) = site_op(dense_S(K, lo), _dq_cd[K & 1], site);
          for(size_t l = k; l < X.outside.size(); l++){
            int const L = X.outside[l];
            if(K < L) X.P_(K, L) = site_op(dense_P(K, L, lo), dq_pair_c(K, L), site);
            X.Q_(K, L) = site_op(dense_Q(K, L, lo), _dq_c[K & 1] + _dq_cd[L & 1], site);
          }
        }
        return X;
      }

    // P^X[K,L] and Q^X[K,L] for either order of two outside modes.
    static Op P_any(const Set& X, int K, int L)
      {
        if(K < L) return X.P_(K, L);
        Op o = X.P_(L, K);
        for(molecular::Sub& b : o) b *= -1.;
        return o;
      }
    static Op Q_any(const Set& X, int K, int L)
      { return K <= L ? X.Q_(K, L) : molecular::transpose(X.Q_(L, K)); }

    // E += sum_u X[u] M[u] (dot second) or sum_u M[u] X[u] (dot first), the
    // X[u] operators of the block part (one charge change), the M[u] dense
    // on the dot site. Expanded in the site's matrix units |r><c|.
    void family(Op& E, const std::vector<const Op*>& X, const std::vector<Dense>& M,
                int dot, bool dot_first, const Enlarged& en) const
      {
        if(X.empty()) return;
        const Basis& sb = h.get_site(dot).basis();
        for(int r = 0; r < 4; r++)
          for(int col = 0; col < 4; col++){
            size_t first = M.size();
            for(size_t u = 0; u < M.size(); u++) if(M[u](r, col) != 0.) { first = u; break; }
            if(first == M.size()) continue;
            // Charge conservation leaves one charge change per matrix unit.
            Op C = molecular::make_op(X[first]->subspaces(), X[first]->dqn);
            for(size_t u = first; u < X.size(); u++){
              if(M[u](r, col) == 0.) continue;
              if(X[u]->dqn != C.dqn) throw std::logic_error("molecular dmtk: mixed charge changes in a family");
              molecular::axpy(C, *X[u], M[u](r, col));
            }
            Dense unit;
            unit(r, col) = 1.;
            Op U = site_op(unit, sb[r].qn() - sb[col].qn(), dot);
            if(dot_first) molecular::kron(E, U, C, 1., en);
            else molecular::kron(E, C, U, 1., en);
          }
      }

    // The enlarged block E of block X and site `dot`, the site after X
    // (dot_first false: E = X dot) or before it (dot_first true: E = dot X).
    void combine(Set& E, const Set& X, int dot, bool dot_first, const Enlarged& en)
      {
        const Set& D = site_set(dot);
        const Set& F = dot_first ? D : X; // first part
        const Set& G = dot_first ? X : D; // second part
        const PackedBasis& eb = en.basis.subspaces();
        E.init(std::min(X.lo, D.lo), std::max(X.hi, D.hi), _nso);
        auto embed = [&](Op& o, const Set& part, const Op& p) {
          if(&part == &F) molecular::embed_first(o, p, 1., en);
          else molecular::embed_second(o, p, 1., en);
        };
        auto finish = [](Op& o) { molecular::prune(o); };

        // Hamiltonian: both parts' plus F + F^T, F the forward half of their coupling.
        {
          Op Fw = molecular::make_op(eb, QN());
          std::vector<const Op*> xs;
          std::vector<Dense> ms;
          if(!dot_first){
            for(int y = D.lo; y < D.hi; y++) molecular::kron(Fw, X.S_(y), D.a_(y), 1., en);
            for(int x = X.lo; x < X.hi; x++) { xs.push_back(&X.a_(x)); ms.push_back(Dense().add(dense_S(x, D.lo), -1.)); }
            family(Fw, xs, ms, dot, false, en);
            for(int y = D.lo; y < D.hi; y++)
              for(int yp = y; yp < D.hi; yp++){
                if(y < yp) molecular::kron(Fw, X.P_(y, yp), D.A_(y, yp), 1., en);
                molecular::kron(Fw, X.Q_(y, yp), D.B_(y, yp), y < yp ? 1. : .5, en);
              }
          } else {
            for(int y = X.lo; y < X.hi; y++) { xs.push_back(&X.a_(y)); ms.push_back(dense_S(y, D.lo)); }
            family(Fw, xs, ms, dot, true, en);
            for(int x = D.lo; x < D.hi; x++) molecular::kron(Fw, D.a_(x), X.S_(x), -1., en);
            for(int x = D.lo; x < D.hi; x++)
              for(int xp = x; xp < D.hi; xp++){
                if(x < xp) molecular::kron(Fw, D.A_(x, xp), X.P_(x, xp), 1., en);
                molecular::kron(Fw, D.B_(x, xp), X.Q_(x, xp), x < xp ? 1. : .5, en);
              }
          }
          E.H = molecular::make_op(eb, QN());
          molecular::embed_first(E.H, F.H, 1., en);
          molecular::embed_second(E.H, G.H, 1., en);
          molecular::axpy(E.H, Fw, 1.);
          molecular::axpy(E.H, molecular::transpose(Fw), 1.);
          finish(E.H);
        }

        // Single-mode and normal pair operators.
        for(int I = E.lo; I < E.hi; I++){
          const Set& pI = F.inside(I) ? F : G;
          E.a_(I) = molecular::make_op(eb, _dq_c[I & 1]); embed(E.a_(I), pI, pI.a_(I)); finish(E.a_(I));
          E.ad_(I) = molecular::make_op(eb, _dq_cd[I & 1]); embed(E.ad_(I), pI, pI.ad_(I)); finish(E.ad_(I));
          for(int J = I; J < E.hi; J++){
            const Set& pJ = F.inside(J) ? F : G;
            if(I < J){
              Op& o = E.A_(I, J);
              o = molecular::make_op(eb, _dq_cd[I & 1] + _dq_cd[J & 1]);
              if(&pI == &pJ) embed(o, pI, pI.A_(I, J));
              else molecular::kron(o, F.ad_(I), G.ad_(J), 1., en);
              finish(o);
            }
            Op& o = E.B_(I, J);
            o = molecular::make_op(eb, _dq_cd[I & 1] + _dq_c[J & 1]);
            if(&pI == &pJ) embed(o, pI, pI.B_(I, J));
            else molecular::kron(o, F.ad_(I), G.a_(J), 1., en);
            finish(o);
          }
        }

        // Complementary operators for the modes outside E.
        std::vector<const Op*> xa, xad;
        for(int u = X.lo; u < X.hi; u++) { xa.push_back(&X.a_(u)); xad.push_back(&X.ad_(u)); }
        std::vector<Dense> m1(X.hi - X.lo), m2(X.hi - X.lo);
        for(size_t k = 0; k < E.outside.size(); k++){
          int const K = E.outside[k];
          // S[K]
          {
            Op& o = E.S_(K);
            o = molecular::make_op(eb, _dq_cd[K & 1]);
            embed(o, F, F.S_(K));
            embed(o, G, G.S_(K));
            if(!dot_first){
              for(int v = D.lo; v < D.hi; v++){
                molecular::kron(o, molecular::transpose(P_any(X, K, v)), D.a_(v), 1., en);
                molecular::kron(o, Q_any(X, v, K), D.ad_(v), 1., en);
              }
              for(int u = X.lo; u < X.hi; u++){
                m1[u - X.lo] = dense_Q(u, K, D.lo);
                m2[u - X.lo] = molecular::transposed(dense_P(K, u, D.lo));
              }
              family(o, xad, m1, dot, false, en);
              family(o, xa, m2, dot, false, en);
            } else {
              for(int v = X.lo; v < X.hi; v++){
                m1[v - X.lo] = molecular::transposed(dense_P(K, v, D.lo));
                m2[v - X.lo] = dense_Q(v, K, D.lo);
              }
              family(o, xa, m1, dot, true, en);
              family(o, xad, m2, dot, true, en);
              for(int u = D.lo; u < D.hi; u++){
                molecular::kron(o, D.ad_(u), Q_any(X, u, K), 1., en);
                molecular::kron(o, D.a_(u), molecular::transpose(P_any(X, K, u)), 1., en);
              }
            }
            finish(o);
          }
          for(size_t l = k; l < E.outside.size(); l++){
            int const L = E.outside[l];
            // P[K,L]: sum over u in F, v in G of (w(K,L,u,v) - w(K,L,v,u)) a_u a_v
            if(K < L){
              Op& o = E.P_(K, L);
              o = molecular::make_op(eb, dq_pair_c(K, L));
              embed(o, F, F.P_(K, L));
              embed(o, G, G.P_(K, L));
              for(int x = X.lo; x < X.hi; x++){
                Dense& m = m1[x - X.lo];
                m = Dense();
                for(int z = D.lo; z < D.hi; z++){
                  int const u = dot_first ? z : x, v = dot_first ? x : z;
                  m.add(c(z), _ints.w(K, L, u, v) - _ints.w(K, L, v, u));
                }
              }
              family(o, xa, m1, dot, dot_first, en);
              finish(o);
            }
            // Q[K,L]: sum over u in F, v in G of g(K,L,u,v) a+_u a_v - g(K,L,v,u) a_u a+_v
            {
              Op& o = E.Q_(K, L);
              o = molecular::make_op(eb, _dq_c[K & 1] + _dq_cd[L & 1]);
              embed(o, F, F.Q_(K, L));
              embed(o, G, G.Q_(K, L));
              for(int x = X.lo; x < X.hi; x++){
                Dense& ma = m1[x - X.lo]; // multiplies X's a_x
                Dense& mad = m2[x - X.lo]; // multiplies X's a+_x
                ma = Dense(); mad = Dense();
                for(int z = D.lo; z < D.hi; z++){
                  if(!dot_first){ // u = x in X, v = z on the dot
                    mad.add(c(z), _ints.g(K, L, x, z));
                    ma.add(cd(z), -_ints.g(K, L, z, x));
                  } else {        // u = z on the dot, v = x in X
                    ma.add(cd(z), _ints.g(K, L, z, x));
                    mad.add(c(z), -_ints.g(K, L, x, z));
                  }
                }
              }
              family(o, xad, m2, dot, dot_first, en);
              family(o, xa, m1, dot, dot_first, en);
              finish(o);
            }
          }
        }
      }

    std::string file_name(int position, int sites) const
      {
        return std::string("qc_") + name() + "_" + std::to_string(sites) + (position == LEFT ? "_l.dat" : "_r.dat");
      }

    void save(Set& X, int position, int sites)
      {
        std::ofstream out(tmp_files.get_filename(file_name(position, sites).c_str()), std::ios::binary);
        if(!out) throw std::runtime_error("molecular dmtk: cannot write block operators");
        X.write(out);
        out.flush();
        if(!out) throw std::runtime_error("molecular dmtk: writing block operators failed (disk full?)");
      }

    void load(Set& X, int position, int sites)
      {
        std::ifstream in(tmp_files.get_filename(file_name(position, sites).c_str()), std::ios::binary);
        if(!in) throw std::runtime_error("molecular dmtk: cannot read block operators");
        X.read(in);
      }

    // Changed here (hook of System::init_iteration): the enlarged blocks'
    // operators and the superblock product terms of this step.
    void init_hamiltonian(bool) override
      {
        int const n1 = int(_b1->lattice().size()), n4 = int(_b4->lattice().size());
        int const s2 = n1, s3 = n1 + 1, s4 = n1 + 2;
        Set xl, xr;
        if(n1 == 1) xl = site_set(0); else load(xl, LEFT, n1);
        if(n4 == 1) xr = site_set(s4);
        else {
          if(s4 + n4 != int(lattice().size())) throw std::logic_error("molecular dmtk: right block does not end the chain");
          load(xr, RIGHT, n4);
        }
        Enlarged enl(_b1->basis(), _b2->basis()), enr(_b3->basis(), _b4->basis());
        combine(_left, xl, s2, false, enl);
        combine(_right, xr, s3, true, enr);

        VectorState<double> v12 = gs.condense(MASK_BLOCK1|MASK_BLOCK2);
        VectorState<double> v = v12.condense(MASK_BLOCK3|MASK_BLOCK4);
        VectorState<double> r(v);
        _terms1.clear();
        _terms2.clear();
        auto add1 = [&](const Op& o, size_t m) {
          auto t = get_product_terms(o, v, r, m, 1., false);
          _terms1.insert(_terms1.end(), t.begin(), t.end());
        };
        auto add2 = [&](const Op& o1, const Op& o2, double coef, bool hc) {
          if(o1.empty() || o2.empty()) return;
          auto t = get_product_terms(o1, o2, v, r, BLOCK1, BLOCK4, coef, hc);
          _terms2.insert(_terms2.end(), t.begin(), t.end());
        };
        Set& L = _left;
        Set& R = _right;
        add1(L.H, BLOCK1);
        add1(R.H, BLOCK4);
        for(int y = R.lo; y < R.hi; y++) add2(L.S_(y), R.a_(y), 1., true);
        for(int x = L.lo; x < L.hi; x++) add2(L.a_(x), R.S_(x), -1., true);
        if(L.hi - L.lo <= R.hi - R.lo){
          for(int x = L.lo; x < L.hi; x++)
            for(int xp = x; xp < L.hi; xp++){
              if(x < xp) add2(L.A_(x, xp), R.P_(x, xp), 1., true);
              add2(L.B_(x, xp), R.Q_(x, xp), 1., x < xp);
            }
        } else {
          for(int y = R.lo; y < R.hi; y++)
            for(int yp = y; yp < R.hi; yp++){
              if(y < yp) add2(L.P_(y, yp), R.A_(y, yp), 1., true);
              add2(L.Q_(y, yp), R.B_(y, yp), 1., y < yp);
            }
        }
      }

    // Changed here (hook of the Lanczos product): H v on the condensed state.
    VectorState<double> hamiltonian_product(const VectorState<double>& vv) override
      {
        DMTKglobals<double> globals = get_globals(double(0));
        VectorState<double> v12 = vv.condense(MASK_BLOCK1|MASK_BLOCK2);
        VectorState<double> v = v12.condense(MASK_BLOCK3|MASK_BLOCK4);
        VectorState<double> r(v);
        r = 0.;
        for(const ProductTerm<double>& t : _terms1) product_term(t, v, r, MASK_PRODUCT_DEFAULT, &globals);
        for(const ProductTerm<double>& t : _terms2) product_term(t, v, r, MASK_PRODUCT_DEFAULT|MASK_PRODUCT_HC, &globals);
        VectorState<double> r12 = r.decondense(MASK_BLOCK3|MASK_BLOCK4, v12);
        return r12.decondense(MASK_BLOCK1|MASK_BLOCK2, vv);
      }

    // The parity of a charge.
    static bool odd(const QN& q) { return q.fermion_sign() == -1; }

    // The normal product of factors on one site, as a dense matrix.
    Dense dense(const std::vector<molecular::Factor>& f, QN& dq) const
      {
        Dense m;
        for(int i = 0; i < 4; i++) m(i, i) = 1.;
        dq = QN();
        for(const molecular::Factor& x : f){
          m = molecular::mul(m, x.dagger ? cd(x.mode) : c(x.mode));
          dq += x.dagger ? _dq_cd[x.mode & 1] : _dq_c[x.mode & 1];
        }
        return m;
      }

    // Changed here (hook of System::rotate): the grown block's operators
    // rotated into the kept states and filed for the later steps.
    void rotate_hamiltonian(int position, Block<double>& b, Basis&, Basis& rho_basis) override
      {
        Set& E = position == LEFT ? _left : _right;
        Set N;
        N.init(E.lo, E.hi, E.nso);
        const PackedBasis& kept = rho_basis.subspaces();
        std::vector<Op*> dst, src;
        N.each([&](Op& o) { dst.push_back(&o); });
        E.each([&](Op& o) { src.push_back(&o); });
        for(size_t k = 0; k < dst.size(); k++) *dst[k] = molecular::rotate(*src[k], rho, kept);
        // The rotation keeps H symmetric only to rounding.
        Op sym = molecular::make_op(kept, QN());
        molecular::axpy(sym, N.H, .5);
        molecular::axpy(sym, molecular::transpose(N.H), .5);
        molecular::prune(sym);
        N.H = sym;
        save(N, position, int(b.lattice().size()));
      }
};

inline std::vector<double>
MolecularSystem::measure_products(const std::vector<std::vector<molecular::Factor> >& products)
{
  using namespace molecular;
  int const ls = int(lattice().size()), last = ls - 3;
  if(int(_b1->lattice().size()) != last || int(_b4->lattice().size()) != 1)
    throw std::logic_error("molecular dmtk: products are measured at the closing configuration of final_sweep");
  if(gs.subspace_begin() == gs.subspace_end())
    throw std::logic_error("molecular dmtk: no final state to measure");
  const StateSpace& first = *gs.subspace_begin();
  QN const target = first[1].qn() + first[2].qn() + first[3].qn() + first[4].qn();

  // The position of each product: its third-lowest site (second for two
  // factors), clamped to [1, ls-3].
  std::vector<double> value(products.size(), 0.);
  std::vector<std::vector<size_t> > at(ls);
  for(size_t i = 0; i < products.size(); i++){
    const std::vector<Factor>& f = products[i];
    if(f.size() != 2 && f.size() != 4)
      throw std::invalid_argument("molecular dmtk: a measured product has two or four factors");
    std::vector<int> s;
    for(const Factor& x : f){
      if(x.mode < 0 || x.mode >= _nso) throw std::invalid_argument("molecular dmtk: factor outside the spin orbitals");
      s.push_back(x.mode >> 1);
    }
    std::sort(s.begin(), s.end());
    at[std::min(std::max(s[std::min<size_t>(2, s.size() - 1)], 1), last)].push_back(i);
  }

  // The closing configuration: left block of ls-3 sites, site ls-3 and the
  // last two sites as the right block.
  Enlarged e2(h.get_site(ls - 2).basis(), h.get_site(ls - 1).basis());
  Basis rbasis = e2.basis;
  Centre C(sectors(_b1->basis().subspaces()), sectors(rbasis.subspaces()), h.get_site(last).basis(), target);
  for(auto it = gs.subspace_begin(); it != gs.subspace_end(); ++it){
    const StateSpace& sp = *it;
    if(sp.start() == -1) continue;
    int const a = containing(C.left, sp[1].begin()), s = sp[2].begin();
    int const r = e2.pos[sp[3].begin() * e2.d2 + sp[4].begin()];
    int const x = a < 0 ? -1 : C.rs[a * C.d + s];
    if(x < 0 || C.left[a].begin != sp[1].begin() || C.left[a].dim != int(sp[1].dim()) || containing(C.right, r) != x)
      throw std::logic_error("molecular dmtk: final state outside the closing configuration");
    auto slice = gs(sp);
    for(int l = 0; l < C.left[a].dim; l++) C.at(a, s)(l, r - C.right[x].begin) = slice(l, 0, 0, 0);
  }
  // Annihilators of the right block's spin orbitals (a+ is the transpose).
  std::vector<Op> rsing(_nso);
  for(int site = ls - 2; site < ls; site++){
    const Set& D = site_set(site);
    for(int I = D.lo; I < D.hi; I++){
      Op E = make_op(rbasis.subspaces(), _dq_c[I & 1]);
      if(site == ls - 2) embed_first(E, D.a_(I), 1., e2);
      else embed_second(E, D.a_(I), 1., e2);
      prune(E);
      rsing[I] = E;
    }
  }

  for(int k = last; k >= 1; k--){
    Set X;
    if(k >= 2) load(X, LEFT, k);
    std::deque<Op> owned; // products built here: the left block's at k = 1, the right's at k = ls-3
    std::map<std::vector<int>, const Op*> built;
    struct RightOp { const Op* op; bool t; };
    struct Job { size_t product; const Op* op; bool t, parity; Dense dot; double sign; };
    std::vector<RightOp> rights;
    std::vector<std::vector<Job> > jobs;
    std::map<std::vector<int>, size_t> right_index;
    auto code = [](const std::vector<Factor>& f) {
      std::vector<int> c;
      for(const Factor& x : f) c.push_back(2 * x.mode + int(x.dagger));
      return c;
    };

    for(size_t i : at[k]){
      const std::vector<Factor>& f = products[i];
      // Group the factors by part (left block, site k, right block); the
      // reordering's sign counts the factors that move past each other.
      std::vector<Factor> fl, fd, fr;
      int moved = 0;
      for(const Factor& x : f){
        int const site = x.mode >> 1;
        if(site < k){ fl.push_back(x); moved += int(fd.size() + fr.size()); }
        else if(site == k){ fd.push_back(x); moved += int(fr.size()); }
        else fr.push_back(x);
      }
      double sign = moved % 2 ? -1. : 1.;
      if((k >= 2 && fl.size() > 2) || (k < last && fr.size() > 1))
        throw std::logic_error("molecular dmtk: product does not fit its position");

      QN dq;
      Dense dot = dense(fd, dq);
      if(fr.size() % 2)
        for(int r = 0; r < 4; r++)
          for(int s = 0; s < C.d; s++) if(odd(C.site[s])) dot(r, s) = -dot(r, s);
      bool nonzero = false;
      for(int r = 0; r < C.d; r++) for(int s = 0; s < C.d; s++) if(dot(r, s) != 0.) nonzero = true;
      if(!nonzero) continue;

      // Left block: its normal product of fl.
      const Op* lop = nullptr;
      bool lt = false;
      if(!fl.empty() && k == 1){
        auto key = code(fl);
        auto found = built.find(key);
        if(found == built.end()){
          QN q;
          Dense m = dense(fl, q);
          owned.push_back(site_op(m, q, 0));
          found = built.emplace(key, &owned.back()).first;
        }
        lop = found->second;
      } else if(fl.size() == 1){
        lop = fl[0].dagger ? &X.ad_(fl[0].mode) : &X.a_(fl[0].mode);
      } else if(fl.size() == 2){
        int const I = fl[0].mode, J = fl[1].mode;
        bool const dI = fl[0].dagger, dJ = fl[1].dagger;
        if(I == J && dI == dJ) continue; // a+ a+ or a a on one mode
        if(dI && dJ){        // a+_I a+_J
          lop = I < J ? &X.A_(I, J) : &X.A_(J, I);
          if(I > J) sign = -sign;
        } else if(dI){       // a+_I a_J
          lop = I <= J ? &X.B_(I, J) : &X.B_(J, I);
          lt = I > J;
        } else if(!dJ){      // a_I a_J = (a+_J a+_I)^T
          lop = J < I ? &X.A_(J, I) : &X.A_(I, J);
          lt = true;
          if(I < J) sign = -sign;
        } else {             // a_I a+_J = -a+_J a_I
          if(I == J) throw std::logic_error("molecular dmtk: a a+ on one mode is not a stored block operator");
          lop = J <= I ? &X.B_(J, I) : &X.B_(I, J);
          lt = J > I;
          sign = -sign;
        }
      }

      // Right block: its normal product of fr.
      auto key = code(fr);
      auto found = right_index.find(key);
      if(found == right_index.end()){
        RightOp R{nullptr, false};
        if(fr.empty()) {}
        else if(k < last) R = {&rsing[fr[0].mode], fr[0].dagger};
        else {
          // Exact on the last two sites: site ls-2's factors, then ls-1's.
          std::vector<Factor> f2, f1;
          int swaps = 0;
          for(const Factor& x : fr){
            if((x.mode >> 1) == ls - 2){ f2.push_back(x); swaps += int(f1.size()); }
            else f1.push_back(x);
          }
          QN q2, q1;
          Dense const d2 = dense(f2, q2), d1 = dense(f1, q1);
          double const sg = swaps % 2 ? -1. : 1.;
          Op o = make_op(rbasis.subspaces(), q2 + q1);
          if(!f2.empty() && !f1.empty()) kron(o, site_op(d2, q2, ls - 2), site_op(d1, q1, ls - 1), sg, e2);
          else if(!f2.empty()) embed_first(o, site_op(d2, q2, ls - 2), sg, e2);
          else embed_second(o, site_op(d1, q1, ls - 1), sg, e2);
          prune(o);
          owned.push_back(o);
          R = {&owned.back(), false};
        }
        found = right_index.emplace(key, rights.size()).first;
        rights.push_back(R);
        jobs.emplace_back();
      }
      jobs[found->second].push_back({i, lop, lt, fl.size() % 2 == 1, dot, sign});
    }

    // Each right operator once: Z = O_R C_s'^T per left sector, Y = C_s Z,
    // then every product's Frobenius products with Y.
    int const d = C.d, nl = int(C.left.size());
    std::map<const Op*, Blocks> lblocks;
    using Part = std::pair<int, Mat>; // (sector of the bra, block)
    for(size_t w = 0; w < rights.size(); w++){
      const RightOp& R = rights[w];
      bool need[4][4] = {};
      bool needp[4] = {};
      for(const Job& j : jobs[w])
        for(int s = 0; s < d; s++)
          for(int sp = 0; sp < d; sp++) if(j.dot(s, sp) != 0.) need[s][sp] = needp[sp] = true;
      std::vector<std::vector<Part> > Z(d);
      std::vector<Blocks::Ref> rrefs;
      if(R.op) rrefs = Blocks(*R.op, C.right).refs;
      for(int sp = 0; sp < d; sp++){
        if(!needp[sp]) continue;
        Z[sp].assign(nl, Part(-1, Mat()));
        for(int a = 0; a < nl; a++){
          int const x = C.rs[a * d + sp];
          if(x < 0) continue;
          const Mat& B = C.at(a, sp);
          if(!R.op){ Z[sp][a] = Part(x, B.transpose()); continue; }
          for(const Blocks::Ref& ref : rrefs){
            if(!R.t && ref.ket == x){ Z[sp][a] = Part(ref.bra, view(*ref.b) * B.transpose()); break; }
            if(R.t && ref.bra == x){ Z[sp][a] = Part(ref.ket, view(*ref.b).transpose() * B.transpose()); break; }
          }
        }
      }
      std::vector<std::vector<std::vector<Part> > > Y(d, std::vector<std::vector<Part> >(d));
      for(int s = 0; s < d; s++)
        for(int sp = 0; sp < d; sp++){
          if(!need[s][sp]) continue;
          Y[s][sp].assign(nl, Part(-1, Mat()));
          for(int ap = 0; ap < nl; ap++){
            const Part& z = Z[sp][ap];
            if(z.first < 0) continue;
            int const a = sector_of(C.left, target - C.right[z.first].qn - C.site[s]);
            if(a < 0 || C.rs[a * d + s] != z.first) continue;
            Y[s][sp][ap] = Part(a, C.at(a, s) * z.second);
          }
        }
      for(const Job& j : jobs[w]){
        const Blocks* lb = nullptr;
        if(j.op){
          auto it = lblocks.find(j.op);
          if(it == lblocks.end()) it = lblocks.emplace(j.op, Blocks(*j.op, C.left)).first;
          lb = &it->second;
        }
        double v = 0.;
        for(int s = 0; s < d; s++)
          for(int sp = 0; sp < d; sp++){
            double const x = j.dot(s, sp);
            if(x == 0.) continue;
            const std::vector<Part>& y = Y[s][sp];
            double t = 0.;
            if(!lb){
              for(int a = 0; a < nl; a++) if(y[a].first == a) t += y[a].second.trace();
            } else for(const Blocks::Ref& ref : lb->refs){
              int const ket = j.t ? ref.bra : ref.ket, bra = j.t ? ref.ket : ref.bra;
              if(y[ket].first != bra) continue;
              double u = j.t ? (view(*ref.b).transpose().array() * y[ket].second.array()).sum()
                             : (view(*ref.b).array() * y[ket].second.array()).sum();
              if(j.parity && odd(C.left[ket].qn)) u = -u;
              t += u;
            }
            v += x * t;
          }
        value[j.product] += j.sign * v;
      }
    }
    if(k == 1) break;

    // Move the centre to site k-1. The right block of sites k..ls-1 gets an
    // orthonormal basis of the centre's rows over (site k, right block), one
    // left sector (one charge) at a time: C = L Q^T by Householder QR, all
    // min(rows, cols) columns of Q kept.
    Enlarged en(h.get_site(k).basis(), rbasis);
    std::vector<Sector> const es = sectors(en.basis.subspaces());
    std::vector<std::pair<int, int> > order; // (enlarged sector, left sector)
    for(int a = 0; a < nl; a++){
      int const e = sector_of(es, target - C.left[a].qn);
      if(e >= 0) order.push_back({e, a});
    }
    std::sort(order.begin(), order.end());
    BMatrix<double> rot;
    std::vector<State> states;
    std::vector<Mat> L(nl);
    std::vector<int> nsec(nl, -1);
    int nsectors = 0;
    for(const auto& [e, a] : order){
      Mat M = Mat::Zero(C.left[a].dim, es[e].dim);
      bool any = false;
      for(int s = 0; s < d; s++){
        int const x = C.rs[a * d + s];
        if(x < 0) continue;
        any = true;
        const Mat& B = C.at(a, s);
        for(int r = 0; r < C.right[x].dim; r++)
          M.col(en.pos[s * en.d2 + C.right[x].begin + r] - es[e].begin) = B.col(r);
      }
      if(!any) continue;
      Eigen::Index const p = std::min(M.rows(), M.cols());
      Eigen::HouseholderQR<Mat> qr(M.transpose());
      Mat const V = qr.householderQ() * Mat::Identity(M.cols(), p);
      Mat const R = qr.matrixQR().topRows(p).triangularView<Eigen::Upper>();
      L[a] = R.transpose();
      int const begin = int(states.size());
      SubMatrix<double> blk(es[e].qn, SubSpace(es[e].qn, begin, begin + int(p) - 1),
                            SubSpace(es[e].qn, es[e].begin, es[e].begin + es[e].dim - 1));
      for(Eigen::Index i = 0; i < p; i++)
        for(Eigen::Index j = 0; j < V.rows(); j++) blk(size_t(i), size_t(j)) = V(j, i);
      rot.push_back(blk);
      for(Eigen::Index i = 0; i < p; i++) states.push_back(State(size_t(begin + i), es[e].qn));
      nsec[a] = nsectors++;
    }
    if(states.empty()) throw std::logic_error("molecular dmtk: the final state has no amplitude");
    Basis nb(states.size());
    for(size_t i = 0; i < states.size(); i++) nb[i] = states[i];
    nb.pack();
    rot.repack(nb);
    const Set& D = site_set(k);
    std::vector<Op> ns(_nso);
    for(int I = lo_of(k); I < _nso; I++){
      Op E = make_op(en.basis.subspaces(), _dq_c[I & 1]);
      if(D.inside(I)) embed_first(E, D.a_(I), 1., en);
      else embed_second(E, rsing[I], 1., en);
      prune(E);
      ns[I] = molecular::rotate(E, rot, nb.subspaces());
    }
    rsing.swap(ns);
    rbasis = nb;

    // The left block of k sites is (block of k-1 sites) x (site k-1) under
    // the final rotation U: C'[(l', s'), r] = sum_l U[(l', s'), l] L[l, r].
    BMatrix<double> U;
    Basis ub;
    read_rho(U, ub, k, LEFT);
    std::vector<Sector> lsec;
    for(size_t i = 0; i < ub.size(); i++){
      State const st = ub[i];
      int const x = sector_of(lsec, st.qn1());
      if(x < 0){ lsec.push_back({st.qn1(), int(st.i1), 1}); continue; }
      int const lo = std::min(lsec[x].begin, int(st.i1));
      int const hi = std::max(lsec[x].begin + lsec[x].dim - 1, int(st.i1));
      lsec[x].begin = lo;
      lsec[x].dim = hi - lo + 1;
    }
    std::sort(lsec.begin(), lsec.end(), [](const Sector& u, const Sector& v) { return u.begin < v.begin; });
    Centre N(lsec, sectors(nb.subspaces()), h.get_site(k - 1).basis(), target);
    for(int a = 0; a < nl; a++){
      if(nsec[a] < 0) continue;
      const SubMatrix<double>* u = U.block(C.left[a].qn);
      if(!u || u->col_range().begin() != C.left[a].begin || int(u->cols()) != C.left[a].dim)
        throw std::logic_error("molecular dmtk: rotation does not match the left block");
      Mat const G = view(*u) * L[a];
      int const r0 = u->row_range().begin();
      for(Eigen::Index j = 0; j < G.rows(); j++){
        State const st = ub[size_t(r0 + j)];
        int const ap = containing(N.left, int(st.i1)), sp = int(st.i2);
        if(ap < 0 || N.rs[ap * N.d + sp] != nsec[a])
          throw std::logic_error("molecular dmtk: rotation rows do not match the centre");
        N.at(ap, sp).row(int(st.i1) - N.left[ap].begin) = G.row(j);
      }
    }
    C = std::move(N);
  }
  return value;
}

} // namespace dmtk

#endif // DMTK_MOLECULAR_SYSTEM_H_
