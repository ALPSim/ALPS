/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
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

// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/lanczos.ipp, the Eigen-based
// rework of src/apps/dmrg/dmrg/dmtk/lanczos.cc on master (Adrian Feiguin's
// dmtk). Adaptation: "dmtk/" include paths. The changes made here are marked
// where they occur.

#ifndef __DMTK_LANCZOS_CC__
#define __DMTK_LANCZOS_CC__

#include "dmtk/output.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "dmtk/range.hpp"
#include "dmtk/vector.hpp"
#include "dmtk/matrix.hpp"

namespace dmtk
{

using namespace std;

#define SIGN(a,b) ((b) < 0 ? -fabs(a) : fabs(a))
#define SWAP(a,b) { itemp=(a);(a)=(b);(b)=itemp; }

template<class T>
class AuxIndex
{
  public:
    T v;
    size_t index;

    AuxIndex():v(T()), index(0) {}
    AuxIndex(const T& _v, const size_t &i): v(_v), index(i) {}
    AuxIndex & operator=(const AuxIndex& o) { v = o.v; index = o.index; return *this; } 

    bool operator>(const AuxIndex &o) const
      {
        if(v > o.v) return true;
        if(v == o.v && index > o.index) return true;
        return false;
      }
    bool operator<(const AuxIndex &o) const
      {
        if(v < o.v) return true;
        if(v == o.v && index < o.index) return true;
        return false;
      }
    bool operator>=(const AuxIndex &o) const
      {
        return (v >= o.v && index >= o.index); 
      }
    bool operator<=(const AuxIndex &o) const
      {
        return(v <= o.v && index <= o.index);
      }
    bool operator==(const AuxIndex &o) const
      {
        return (v == o.v);
      }
    bool operator!=(const AuxIndex &o) const
      {
        return (v != o.v);
      }
};


template <class T>
int sort2(size_t n, Vector<T>& arr, bool ascending = true)
{
  T itemp;
  int nswaps = 0;
  bool changed = true;

  if(n == 0){
    dmtk::out() << "*** WARNING: sort width (n == 0) \n";
    return 0;
  }

  Vector<AuxIndex<T> > aux(n);
  AuxIndex<T> *ra = &aux[0];
  for(size_t i = 0; i < n; i++) { ra->v = arr[i]; ra->index = i; ra++; }

  while(changed){
    changed = false;
    for(size_t l = 0; l < n-1; l++){
      if(aux[l] > aux[l+1]){
        SWAP(aux[l],aux[l+1]);
        nswaps++;
        changed = true;
      }
    }
  }

  if(!ascending && n > 1)
    for(size_t j = 0; j <= n/2-1; j++)
      if(aux[j] != aux[n-j-1]) {
        SWAP(aux[j],aux[n-j-1]);
        nswaps++;
      }

  for(size_t i = 0; i < n; i++) arr[i] = aux[i].v;
  return nswaps;
}

template <class T, class A>
int indexx2(size_t n, const A& arr, Vector<size_t>& indx, bool ascending = true)
{
  size_t itemp;
  bool changed = true;
  int nswaps = 0;

  Vector<AuxIndex<T> > aux(n);
  AuxIndex<T> *ra = &aux[0];
  for(size_t i = 0; i < n; i++) { ra->v = arr[i]; ra->index = i; ra++; }

  for (size_t j = 0; j < n; j++) indx[j] = j;
  while(changed){
    changed = false;
    for(size_t l = 0; l < n-1; l++){
      if(aux[indx[l]] > aux[indx[l+1]]){
        SWAP(indx[l],indx[l+1]);
        nswaps++;
        changed = true;
      }
    }
  }

  if(!ascending && n > 1)
    for(size_t j = 0; j <= n/2-1; j++){
      SWAP(indx[j],indx[n-j-1]);
      nswaps++;
    }

  return nswaps;
}

template <class T>
void sort(size_t n, Vector<T>& arr, bool ascending = true)
{
  if (n < 2) return;

  T itemp;
                                                                              
  std::vector<AuxIndex<T> > aux(n);
  AuxIndex<T> *ra = &aux[0];
  for(int i = 0; i < n; i++) { ra->v = arr[i]; ra->index = i; ra++; } 

  std::sort(aux.begin(), aux.end());

  ra = &aux[0];
  for(int i = 0; i < n; i++, ra++) arr[i] = ra->v;  

  if(!ascending && n > 1)
    for(int j = 0; j <= n/2-1; j++)
      if(arr[j] != arr[n-j-1]) SWAP(arr[j],arr[n-j-1]);
}

template<class T, class A>
void indexx(size_t n, const A& arr, Vector<size_t>& indx, bool ascending = true)
{
  if (n == 1) { indx[0] = 0; return; }

  size_t itemp;
                                                                              
  std::vector<AuxIndex<T> > aux(n);
  AuxIndex<T> *ra = &aux[0];
  for(int i = 0; i < n; i++) { ra->v = arr[i]; ra->index = i; ra++; } 
  std::sort(aux.begin(), aux.end());

  ra = &aux[0];
  for(int i = 0; i < n; i++, ra++) indx[i] = ra->index;  

  if(!ascending && n > 1)
    for(int j = 0; j <= n/2-1; j++)
      SWAP(indx[j],indx[n-j-1]);
}

void tqli(dmtk::Vector<double>& d, dmtk::Vector<double>& e,
          int n, dmtk::Matrix<double>& z, bool calc_vectors = false)
{
  // Symmetric tridiagonal eigendecomposition via Eigen, replacing the
  // old dstevd call. Eigen has no specialised tridiagonal solver, so we
  // assemble a dense n×n symmetric tridiagonal matrix from d (diagonal)
  // and e (sub-/super-diagonal) and use SelfAdjointEigenSolver. For
  // small n (Lanczos sub-iteration sizes ≤ ~100) this is fine.
  if(n <= 0) return;

  using EMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  EMat T = EMat::Zero(n, n);
  for(int i = 0; i < n; ++i) T(i, i) = d(i);
  for(int i = 0; i < n - 1; ++i) {
    T(i, i + 1) = e(i);
    T(i + 1, i) = e(i);
  }

  int options = calc_vectors ? Eigen::ComputeEigenvectors
                             : Eigen::EigenvaluesOnly;
  Eigen::SelfAdjointEigenSolver<EMat> solver(T, options);
  if(solver.info() != Eigen::Success) {
    dmtk::out() << "ERROR in dstevd\n";
    return;
  }

  // Write back eigenvalues into d (overwrites diagonal, matching the
  // in-place behaviour of dstevd).
  auto const& evals = solver.eigenvalues();
  for(int i = 0; i < n; ++i) d(i) = evals(i);

  // Write eigenvectors into z (column-major) if requested. The dmtk
  // Matrix z is column-major flat storage of dimension (z.rows × z.cols);
  // we copy the n×n top-left block via an Eigen::Map with the leading
  // dimension equal to z.rows() (which the caller sets to the full
  // workspace stride).
  if(calc_vectors) {
    using MatMap = Eigen::Map<EMat, 0, Eigen::OuterStride<>>;
    MatMap z_view(z.array(), n, n,
                  Eigen::OuterStride<>(static_cast<Eigen::Index>(z.rows())));
    z_view = solver.eigenvectors();
  }
}

//////////////////////////////////////////////////////////////////////

template<class T, class A, class B>
void lanczos(A& m,
             Vector<B>& gs, const B& seed, 
             Vector<double>& ener, int nvectors,
             dmtk::Vector<double>& a, dmtk::Vector<double>& b, 
             int &maxiter, double tol,
             bool use_seed, bool calc_gs, const char* vectors_file, bool force_maxiter,
             bool reorthogonalize, double residual_tol)
{  
  // Stop rule (written for this tree; the legacy engine stopped once the
  // lowest Ritz energy had changed by at most tol, the alpsx rework once
  // every requested Ritz residual was at most tol, about twice the
  // iterations): the Krylov space is accepted once every one of the nvectors
  // lowest Ritz energies has changed by at most tol since the previous
  // iteration and, when residual_tol >= 0, every one of their residuals
  // |b(iter) z(n, iter-1)| is at most residual_tol. The caller chooses the
  // residual bound per solve: sqrt(tol) for the solves that only feed a
  // density matrix, tol for the solves whose states are measured.
  B x1(gs(0)), x2(gs(0)), aux(gs(0));
  dmtk::Vector<double> d(101), e(101), xc(101);
  dmtk::Matrix<double> z(101,101);
  dmtk::Vector<double> previous(std::max(nvectors, 1));
  uint j, col, iter = 0;
  // uint i;
  double e0;
  int control_max = maxiter;
  maxiter = 0;

  if(control_max == 0) { gs[0] = T(1); maxiter = 1; return; }
#ifdef USE_LANCZOS_VECTORS
  if(reorthogonalize)
    throw std::invalid_argument("Lanczos reorthogonalization requires stored basis vectors");
#endif
#ifndef USE_LANCZOS_VECTORS
  ofstream outputfile(vectors_file,std::ios::out|std::ios::binary);
  if(reorthogonalize) outputfile.exceptions(std::ios::failbit | std::ios::badbit);
  if(!outputfile) dmtk::err() << "*** ERROR 1: Lanczos could not open vectors.dat\n";
  ifstream basis;
  if(reorthogonalize) {
    basis.exceptions(std::ios::failbit | std::ios::badbit);
    basis.open(vectors_file, std::ios::in | std::ios::binary);
  }
  // Every Krylov vector has the same basis metadata, so the temporary file
  // stores only the numerical values of each one (alpsx did so for the
  // reorthogonalized solves; here for every solve): no basis metadata is
  // filed per iteration, and no constrained space is reconstructed per read.
  const auto read_basis = [&](B& value, std::istream& input) {
    size_t dimension = 0;
    input.read(reinterpret_cast<char*>(&dimension), sizeof(dimension));
    if(dimension != value.size())
      throw std::runtime_error("Lanczos stored basis dimension mismatch");
    input.read(reinterpret_cast<char*>(value.array()),
               static_cast<std::streamsize>(dimension * sizeof(T)));
  };
  const auto write_basis = [&](const B& value) {
    const size_t dimension = value.size();
    outputfile.write(reinterpret_cast<const char*>(&dimension), sizeof(dimension));
    outputfile.write(reinterpret_cast<const char*>(value.array()),
                     static_cast<std::streamsize>(dimension * sizeof(T)));
  };
#endif
  // The recurrence runs in place on the raw arrays: the materializing
  // VectorState operators would copy the vector and its basis metadata for
  // every temporary of every iteration.
  const size_t dimension = x1.size();
  const auto axpy = [dimension](const T alpha, const B& x, B& y) {
    const T* px = x.array();
    T* py = y.array();
    for(size_t c = 0; c < dimension; ++c) py[c] += alpha*px[c];
  };
 
  e0 = std::numeric_limits<double>::infinity();
  previous = std::numeric_limits<double>::infinity();
  x1 = x2 = aux = T(0);

  while(true) // Two iterations: 1) calculate energy 2) build gs
  { 
    a = 0.0f;
    b = 0.0f;
    if(use_seed){
      x1 = seed;
    } else { 
#ifdef USE_LANCZOS_BASIC
      x1 = 1.;
//      x1 = 0;
//      x1[0] = 1.;
#else
      x1.randomize();
#endif
    } 

    iter = 0;

    b(0) = sqrt(real(dot_product(x1,x1)));
    x1 /= T(b(0));
    x2 = T(0);
    b(0) = 1.;


    uint nmax = std::min(99, static_cast<int>(gs(0).size()));
    bool previous_candidate_rejected = false;
    for(iter = 1; iter <= nmax; iter++){ // Lanczos iteration

      if(b(iter - 1) != 0.){
        // (x1, x2) <- (-b x2, x1 / b)
        const T beta = T(b(iter-1));
        T* p1 = x1.array();
        T* p2 = x2.array();
        for(size_t c = 0; c < dimension; ++c) {
          const T previous_x1 = p1[c];
          p1[c] = -beta*p2[c];
          p2[c] = previous_x1/beta;
        }
      }

      aux = product(m,x2);
      axpy(T(1), aux, x1);
      T xa = dot_product(x1, x2);
      a(iter) = real(xa);
      axpy(-T(a(iter)), x2, x1);
      // Cancellation near an eigenvector can leave a tiny parallel component
      // which division by the small next beta would amplify. Project once more
      // locally; retain its real diagonal correction in the recurrence.
      const T correction = dot_product(x2, x1);
      axpy(-correction, x2, x1);
      a(iter) += real(correction);
#ifndef USE_LANCZOS_VECTORS
      if(reorthogonalize && iter > 1) {
        // Projected excited-state solves have a large spectral penalty. Keep
        // their basis orthogonal without retaining all Krylov vectors in RAM.
        outputfile.flush();
        for(int pass = 0; pass < 2; ++pass) {
          basis.seekg(0);
          for(uint previous = 1; previous < iter; ++previous) {
            read_basis(aux, basis);
            const T overlap = dot_product(aux, x1);
            axpy(-overlap, aux, x1);
          }
          const T parallel = dot_product(x2, x1);
          axpy(-parallel, x2, x1);
          a(iter) += real(parallel);
        }
      }
#endif
      b(iter) = sqrt(real(dot_product(x1, x1)));
   
//      cout << "Iter=" << iter << " " << xa << " " << a(iter) << " " << b(iter) << endl;

#ifndef USE_LANCZOS_VECTORS
      if(calc_gs || reorthogonalize) write_basis(x2);
#endif

      if(maxiter > 0) {  // We are building the ground state;

#ifdef USE_LANCZOS_VECTORS
        gs = gs + xc(iter-1) * x2;
        if(iter == maxiter) return; 
#endif
      } else{  // we are calculating energy

        // A finite space includes its last Krylov vector. Breakdown at the
        // first product must also be resolved before attempting another step.
        const bool exhausted = iter == nmax || iter >= control_max;
        const bool breakdown = b(iter) == 0. || (iter >= 2 && fabs(b(iter)) < tol);
        if(iter >= 2 || exhausted || breakdown){
          d(Range(0,iter-1)) = a(Range(1,iter));
          e(Range(0,iter-1)) = b(Range(1,iter));
          // Once the energies have settled, residual checks usually need
          // vectors on consecutive iterations. Compute them once at that
          // point, and at mandatory reconstruction stops, rather than
          // solving twice.
          bool ritz_vectors_ready = previous_candidate_rejected ||
              (calc_gs && (exhausted || breakdown));
          tqli(d, e, iter, z, ritz_vectors_ready);
  
          e0 = d(0);
          col = 0;
          for(j = 1; j < iter; j++){
            if(d(j) < e0) {
              e0 = d(j);
              col = j;
            }
          }

#ifdef DMTK_DEBUG
          dmtk::out() << setprecision(12) << "Iter = " << iter << "  Ener = " << e0 << endl;
#endif // DMTK_DEBUG
          // tqli returns ascending eigenvalues, so d(n) is the n-th target.
          bool converged = !force_maxiter && iter >= static_cast<uint>(nvectors);
          for(int n = 0; n < nvectors && n < static_cast<int>(iter); ++n) {
            if(converged && !(fabs(previous(n)-d(n)) <= tol)) converged = false;
            previous(n) = d(n);
          }
          previous_candidate_rejected = false;
          if(converged && residual_tol >= 0. && !exhausted && !breakdown) {
            // Settled Ritz energies alone do not certify their states.
            // Compute vectors only at candidate stops and check every
            // requested Ritz residual; z(column,row) holds the vectors.
            if(!ritz_vectors_ready) {
              d(Range(0,iter-1)) = a(Range(1,iter));
              e(Range(0,iter-1)) = b(Range(1,iter));
              tqli(d, e, iter, z, true);
              ritz_vectors_ready = true;
            }
            for(int n = 0; n < nvectors; ++n) {
              if(!(fabs(b(iter)*z(n,iter-1)) <= residual_tol)) {
                converged = false;
                break;
              }
            }
            previous_candidate_rejected = !converged;
          }
          if(exhausted || breakdown || converged)
          { 
             // Convergence, breakdown or the available iteration limit.
             if(fabs(b(iter)) < tol) dmtk::out() << "b(iter) < tol " << b(iter) << endl;
             dmtk::out() << setprecision(12) << "E0 = " << e0 << endl;
             maxiter = iter;
             if(nvectors > maxiter)
               throw std::runtime_error("Lanczos Krylov space has fewer vectors than requested targets");
             if(!calc_gs) {
               Vector<size_t> indices(maxiter);
               indexx<double,Vector<double> >(maxiter, d, indices, true);
               for(int n = 0; n < nvectors; n++) ener(n) = d(indices(n));
               return;
             }
             if(!ritz_vectors_ready) {
               d(Range(0,iter-1)) = a(Range(1,iter));
               e(Range(0,iter-1)) = b(Range(1,iter));
               tqli(d, e, iter, z, true);
             }
             xc = z(col, Range(0,iter));

#ifdef USE_LANCZOS_VECTORS
             break; // Exit Lanczos iteration. Re-start for ground state.
#endif
             outputfile.flush();
             outputfile.close();
         
             Vector<size_t> indx(maxiter); 
//             d(Range(0,maxiter-1)) = d(Range(1,maxiter));
             indexx<double,Vector<double> >(maxiter, d, indx, true);

#ifdef DMTK_DEBUG
             for(int n = 0; n < maxiter; n++)
               dmtk::out() << "LEVEL " << n << " " << " " << d(indx(n)) << endl;
#endif // DMTK_DEBUG

             for(int n = 0; n < nvectors; n++){
               ifstream inputfile;
               if(reorthogonalize) inputfile.exceptions(std::ios::failbit | std::ios::badbit);
               inputfile.open(vectors_file,std::ios::in|std::ios::binary);
               if(!inputfile) dmtk::err() << "*** ERROR 2: Lanczos could not open vectors.dat\n";
               gs(n) = T(0);
               ener(n) = d(indx(n));
//               cout << n << " " << ener(n) << endl;
               for(int i = 0; i < maxiter; i++){
                 read_basis(x2, inputfile);
                 axpy(T(z(indx(n),i)), x2, gs(n));
               }
               inputfile.close();

             }

             return;
          }
        } // diagonalization of tridiagonal matrix

      } 
    } // Lanczos iteration
  } // main iteration

}

template<class T, class A, class B>
void verif_hamiltonian(A& m, const B& gs)
{
  Matrix<T> h(gs.size(),gs.size());
  B aux(gs);

  dmtk::out() << "***** VERIFYING HAMILTONIAN *****\n";
  for(int col = 0; col < gs.size(); col++){
    aux = T(0);
    aux(col) = T(1);
    h.column(col) = product(m,aux);
    for(int row = 0; row < gs.size(); row++){
      dmtk::out() << col << " " << row << " " << h.column(col)(row) << endl;
    }
  }
  for(int col = 0; col < gs.size(); col++){
    for(int row = col+1; row < gs.size(); row++){
      if(fabs(h(row,col)-h(col,row)) > 1.e-5) 
        dmtk::out() << "ERROR: " << col << " " << row << " " << h(col,row) << " " << h(row,col) << endl;
    }
  }
}
 
#undef SWAP
#undef SIGN

} // namespace dmtk

#endif // __DMTK_LANCZOS_CC__
