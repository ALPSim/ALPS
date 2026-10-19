/*****************************************************************************
*
* ALPS Project Applications
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Helpers for SAVE_MPS: assembling the final DMRG state as a matrix product
// state. DMRGTask<T>::extract_mps (dmrg.h) reads the block transformations
// and the two-site wave function that the final sweep leaves in the
// temporary files and uses split_two_site below for the center.

#ifndef ALPS_DMRG_MPS_EXPORT_H
#define ALPS_DMRG_MPS_EXPORT_H

#include "dmtk/lapack_interface.h"

#include <algorithm>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

extern "C" {

void FORTRAN_ID(dgesvd)(
    const char &jobu, const char &jobvt, const int &m, const int &n,
    double *a, const int &lda, double *s,
    double *u, const int &ldu, double *vt, const int &ldvt,
    double *work, const int &lwork, int &info);

void FORTRAN_ID(zgesvd)(
    const char &jobu, const char &jobvt, const int &m, const int &n,
    std::complex<double> *a, const int &lda, double *s,
    std::complex<double> *u, const int &ldu, std::complex<double> *vt, const int &ldvt,
    std::complex<double> *work, const int &lwork, double *rwork, int &info);

}

namespace dmrg_mps {

// Thin SVD of the column-major m x n matrix a (overwritten): a = u diag(s) vt
inline void gesvd(int m, int n, std::vector<double> &a, std::vector<double> &s,
                  std::vector<double> &u, std::vector<double> &vt)
{
  int k = std::min(m, n), info = 0;
  s.resize(k); u.resize(std::size_t(m) * k); vt.resize(std::size_t(k) * n);
  double query;
  FORTRAN_ID(dgesvd)('S', 'S', m, n, &a[0], m, &s[0], &u[0], m, &vt[0], k, &query, -1, info);
  std::vector<double> work(static_cast<std::size_t>(query));
  FORTRAN_ID(dgesvd)('S', 'S', m, n, &a[0], m, &s[0], &u[0], m, &vt[0], k, &work[0], int(work.size()), info);
  if (info != 0)
    throw std::runtime_error("SAVE_MPS: dgesvd failed to split the center of the MPS");
}

inline void gesvd(int m, int n, std::vector<std::complex<double> > &a, std::vector<double> &s,
                  std::vector<std::complex<double> > &u, std::vector<std::complex<double> > &vt)
{
  int k = std::min(m, n), info = 0;
  s.resize(k); u.resize(std::size_t(m) * k); vt.resize(std::size_t(k) * n);
  std::vector<double> rwork(5 * std::size_t(k));
  std::complex<double> query;
  FORTRAN_ID(zgesvd)('S', 'S', m, n, &a[0], m, &s[0], &u[0], m, &vt[0], k, &query, -1, &rwork[0], info);
  std::vector<std::complex<double> > work(static_cast<std::size_t>(query.real()));
  FORTRAN_ID(zgesvd)('S', 'S', m, n, &a[0], m, &s[0], &u[0], m, &vt[0], k, &work[0], int(work.size()), &rwork[0], info);
  if (info != 0)
    throw std::runtime_error("SAVE_MPS: zgesvd failed to split the center of the MPS");
}

// Split the row-major two-site tensor psi[a][s][t][b] exactly into
// left[a][s][x] (orthonormal columns) and right[x][t][b] = S Vh, dropping
// only singular values below round-off (1e-14 of the largest).
template <class T>
void split_two_site(const std::vector<T> &psi, std::size_t da, std::size_t ds,
                    std::size_t dt, std::size_t db,
                    std::vector<T> &left, std::vector<std::size_t> &left_shape,
                    std::vector<T> &right, std::vector<std::size_t> &right_shape)
{
  int m = int(da * ds), n = int(dt * db);
  std::vector<T> a(std::size_t(m) * n);   // column-major copy of psi as (a s) x (t b)
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      a[i + std::size_t(j) * m] = psi[std::size_t(i) * n + j];

  std::vector<double> s;
  std::vector<T> u, vt;
  gesvd(m, n, a, s, u, vt);
  int k = int(s.size());
  std::size_t keep = 1;
  while (keep < s.size() && s[keep] > 1e-14 * s[0])
    ++keep;

  left.assign(std::size_t(m) * keep, T(0));
  for (int i = 0; i < m; ++i)
    for (std::size_t x = 0; x < keep; ++x)
      left[i * keep + x] = u[i + x * std::size_t(m)];
  right.assign(keep * std::size_t(n), T(0));
  for (std::size_t x = 0; x < keep; ++x)
    for (int j = 0; j < n; ++j)
      right[x * n + j] = s[x] * vt[x + std::size_t(j) * k];

  left_shape.assign(3, 0);
  left_shape[0] = da; left_shape[1] = ds; left_shape[2] = keep;
  right_shape.assign(3, 0);
  right_shape[0] = keep; right_shape[1] = dt; right_shape[2] = db;
}

} // namespace dmrg_mps

#endif // ALPS_DMRG_MPS_EXPORT_H
