/*
 * Copyright (C) 1998-2018 ALPS Collaboration.
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

// Adapted from CQMP/Maxent src/eigen_lapack.hpp at
// a3992534228b0de732e033ae23154969724814be; see ../README.md.
// Updated 2026-10-14 to size thin factors and report LAPACK errors.
#pragma once

#include <algorithm>
#include <stdexcept>
#include <string>

#ifdef HAVE_LAPACK
extern "C" void dgesvd_( const char* jobu, const char* jobvt,
                                              const int* m, const int* n, double* a, const int* lda,
                                              double* s, double* u, const int* ldu,
                                              double* vt, const int* ldvt,
                                              double* work, const int* lwork, int* info);

///performs a SVD on input matrix K
/// returns K =  U (S) V^T where S is a vector of the diagonal values of
/// matrix Sigma
void lapack_svd(matrix_type &K, vector_type &S, matrix_type &Vt, matrix_type &U){

   /*boost::numeric::bindings::lapack::gesvd('S', 'S', Kt, S, U_, Vt_);*/
  //use 'S' for thin U,Vt matrices
  char jobu = 'S';
  char jobvt = 'S';
  int m = K.rows();
  int n = K.cols();
  // Thin SVD writes exactly min(m,n) modes, including when NDAT > NFREQ.
  const int rank = std::min(m, n);
  S.resize(rank);
  U.resize(m, rank);
  Vt.resize(rank, n);
  double *a = K.data();
  int lda = K.outerStride();
  int ldu = U.outerStride();
  int ldvt = Vt.outerStride();

  double *vt = Vt.data();
  double *u = U.data();
  double *s =S.data();

  double dummywork;
  int LWORK=-1;
  int INFO=0;
  //get best workspace
  dgesvd_(&jobu,&jobvt,&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,&dummywork,&LWORK,&INFO);
  if (INFO != 0)
    throw std::runtime_error("LAPACK SVD workspace query failed: " + std::to_string(INFO));

  LWORK=int(dummywork)+32;
  vector_type WORK(LWORK);

  dgesvd_(&jobu,&jobvt,&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,WORK.data(),&LWORK,&INFO);
  if(INFO!=0)
    throw std::runtime_error("LAPACK SVD failed: " + std::to_string(INFO));



}
#endif
