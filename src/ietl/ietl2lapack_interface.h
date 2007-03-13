/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2003 by Prakash Dayal <prakash@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>,
*                            Rene Villiger <rvilliger@smile.ch>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

/* $Id: ietl2lapack_interface.h,v 1.5 2003/09/05 08:12:38 troyer Exp $ */

#ifndef IETL_LAPACK_INTERFACE_H
#define IETL_LAPACK_INTERFACE_H
#include <complex>

#if defined(_AIX) || defined(__blrts__) || defined(hpux) || defined(__IBMCPP__)
#define MTL_FCALL(x) x
#else
#define MTL_FCALL(x) x##_
#endif

extern "C" {
void MTL_FCALL(sstev)(const char& jobz, const int& n, float sd[],
            float se[], float sz[], const int& ldz, 
            float swork[], int& info);

void MTL_FCALL(dstev)(const char& jobz, const int& n, double dd[],
            double de[], double dz[], const int& ldz, 
            double swork[], int& info);
  
void MTL_FCALL(csteqr)(const char& jobz, const int& n, float sd[],
            float se[], std::complex<float> sz[], const int& ldz, 
            float swork[], int& info);

void MTL_FCALL(zsteqr)(const char& jobz, const int& n, double dd[],
            double de[], std::complex<double> dz[], const int& ldz, 
            double dwork[], int& info);  

void MTL_FCALL(dsyevx)(const char& jobz,            const char& range,         const char& uplo,
                             const int& n,                double da[],               const int& lda,
                             const double& vl,            const double& vu,          const int& il,
                             const int& iu,               const double& abstol,      int& m,
                             double w[],                  double z[],                const int& ldz,
                             double work[],               const int& lwork,          int iwork[],
                             int ifail[],                 int& info);                                             

void MTL_FCALL(ssyevx)(const char& jobz,            const char& range,         const char& uplo,
                             const int& n,                float da[],                const int& lda,
                             const float& vl,             const float& vu,           const int& il,
                             const int& iu,               const float& abstol,       int& m,
                             float w[],                   float z[],                 const int& ldz,
                             float work[],                const int& lwork,          int iwork[],
                             int ifail[],                 int& info);                                             
                             
void MTL_FCALL(zheevx)(const char& jobz,            const char& range,         const char& uplo,
                             const int& n,                std::complex<double> da[], const int& lda,
                             const double& vl,            const double& vu,          const int& il,
                             const int& iu,               const double& abstol,      int& m,
                             double w[],                  std::complex<double> z[],  const int& ldz,
                             std::complex<double> work[], const int& lwork,          double rwork[],
                             int iwork[],                 int ifail[],               int& info);                  

void MTL_FCALL(cheevx)(const char& jobz,            const char& range,         const char& uplo,
                             const int& n,                std::complex<float> da[],  const int& lda,
                             const float& vl,             const float& vu,           const int& il,
                             const int& iu,               const float& abstol,       int& m,
                             float w[],                   std::complex<float> z[],   const int& ldz,
                             std::complex<float> work[],  const int& lwork,          float rwork[],
                             int iwork[],                 int ifail[],               int& info);                  
                             
void MTL_FCALL(dsysv)(const char& uplo,             const int& n,              const int& nrhs,
                            double a[],                   const int& lda,            int ipiv[],
                            double b[],                   const int& ldb,            double work[],
                            const int& lwork,             int& info);                                             
                            
void MTL_FCALL(ssysv)(const char& uplo,             const int& n,              const int& nrhs,
                            float a[],                    const int& lda,            int ipiv[],
                            float b[],                    const int& ldb,            float work[],
                            const int& lwork,             int& info);                                             
                        
void MTL_FCALL(chesv)(const char& uplo,             const int& n,              const int& nrhs,
                            std::complex<float> a[],      const int& lda,            int ipiv[],
                            std::complex<float> b[],      const int& ldb,            std::complex<float> work[],
                            const int& lwork,             int& info);                                             
                            
void MTL_FCALL(zhesv)(const char& uplo,             const int& n,              const int& nrhs,
                            std::complex<double> a[],     const int& lda,            int ipiv[],
                            std::complex<double> b[],     const int& ldb,            std::complex<double> work[],
                            const int& lwork,             int& info);                                             

void MTL_FCALL(dsyev)(const char& jobz,            const char& uplo, const int& n,
                            double da[],                 const int& lda,   double w[],
                            double work[],               const int& lwork, int& info);                 

void MTL_FCALL(ssyev)(const char& jobz,            const char& uplo, const int& n,
                            float da[],                  const int& lda,   float w[],
                            float work[],                const int& lwork, int& info);                  

void MTL_FCALL(zheev)(const char& jobz,            const char& uplo, const int& n,
                            std::complex<double> da[],   const int& lda,   double w[],
                            std::complex<double> work[], const int& lwork, double rwork[], int& info); 

void MTL_FCALL(cheev)(const char& jobz,            const char& uplo, const int& n,
                            std::complex<float> da[],    const int& lda,   float w[],
                            std::complex<float> work[],  const int& lwork, float rwork[],  int& info);

}
#endif
