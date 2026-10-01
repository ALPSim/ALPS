// SPDX-License-Identifier: MIT
#include <alps/config.h>
#include <boost/numeric/bindings/detail/config/fortran.hpp>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>

extern "C" {
void FORTRAN_ID(ilaver)(void*, void*, void*);
void FORTRAN_ID(dgesv)(const fortran_int_t*, const fortran_int_t*, double*,
                      const fortran_int_t*, fortran_int_t*, double*,
                      const fortran_int_t*, fortran_int_t*);
void FORTRAN_ID(dgemm)(const char*, const char*, const fortran_int_t*,
                      const fortran_int_t*, const fortran_int_t*, const double*,
                      const double*, const fortran_int_t*, const double*,
                      const fortran_int_t*, const double*, double*, const fortran_int_t*);
void FORTRAN_ID(spotrf)(const char*, const fortran_int_t*, float*,
                       const fortran_int_t*, fortran_int_t*);
}

int main() {
    // ILAVER only writes integers. Oversized, aligned buffers safely detect both
    // 32- and 64-bit providers before any routine can misread a dimension.
    alignas(std::int64_t) std::array<unsigned char, 8> major, minor, patch;
    major.fill(0xff); minor.fill(0xff); patch.fill(0xff);
    FORTRAN_ID(ilaver)(major.data(), minor.data(), patch.data());
    for (const auto& bytes : {major, minor, patch}) {
        fortran_int_t value;
        std::memcpy(&value, bytes.data(), sizeof(value));
        if (value < 0 || value > 1000)
            throw std::runtime_error("LAPACK integer ABI differs from ALPS");
        for (std::size_t i = sizeof(value); i < bytes.size(); ++i)
            if (bytes[i] != 0xff)
                throw std::runtime_error("LAPACK integer writes exceed the ALPS integer ABI");
    }

    const fortran_int_t n = 2, nrhs = 1;
    double a[] = {3., 1., 1., 2.}; // column-major, symmetric positive definite
    double b[] = {9., 8.};
    fortran_int_t pivot[2] = {}, info = -1;
    FORTRAN_ID(dgesv)(&n, &nrhs, a, &n, pivot, b, &n, &info);
    // Negated <= comparisons also reject NaN results.
    if (info != 0 || !(std::abs(b[0] - 2.) <= 1e-12) || !(std::abs(b[1] - 3.) <= 1e-12))
        throw std::runtime_error("LAPACK solve failed");

    const double matrix[] = {3., 1., 1., 2.}, one = 1., zero = 0.;
    double result[2] = {};
    const char normal = 'N';
    FORTRAN_ID(dgemm)(&normal, &normal, &n, &nrhs, &n, &one,
                     matrix, &n, b, &n, &zero, result, &n);
    if (!(std::abs(result[0] - 9.) <= 1e-12) || !(std::abs(result[1] - 8.) <= 1e-12))
        throw std::runtime_error("BLAS matrix multiplication failed");

    // A mixed f2c/native BLAS ABI can return success but the wrong factor:
    // CLAPACK expects SDOT to return double; native BLAS returns float.
    // For A = [[4,2],[2,3]], L must be [[2,0],[1,sqrt(2)]].
    float cholesky[] = {4.f, 2.f, 2.f, 3.f};
    const char lower = 'L';
    FORTRAN_ID(spotrf)(&lower, &n, cholesky, &n, &info);
    if (info != 0 || !(std::abs(cholesky[0] - 2.f) <= 1e-6f)
        || !(std::abs(cholesky[1] - 1.f) <= 1e-6f)
        || !(std::abs(cholesky[3] - std::sqrt(2.f)) <= 1e-6f))
        throw std::runtime_error("LAPACK single-precision Cholesky failed; check BLAS/LAPACK ABI compatibility");
}
