// SPDX-License-Identifier: MIT
#include <alps/numeric/diagonal_matrix.hpp>
#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>
#include <alps/numeric/matrix/gemm.hpp>
#include <cmath>
#include <complex>
#include <stdexcept>

int main() {
    // Exercise the installed adapters and BLAS/LAPACK linkage through ALPS APIs.
    alps::numeric::matrix<double> a(2, 2), square(2, 2);
    a(0, 0) = 3.; a(0, 1) = 1.;
    a(1, 0) = 1.; a(1, 1) = 2.;
    alps::numeric::gemm(a, a, square);
    if (!(std::abs(square(0, 0) - 10.) < 1e-12)
        || !(std::abs(square(0, 1) - 5.) < 1e-12)
        || !(std::abs(square(1, 0) - 5.) < 1e-12)
        || !(std::abs(square(1, 1) - 5.) < 1e-12))
        throw std::runtime_error("Installed SDK matrix multiplication failed");

    using ComplexMatrix = alps::numeric::matrix<std::complex<double>>;
    ComplexMatrix h(2, 2), eigenvectors;
    h(0, 0) = 2.; h(0, 1) = {0., 1.};
    h(1, 0) = {0., -1.}; h(1, 1) = 2.;
    alps::numeric::associated_real_vector<ComplexMatrix>::type eigenvalues(2);
    alps::numeric::heev(h, eigenvectors, eigenvalues);
    if (!(std::abs(eigenvalues[0] - 3.) < 1e-12)
        || !(std::abs(eigenvalues[1] - 1.) < 1e-12))
        throw std::runtime_error("Installed SDK Hermitian eigenvalues failed");
    for (int column = 0; column < 2; ++column)
        for (int row = 0; row < 2; ++row) {
            const auto residual = h(row, 0) * eigenvectors(0, column)
                + h(row, 1) * eigenvectors(1, column)
                - eigenvalues[column] * eigenvectors(row, column);
            if (!(std::abs(residual) < 1e-12))
                throw std::runtime_error("Installed SDK Hermitian eigenvectors failed");
        }
}
