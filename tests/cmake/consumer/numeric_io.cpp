// SPDX-License-Identifier: MIT
#include <alps/hdf5/matrix.hpp>
#include <alps/hdf5/numeric_vector.hpp>
#include <complex>
#include <cstdio>
#include <stdexcept>

int main() {
    using Complex = std::complex<double>;
    alps::numeric::matrix<Complex> source(2, 3);
    for (std::size_t row = 0; row != 2; ++row)
        for (std::size_t col = 0; col != 3; ++col)
            source(row, col) = Complex(10 * row + col, row + col + 0.5);
    alps::numeric::vector<Complex> vector(2);
    vector[0] = {1., -2.}; vector[1] = {3., 4.};
    char const* filename = "numeric_io_contract.h5";
    {
        alps::hdf5::archive archive(filename, "w");
        archive["matrix"] << source;
        archive["vector"] << vector;
    }
    {
        alps::hdf5::archive archive(filename, "r");
        alps::numeric::matrix<Complex> result;
        alps::numeric::vector<Complex> restored;
        archive["matrix"] >> result;
        archive["vector"] >> restored;
        if (result.num_rows() != 2 || result.num_cols() != 3 || restored.size() != 2
            || !archive.is_complex("matrix") || !archive.is_complex("vector"))
            throw std::runtime_error("Numeric archive shape/type changed");
        for (std::size_t row = 0; row != 2; ++row)
            for (std::size_t col = 0; col != 3; ++col)
                if (result(row, col) != source(row, col))
                    throw std::runtime_error("Numeric archive matrix value changed");
        if (restored[0] != vector[0] || restored[1] != vector[1])
            throw std::runtime_error("Numeric archive vector value changed");
    }
    std::remove(filename);
}
