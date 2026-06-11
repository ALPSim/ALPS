/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2026 by the ALPS collaboration
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

// Layout test for alps::numeric::kron().
//
// Standard Kronecker product layout:
//   (A (x) B)(i1*rows(B) + i2, j1*cols(B) + j2) = A(i1,j1) * B(i2,j2)
// i.e. the INNER factor B's dimensions stride the output. The historical
// implementation used A's dimensions in the index formula, which coincides
// with the standard layout only for equal-size operands; for unequal shapes
// it scattered entries (and indexed out of bounds when A is the larger
// factor). The equal-size arm below doubles as a no-behavior-change check
// for callers passing equal-size square operands.

#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>

#include <cstdio>
#include <cstddef>

namespace {

int failures = 0;

// deterministic, position-dependent fill so misplaced entries can't
// accidentally match
double fill(std::size_t i, std::size_t j, double salt)
{
    return salt + 10.0 * static_cast<double>(i) + static_cast<double>(j);
}

void run_case(std::size_t ar, std::size_t ac, std::size_t br, std::size_t bc,
              const char* label)
{
    typedef alps::numeric::matrix<double> matrix_type;
    matrix_type A(ar, ac), B(br, bc);
    for (std::size_t i = 0; i < ar; ++i)
        for (std::size_t j = 0; j < ac; ++j) A(i, j) = fill(i, j, 1.25);
    for (std::size_t i = 0; i < br; ++i)
        for (std::size_t j = 0; j < bc; ++j) B(i, j) = fill(i, j, -3.5);

    matrix_type MM = alps::numeric::kron(A, B);

    if (num_rows(MM) != ar * br || num_cols(MM) != ac * bc) {
        std::fprintf(stderr,
                     "FAIL %s: result is %lux%lu, expected %lux%lu\n", label,
                     static_cast<unsigned long>(num_rows(MM)),
                     static_cast<unsigned long>(num_cols(MM)),
                     static_cast<unsigned long>(ar * br),
                     static_cast<unsigned long>(ac * bc));
        ++failures;
        return;
    }

    for (std::size_t i1 = 0; i1 < ar; ++i1)
        for (std::size_t j1 = 0; j1 < ac; ++j1)
            for (std::size_t i2 = 0; i2 < br; ++i2)
                for (std::size_t j2 = 0; j2 < bc; ++j2) {
                    const double expected = A(i1, j1) * B(i2, j2);
                    const double actual = MM(i1 * br + i2, j1 * bc + j2);
                    if (expected != actual) {
                        std::fprintf(stderr,
                                     "FAIL %s: (A kron B)(%lu,%lu) = %.17g, "
                                     "expected A(%lu,%lu)*B(%lu,%lu) = %.17g\n",
                                     label,
                                     static_cast<unsigned long>(i1 * br + i2),
                                     static_cast<unsigned long>(j1 * bc + j2),
                                     actual,
                                     static_cast<unsigned long>(i1),
                                     static_cast<unsigned long>(j1),
                                     static_cast<unsigned long>(i2),
                                     static_cast<unsigned long>(j2),
                                     expected);
                        ++failures;
                        return; // one report per case is enough
                    }
                }
}

} // namespace

int main()
{
    // Equal square operands: old and new index formulas coincide here, so
    // this arm characterizes that existing equal-size callers see no change.
    run_case(3, 3, 3, 3, "kron 3x3 (x) 3x3 (equal-size, current callers)");

    // Unequal shapes: pin the standard Kronecker layout (inner factor B
    // strides the output). These fail on the historical A-dims formula.
    run_case(2, 3, 4, 2, "kron 2x3 (x) 4x2");
    run_case(4, 2, 2, 5, "kron 4x2 (x) 2x5");
    run_case(1, 1, 3, 2, "kron 1x1 (x) 3x2");

    if (failures) {
        std::fprintf(stderr, "matrix_kron: %d failure(s)\n", failures);
        return 1;
    }
    std::printf("matrix_kron: PASS\n");
    return 0;
}
