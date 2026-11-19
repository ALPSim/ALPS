// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// A two-component observable with known covariance [[1, 0.5], [0.5, 2]]:
// x1 = z1, x2 = 0.5 z1 + sqrt(1.75) z2 for independent standard normal z.
#include "streams.hpp"

#include <alps/alea.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

namespace {

const std::size_t n = std::size_t(1) << 16;
const double sigma[2][2] = {{1.0, 0.5}, {0.5, 2.0}};

std::vector<std::vector<double> > samples() {
    statistics_test::normal xi(31337);
    std::vector<std::vector<double> > values(n, std::vector<double>(2));
    for (auto & value : values) {
        const double z1 = xi(), z2 = xi();
        value[0] = z1;
        value[1] = 0.5 * z1 + std::sqrt(1.75) * z2;
    }
    return values;
}

// A variance estimated from n samples scatters by sqrt(2 / n) = 0.55%, the
// error by half of that: 1.5% is more than four standard deviations.
void expect_errors(double e0, double e1, const char * library) {
    EXPECT_NEAR(std::sqrt(sigma[0][0] / n), e0, 0.015 * std::sqrt(sigma[0][0] / n)) << library;
    EXPECT_NEAR(std::sqrt(sigma[1][1] / n), e1, 0.015 * std::sqrt(sigma[1][1] / n)) << library;
}

}

TEST(VectorCovariance, ComponentErrors)
{
    alps::alea::var_acc<double> v2(2);
    for (const std::vector<double> & value : samples())
        v2 << value;
    const alps::alea::var_result<double> result = v2.finalize();
    expect_errors(result.stderror()[0], result.stderror()[1], "ALEA v2");
}

// The full covariance. An entry estimated from n
// samples scatters by sqrt((S_ii S_jj + S_ij^2) / n), at most 0.8%.
TEST(VectorCovariance, AleaCovariance)
{
    alps::alea::cov_acc<double> cov(2);
    for (const std::vector<double> & value : samples())
        cov << value;
    const alps::alea::cov_result<double> result = cov.finalize();
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            EXPECT_NEAR(sigma[i][j], result.cov()(i, j), 0.035 * std::sqrt(sigma[i][i] * sigma[j][j]))
                << i << ", " << j;
}
