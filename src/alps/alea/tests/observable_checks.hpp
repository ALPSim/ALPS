// SPDX-License-Identifier: MIT
#pragma once

#include <gtest/gtest.h>
#include <initializer_list>
#include <type_traits>

namespace alps_test {

// Reference values come from the pre-migration stdout fixtures. Bounds match
// their printed precision, so formatting is no longer part of the numeric
// contract. In particular, these are regression checks, not statistical tests.
struct Estimate {
    double mean;
    double error;
    double mean_tolerance;
    double error_tolerance;
};

inline void expect_estimate_values(double mean, double error, const Estimate& reference)
{
    EXPECT_NEAR(mean, reference.mean, reference.mean_tolerance);
    EXPECT_NEAR(error, reference.error, reference.error_tolerance);
}

template<class Observable>
void expect_estimate(const Observable& observable,
                     std::initializer_list<Estimate> references, const char* expression)
{
    SCOPED_TRACE(expression);
    const auto mean = observable.mean();
    const auto error = observable.error();
    if constexpr (std::is_arithmetic_v<decltype(mean)>) {
        ASSERT_EQ(references.size(), 1u);
        expect_estimate_values(mean, error, *references.begin());
    } else {
        ASSERT_EQ(mean.size(), references.size());
        ASSERT_EQ(error.size(), references.size());
        std::size_t i = 0;
        for (const auto& reference : references) {
            SCOPED_TRACE(::testing::Message() << "component " << i);
            expect_estimate_values(mean[i], error[i], reference);
            ++i;
        }
    }
}

template<class Observable>
void expect_uniform_estimate(const Observable& observable, const Estimate& reference,
                             const char* expression)
{
    SCOPED_TRACE(expression);
    const auto mean = observable.mean();
    const auto error = observable.error();
    ASSERT_EQ(mean.size(), 10u);
    ASSERT_EQ(error.size(), 10u);
    for (std::size_t i = 0; i < mean.size(); ++i) {
        SCOPED_TRACE(::testing::Message() << "component " << i);
        expect_estimate_values(mean[i], error[i], reference);
    }
}
} // namespace alps_test
