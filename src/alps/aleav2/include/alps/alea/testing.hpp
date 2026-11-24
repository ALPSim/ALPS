/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */

#pragma once

#include <alps/alea/core.hpp>
#include <alps/alea/util.hpp>
#include <alps/alea/variance.hpp>
#include <alps/alea/covariance.hpp>
#include <alps/alea/fisher_f.hpp>

#include <alps/alea/internal/joined.hpp>
#include <alps/alea/internal/pooling.hpp>

#include <algorithm>
#include <cmath>

namespace alps { namespace alea {
    class t2_result;
    template <typename T> struct diag_diffs;

    template <typename R1, typename R2, typename T>
    var_result<T> pooled_var(const R1 &, const R2 &);

    template <typename R1, typename R2, typename T>
    cov_result<T> pooled_cov(const R1 &, const R2 &);
}}

namespace alps { namespace alea {

/**
 * Result of Hotelling's T^2 test
 *
 * @see test_mean(), t2_test()
 */
class t2_result
{
public:
    typedef fisher_f_distribution dist_type;

public:
    /** Initializes t2 test */
    t2_result(double score, double f_size, double f_dof)
        : score_(score), dist_(f_size, f_dof)
    { }

    /** Returns whether the lower alternate hypothesis can be tested.
     *
     * For low values of size, the F distribution does not have a mode,
     * because there is no way to distinguish the lower alternate hypothesis
     * (too large error bars) from an accidental "hit".
     */
    bool has_plower() const { return dist_.degrees_of_freedom1() > 3; }

    /** p-value in favour of the lower alternate Hypothesis */
    double pvalue_lower() const { return has_plower() ? dist_.cdf(score_) : 1; }

    /** p-value in favour of the upper alternate Hypothesis */
    double pvalue_upper() const { return tested() ? dist_.ccdf(score_) : 1; }

    /**
     * Two-sided p-value, i.e., p-value in favour of all possible alternate hypotheses.
     */
    double pvalue() const
    {
        // std::min below would turn a NaN into 1
        if (std::isnan(score_))
            return NAN;
        if (!has_plower())
            return pvalue_upper();
        return std::min(1.0, 2 * std::min(pvalue_lower(), pvalue_upper()));
    }

    /** Whether any component was tested */
    bool tested() const { return dist_.degrees_of_freedom1() != 0; }

    /** t2 statistic and argument to `dist()` */
    double score() const { return score_; }

    /** distribution */
    const dist_type &dist() const { return dist_; }

private:
    double score_;
    dist_type dist_;
};

/**
 * Perform Hotelling's T^2 test given a set of uncorrelated differences.
 *
 * Takes the difference `diff` between the two results under test, which must
 * be independently Gaussian distributed with variances `var`, and the number
 * of uncorrelated observations `nmeas`, and performs Hotelling's T^2 test on
 * it. `pools` is the number of means estimated from the data (1 for one
 * sample, and for the Welch difference of two samples, which carries its
 * degrees of freedom in `nmeas`; 2 for two samples with one pooled
 * variance), and `atol` is an upper tolerance for exact (zero variance)
 * results, compared with `var` and the squared `diff`. If no component is
 * left, the result has score 0 and pvalue() 1. test_mean() handles exact
 * components itself (see its `rtol` and `atol`) and calls this with 0.
 *
 * @see test_mean()
 */
template <typename T>
t2_result t2_test(const column<T> &diff,
                  const column<typename make_real<T>::type> &var,
                  double nmeas, size_t pools=1, double atol=1e-14
                  );

/**
 * Test mean of stochastic result `result` against known result `expected`.
 *
 * Exact components are those that the measurement computes exactly, such as
 * constants and components related by a symmetry or a sum rule. They have
 * no (or no independent) variance, so the test leaves them out where they
 * agree with `expected`, and rejects where they do not. Exact means within
 * `atol + rtol |x|`, in the units of the data, with |x| = sqrt(var + mean^2)
 * of the component: `rtol` is the relative accuracy of the measurement routine
 * (a few eps if it is careful, larger if, e.g., symmetric components are
 * computed along different paths), and `atol` a floor for components that
 * are zero, such as an expected value of zero computed with rounding. The
 * rounding of the accumulation is added, and the variance is resolved only to
 * its rounding.
 */
template <typename Result, typename Derived>
t2_result test_mean(const Result &result,
                    const Eigen::MatrixBase<Derived> &expected,
                    double rtol=internal::default_rtol, double atol=0)
{
    static_assert(is_alea_result<Result>::value, "Result is not alea result");
    static_assert(traits<Result>::HAVE_VAR, "Result1 must have variance");

    using diff_scalar = internal::diff_scalar_type<Result, Derived>;
    internal::diff_stats<diff_scalar> diff = internal::make_diff(result, expected, rtol, atol);
    return t2_test(diff.mean, diff.var, diff.nmeas, 1, 0);
}

/**
 * Test mean of stochastic result `result` against known result `expected`.
 */
template <typename Result, typename Derived>
t2_result test_mean(const Eigen::MatrixBase<Derived> &expected,
                    const Result &result,
                    double rtol=internal::default_rtol, double atol=0)
{
    return test_mean(result, expected, rtol, atol);
}

/**
 * Test mean of two independent stochastic results against each other.
 *
 * The results may have different variances and numbers of observations
 * (Welch test, see internal::pool_var()). `rtol` and `atol` apply to each
 * result, as in test_mean() against a known result.
 */
template <typename Result1, typename Result2>
typename std::enable_if<
    is_alea_result<Result1>::value && is_alea_result<Result2>::value,
    t2_result>::type
test_mean(const Result1 &result1, const Result2 &result2,
          double rtol=internal::default_rtol, double atol=0)
{
    static_assert(traits<Result1>::HAVE_VAR, "Result1 must have variance");
    static_assert(traits<Result2>::HAVE_VAR, "Result2 must have variance");

    // pool_var folds the two samples into one difference with Welch degrees
    // of freedom, so it is tested as one sample.
    using diff_scalar = internal::joined_value_type<Result1, Result2>;
    internal::diff_stats<diff_scalar> diff = internal::pool_var(result1, result2, rtol, atol);
    return t2_test(diff.mean, diff.var, diff.nmeas, 1, 0);
}


}} /* namespace alps::alea */
