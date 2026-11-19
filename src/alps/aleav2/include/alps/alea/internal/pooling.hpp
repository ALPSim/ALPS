/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */
#pragma once

#include <alps/alea/core.hpp>
#include <alps/alea/util.hpp>

#include <alps/alea/variance.hpp>
#include <alps/alea/covariance.hpp>
#include <alps/alea/autocorr.hpp>

#include <alps/alea/internal/joined.hpp>

#include <limits>
#include <vector>

namespace alps { namespace alea { namespace internal {

/** Returns covariance, or construct from variances if not available */
template <typename Result>
typename eigen<typename traits<Result>::cov_type>::matrix
get_cov(const Result &result)
{
    return result.cov();
}

/** Returns covariance, or construct from variances if not available */
template <typename T, typename Str>
typename eigen<typename traits<cov_result<T,Str>>::cov_type>::matrix
get_cov(const var_result<T,Str> &result)
{
    return result.var().asDiagonal();
}

/** Returns covariance, or construct from variances if not available */
template <typename T>
typename eigen<typename traits<cov_result<T>>::cov_type>::matrix
get_cov(const autocorr_result<T> &result)
{
    return result.var().asDiagonal();
}

template <typename Result, typename Derived>
using diff_scalar_type = add_scalar_type<
                            typename traits<Result>::value_type,
                            typename Eigen::internal::traits<Derived>::Scalar>;

/**
 * Difference of a result from a reference, rotated so that its components are
 * uncorrelated: `mean` and `var` are the difference and the per-observation
 * variance along each principal axis, `nmeas` the number of observations that
 * enter the F distribution of Hotelling's T^2 (one sample).
 */
template <typename T>
struct diff_stats
{
    column<T> mean;
    column<make_real_type<T>> var;
    double nmeas;
};

/**
 * Rounding level of the eigenvalues of a covariance with diagonal `var`
 * around `mean`. Exactly correlated components (duplicates, sums,
 * symmetries) make the covariance singular. Its null eigenvalues are then not
 * zero but rounding: accumulating sum(x x^T) and subtracting count mean
 * mean^T leaves errors of order eps (var + mean^2), of either sign.
 */
template <typename Mean, typename Var>
double rounding_variance(const Eigen::MatrixBase<Mean> &mean,
                         const Eigen::MatrixBase<Var> &var)
{
    if (mean.rows() == 0)
        return 0;
    return 64 * mean.rows() * std::numeric_limits<double>::epsilon()
         * (var.array() + mean.array().abs2()).maxCoeff();
}

/**
 * Drop the null directions of a rotated difference: a direction whose
 * variance is at rounding level carries no information, and is removed
 * (difference and variance set to zero, which t2_test() skips) if the
 * difference along it is rounding as well. A difference beyond rounding
 * there keeps a zero variance, so that the test rejects.
 */
template <typename T, typename Real>
void drop_null_directions(column<T> &diff, column<Real> &var,
                          double var_cut, double diff_cut)
{
    for (int i = 0; i != var.rows(); ++i) {
        if (var[i] > var_cut)
            continue;
        var[i] = 0;
        if (std::abs(diff[i]) <= diff_cut)
            diff[i] = 0;
    }
}

/**
 * Get difference of `result` from the exact `expected`
 */
template <typename Result, typename Derived,
          typename T=diff_scalar_type<Result, Derived>>
diff_stats<T> make_diff(const Result &result,
                        const Eigen::MatrixBase<Derived> &expected)
{
    if ((int)result.size() != expected.size())
        throw size_mismatch();

    diff_stats<T> diff;
    // observations() is the count consistent with var() (for autocorr_result
    // the number of bins at the binning level used), not count().
    diff.nmeas = result.observations();
    diff.mean = result.mean() - expected;

    if (traits<Result>::HAVE_COV) {
        const typename eigen<T>::matrix cov = get_cov(result);
        Eigen::SelfAdjointEigenSolver<typename eigen<T>::matrix> ecov(cov);
        diff.mean = ecov.eigenvectors().adjoint() * diff.mean;
        diff.var = ecov.eigenvalues();
        const double eps = 64 * result.size() * std::numeric_limits<double>::epsilon();
        drop_null_directions(diff.mean, diff.var,
                             rounding_variance(result.mean(), result.var()),
                             eps * (result.mean().cwiseAbs().maxCoeff()
                                    + expected.cwiseAbs().maxCoeff()));
    } else {
        // Only components with zero variance are exact; keep them out of the
        // test if they agree up to rounding.
        diff.var = result.var();
        const double eps = 64 * result.size() * std::numeric_limits<double>::epsilon();
        drop_null_directions(diff.mean, diff.var, 0,
                             eps * (result.mean().cwiseAbs().maxCoeff()
                                    + expected.cwiseAbs().maxCoeff()));
    }
    return diff;
}

/**
 * Difference of two independent results, for a test of equal means.
 *
 * The covariance of the difference of the means is S = C1/N1 + C2/N2 (Welch),
 * which does not assume equal variances: two runs usually differ in length,
 * update scheme and autocorrelation. The degrees of freedom are those of
 * Krishnamoorthy and Yu (2004), the affine-invariant form of Nel and van der
 * Merwe, which for one component is Welch-Satterthwaite:
 *
 *     nu = (r + r^2) / sum_k [tr(W_k^2) + tr(W_k)^2] / (N_k - 1),
 *     W_k = S^{-1/2} (C_k / N_k) S^{-1/2},
 *
 * with r the number of directions of S that are not null. The result is
 * returned as a one-sample difference with nu + 1 observations, so that
 * Hotelling's T^2 is referred to F(r, nu + 1 - r).
 */
template <typename Result1, typename Result2,
          typename T=joined_value_type<Result1, Result2>>
diff_stats<T> pool_var(const Result1 &r1, const Result2 &r2)
{
    if (r1.size() != r2.size())
        throw size_mismatch();

    typedef typename eigen<T>::matrix matrix;
    const double obs1 = r1.observations();
    const double obs2 = r2.observations();
    const matrix c1 = get_cov(r1) / obs1;
    const matrix c2 = get_cov(r2) / obs2;

    Eigen::SelfAdjointEigenSolver<matrix> es(c1 + c2);
    const matrix &v = es.eigenvectors();
    const auto &lambda = es.eigenvalues();

    // Whitened covariances of the two means, restricted to the directions
    // above rounding (see rounding_variance()).
    const double cut = rounding_variance(r1.mean(), r1.var()) / obs1
                     + rounding_variance(r2.mean(), r2.var()) / obs2;
    std::vector<int> keep;
    for (int i = 0; i != lambda.rows(); ++i)
        if (lambda[i] > cut)
            keep.push_back(i);

    double denom = 0;
    const matrix w1 = v.adjoint() * c1 * v, w2 = v.adjoint() * c2 * v;
    for (int k = 0; k != 2; ++k) {
        const matrix &wk = k == 0 ? w1 : w2;
        const double nk = k == 0 ? obs1 : obs2;
        double tr = 0, tr2 = 0;
        for (int i : keep) {
            tr += std::real(wk(i, i)) / lambda[i];
            for (int j : keep)
                tr2 += std::norm(wk(i, j)) / (lambda[i] * lambda[j]);
        }
        denom += (tr2 + tr * tr) / (nk - 1);
    }
    const double r = keep.size();
    // Both results exact: nothing to estimate, t2_test then only compares
    // the difference with zero.
    const double nu = keep.empty() ? obs1 + obs2 - 2 : (r + r * r) / denom;

    diff_stats<T> diff;
    diff.nmeas = nu + 1;
    diff.mean = v.adjoint() * (r1.mean() - r2.mean());
    diff.var = lambda;
    const double eps = 64 * r1.size() * std::numeric_limits<double>::epsilon();
    drop_null_directions(diff.mean, diff.var, cut,
                         eps * (r1.mean().cwiseAbs().maxCoeff()
                                + r2.mean().cwiseAbs().maxCoeff()));
    // t2_test expects per-observation variances, stderror^2 * nmeas
    diff.var *= diff.nmeas;
    return diff;
}

}}} /* namespace alps::alea::internal */
