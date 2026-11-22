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
 * Difference of a result from a reference, in coordinates where its
 * components are uncorrelated: `mean` and `var` are the difference and the
 * per-observation variance along each principal axis of the covariance
 * scaled per component (see whiten()), `nmeas` the number of observations
 * that enter the F distribution of Hotelling's T^2 (one sample). T^2 does not
 * depend on the scaling.
 */
template <typename T>
struct diff_stats
{
    column<T> mean;
    column<make_real_type<T>> var;
    double nmeas;
};

/**
 * Bound on the relative rounding error of a mean, variance or covariance
 * accumulated over `count` samples. The sums are accumulated one sample (or
 * batch) after the other, so the error of a mean is bounded by
 * count * eps/2 * mean(|x|) and that of a (co)variance by about
 * 3 count * eps/2 (var + mean^2), in units of s = sqrt(var + mean^2) per
 * component. The bound is reached in practice: a constant such as 0.1 is
 * rounded the same way at every addition, and its mean and variance are off
 * by about 0.1 count eps. This returns 2 count eps.
 */
inline double accumulation_error(double count)
{
    return 2 * (count + 4) * std::numeric_limits<double>::epsilon();
}

/**
 * Covariance and difference in principal axes, scaled per component.
 *
 * Exactly correlated components (duplicates, sums, symmetries) make a
 * covariance singular, and components that are constant (or fluctuate below
 * the precision of the accumulation) make its diagonal vanish. The computed
 * null eigenvalues are then not zero but rounding, of either sign, and a
 * difference along them is rounding as well. To tell rounding from data, the
 * covariance is first scaled to R = D^-1 C D^-1 with D = diag(scale), where
 * scale is sqrt(var + mean^2) of each component. The rounding of each entry
 * of R is then below `rel_error` whatever the units of the components,
 * and that of an eigenvalue below n rel_error. T^2 does not change, since
 * d^T C^-1 d = (D^-1 d)^T R^-1 (D^-1 d).
 *
 * Along each principal axis, with eigenvalue lambda and difference delta:
 *  - lambda at rounding level (below its bound) and delta within the bound
 *    of its rounding (`diff_error`, scaled and rotated): no information,
 *    the axis is dropped (var and diff set to zero, which t2_test() skips);
 *  - otherwise the axis is tested with the variance lambda plus the bound of
 *    its rounding. Rounding can thus not make a variance too small, and a
 *    difference beyond rounding along a null axis is tested against the
 *    rounding level instead of zero: it rejects if it is larger than the
 *    data can resolve, and does not if it is noise of a variance below the
 *    resolution.
 * NaN in the data is kept and makes the test NaN. `keep` receives the axes
 * whose variance is above rounding.
 */
template <typename T, typename Real>
void whiten(const typename eigen<T>::matrix &cov, bool diagonal,
            column<T> &diff, column<Real> &var, column<Real> diff_error,
            column<Real> scale, double rel_error,
            typename eigen<T>::matrix &vecs, std::vector<int> &keep)
{
    const int n = diff.rows();
    const double eps = std::numeric_limits<double>::epsilon();
    keep.clear();
    if (n == 0) {
        var.resize(0);
        vecs.resize(0, 0);
        return;
    }
    for (int i = 0; i != n; ++i) {
        if (scale[i] == 0)
            scale[i] = 1;       // component identically zero
        diff[i] /= scale[i];
        diff_error[i] /= scale[i];
    }
    const typename eigen<T>::matrix r =
        scale.cwiseInverse().asDiagonal() * cov * scale.cwiseInverse().asDiagonal();

    double var_error;
    if (diagonal) {
        vecs = eigen<T>::matrix::Identity(n, n);
        var = r.diagonal().real();
        var_error = rel_error;
    } else {
        Eigen::SelfAdjointEigenSolver<typename eigen<T>::matrix> es(r);
        vecs = es.eigenvectors();
        var = es.eigenvalues();
        // rounding of the entries plus that of the eigensolver
        var_error = n * (rel_error + 8 * n * eps);
    }

    // Rounding of the difference along each axis: that of the components,
    // and of the rotation itself.
    const column<T> scaled = diff;
    diff = vecs.adjoint() * scaled;
    for (int k = 0; k != n; ++k) {
        double bound = 0;
        for (int i = 0; i != n; ++i)
            bound += std::abs(vecs(i, k)) * (diff_error[i] + n * eps * std::abs(scaled[i]));

        if (var[k] <= var_error && std::abs(diff[k]) <= bound) {
            var[k] = 0;
            diff[k] = 0;
            continue;
        }
        if (!(var[k] <= var_error))
            keep.push_back(k);      // also NaN, so that it propagates
        if (var[k] < 0)
            var[k] = 0;
        var[k] += var_error;
    }
}

/**
 * Components with variance and squared difference both below `atol` (in the
 * units of the data) are taken as exact and agreeing: their difference and
 * covariance are set to zero, so that whiten() drops them.
 */
template <typename T, typename Real>
void apply_atol(typename eigen<T>::matrix &cov, column<T> &diff,
                const column<Real> &var, double atol)
{
    if (!(atol > 0))
        return;
    for (int i = 0; i != diff.rows(); ++i) {
        if (std::abs(var[i]) <= atol && std::norm(diff[i]) <= atol) {
            diff[i] = 0;
            cov.row(i).setZero();
            cov.col(i).setZero();
        }
    }
}

/**
 * Get difference of `result` from the exact `expected`
 */
template <typename Result, typename Derived,
          typename T=diff_scalar_type<Result, Derived>>
diff_stats<T> make_diff(const Result &result,
                        const Eigen::MatrixBase<Derived> &expected,
                        double atol=0)
{
    if ((int)result.size() != expected.size())
        throw size_mismatch();

    typedef make_real_type<T> real;
    typedef typename eigen<T>::matrix matrix;
    const double eps = std::numeric_limits<double>::epsilon();

    diff_stats<T> diff;
    // observations() is the count consistent with var() (for autocorr_result
    // the number of bins at the binning level used), not count().
    diff.nmeas = result.observations();
    diff.mean = result.mean() - expected;

    // Scale per component; a component that is identically zero is scaled
    // by the expected value, so that a nonzero one rejects.
    const column<real> var = result.var();
    const column<real> mean_abs = result.mean().cwiseAbs();
    const column<real> exp_abs = expected.cwiseAbs().template cast<real>();
    column<real> scale = (var + mean_abs.cwiseAbs2()).cwiseSqrt();
    for (int i = 0; i != scale.rows(); ++i)
        if (scale[i] == 0)
            scale[i] = exp_abs[i];

    // Rounding of the mean, of the expected value and of the subtraction
    const double rel_error = accumulation_error(result.count());
    const column<real> diff_error = rel_error * scale + 8 * eps * exp_abs;

    matrix cov = get_cov(result);
    apply_atol(cov, diff.mean, var, atol);
    matrix vecs;
    std::vector<int> keep;
    whiten<T>(cov, !traits<Result>::HAVE_COV, diff.mean, diff.var,
              diff_error, scale, rel_error, vecs, keep);
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
 * with r the number of directions of S above rounding (see whiten(), here
 * with the scale sqrt((var_1 + mean_1^2)/N_1 + (var_2 + mean_2^2)/N_2)). The
 * result is returned as a one-sample difference with nu + 1 observations, so
 * that Hotelling's T^2 is referred to F(r, nu + 1 - r).
 */
template <typename Result1, typename Result2,
          typename T=joined_value_type<Result1, Result2>>
diff_stats<T> pool_var(const Result1 &r1, const Result2 &r2, double atol=0)
{
    if (r1.size() != r2.size())
        throw size_mismatch();

    typedef make_real_type<T> real;
    typedef typename eigen<T>::matrix matrix;
    const double obs1 = r1.observations();
    const double obs2 = r2.observations();

    const column<real> var1 = r1.var(), var2 = r2.var();
    const column<real> s1 = (var1 + r1.mean().cwiseAbs2()).cwiseSqrt();
    const column<real> s2 = (var2 + r2.mean().cwiseAbs2()).cwiseSqrt();
    const column<real> scale = (s1.cwiseAbs2() / obs1 + s2.cwiseAbs2() / obs2).cwiseSqrt();

    // Rounding of each mean, and of the entries of C_k/N_k relative to
    // scale_i scale_j (Cauchy-Schwarz bounds the sum by the larger one).
    const double err1 = accumulation_error(r1.count());
    const double err2 = accumulation_error(r2.count());
    const column<real> diff_error = err1 * s1 + err2 * s2;

    matrix c1 = get_cov(r1) / obs1;
    matrix c2 = get_cov(r2) / obs2;
    diff_stats<T> diff;
    diff.mean = r1.mean() - r2.mean();
    {
        // atol in the units of the data: both results exact there
        column<real> var_max = var1.cwiseAbs().cwiseMax(var2.cwiseAbs());
        column<T> mean = diff.mean;
        apply_atol(c1, mean, var_max, atol);
        apply_atol(c2, diff.mean, var_max, atol);
    }

    matrix v;
    std::vector<int> keep;
    const bool diagonal = !traits<Result1>::HAVE_COV && !traits<Result2>::HAVE_COV;
    whiten<T>(matrix(c1 + c2), diagonal, diff.mean, diff.var, diff_error,
              scale, std::max(err1, err2), v, keep);

    // Whitened covariances of the two means along the axes above rounding
    const column<real> sinv = scale.unaryExpr(
                    [](real s) { return s == 0 ? real(1) : 1 / s; });
    double denom = 0;
    for (int k = 0; k != 2; ++k) {
        const matrix &ck = k == 0 ? c1 : c2;
        const matrix wk = v.adjoint() * sinv.asDiagonal() * ck * sinv.asDiagonal() * v;
        const double nk = k == 0 ? obs1 : obs2;
        double tr = 0, tr2 = 0;
        for (int i : keep) {
            tr += std::real(wk(i, i)) / diff.var[i];
            for (int j : keep)
                tr2 += std::norm(wk(i, j)) / (diff.var[i] * diff.var[j]);
        }
        denom += (tr2 + tr * tr) / (nk - 1);
    }
    const double r = keep.size();
    // Both results exact: no variance to estimate, the difference is
    // compared with the rounding level only.
    const double nu = keep.empty() ? obs1 + obs2 - 2 : (r + r * r) / denom;

    diff.nmeas = nu + 1;
    // t2_test expects per-observation variances, stderror^2 * nmeas
    diff.var *= diff.nmeas;
    return diff;
}

}}} /* namespace alps::alea::internal */
