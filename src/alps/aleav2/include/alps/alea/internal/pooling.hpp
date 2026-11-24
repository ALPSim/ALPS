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

#include <cmath>
#include <limits>
#include <type_traits>
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
 * Default relative accuracy to which exact components (symmetries, sum
 * rules, constants) are assumed to be computed by the measurement, see
 * test_mean().
 */
constexpr double default_rtol = 64 * std::numeric_limits<double>::epsilon();

/**
 * Components of a difference for the test of Srivastava and Du (2008), used
 * when there are too few observations for Hotelling's T^2 (see
 * srivastava_du_test()): the squared differences in units of their standard
 * errors, the degrees of freedom of each standard error (infinite for a
 * variance that is known, not estimated), and the correlation matrix of the
 * components, empty if they are taken as independent (no covariance).
 */
struct component_stats
{
    Eigen::VectorXd z2;
    Eigen::VectorXd nu;
    Eigen::MatrixXd corr;
    bool complex_data;
};

/**
 * Difference of a result from a reference, in coordinates where its
 * components are uncorrelated: `mean` and `var` are the difference and the
 * per-observation variance along each principal axis of the covariance
 * scaled per component (see whiten()), `nmeas` the number of observations
 * that enter the F distribution of Hotelling's T^2 (one sample). T^2 does not
 * depend on the scaling. `components` is the same difference per component,
 * for the test used when there are too few observations for T^2: when
 * `nmeas` is not more than the number of tested directions plus one, or for
 * two samples, when the smaller one has fewer than 1.5 observations per
 * tested direction (`min_obs`, infinite for one sample; see pool_var()).
 */
template <typename T>
struct diff_stats
{
    column<T> mean;
    column<make_real_type<T>> var;
    double nmeas;
    double min_obs;
    component_stats components;
};

/**
 * Bound on the relative rounding error of a mean, variance or covariance
 * accumulated over `count` samples, relative to s = sqrt(var + mean^2) of
 * each component. The sums are accumulated one sample (or batch) after the
 * other, so the error of a mean is bounded by count * eps/2 * mean(|x|) and
 * that of a (co)variance by about 3 count * eps/2 * s^2. The bound is reached
 * in practice: a constant such as 0.1 is rounded the same way at every
 * addition, and its mean and variance are off by about 0.1 count eps. This
 * returns 2 (count + 4) eps; the 4 covers the few roundings of finalize()
 * (division by the count, subtraction of the squared mean), which do not grow
 * with the count.
 */
inline double accumulation_error(double count)
{
    return 2 * (count + 4) * std::numeric_limits<double>::epsilon();
}

/**
 * Covariance and difference in principal axes, scaled per component.
 *
 * Components can be exact: constant, or related by a symmetry or sum rule
 * (duplicates, sums), which makes the covariance singular. How exact they are
 * is set by the measurement, not by the sampling, so the caller gives the
 * tolerances per component, in the units of the data: `tol_var` for the
 * fluctuation (the square root of a diagonal entry of `cov`) and `tol_diff`
 * for the difference. The covariance is first scaled to R = D^-1 C D^-1 with
 * D = diag(scale), so that tolerances and rounding are relative and do not
 * depend on the units of each component: a cut relative to the largest
 * eigenvalue of C would take a component with a small but real variance next
 * to a large one for exact. T^2 does not change, since
 * d^T C^-1 d = (D^-1 d)^T R^-1 (D^-1 d).
 *
 * R is then diagonalized. Its eigenvalues below the cut, the squared
 * tolerance rotated to that axis plus the rounding of an eigenvalue
 * (`var_floor` relative to scale^2, plus that of the eigensolver), are the
 * exact directions; their number is the numerical rank deficiency of R for
 * that tolerance (Golub and Van Loan, Matrix Computations, Sec. 5.4). Along
 * each principal axis k, with eigenvalue lambda_k and difference delta_k:
 *  - exact, and delta_k within the tolerance: the reference agrees there,
 *    the axis is dropped (var and diff set to zero, which t2_test() skips);
 *  - exact, and delta_k beyond the tolerance: tested with the cut as its
 *    variance, the largest variance the data cannot tell from zero. A
 *    violated symmetry thus rejects, and noise below the resolution of the
 *    variance does not;
 *  - otherwise: tested with the variance lambda_k, and listed in `keep`.
 * NaN in the data is kept and makes the test NaN.
 */
template <typename T, typename Real>
void whiten(const typename eigen<T>::matrix &cov, bool diagonal,
            column<T> &diff, column<Real> &var,
            const column<Real> &tol_var, const column<Real> &tol_diff,
            column<Real> scale, double var_floor,
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
    for (int i = 0; i != n; ++i)
        if (scale[i] == 0)
            scale[i] = 1;       // component identically zero
    const column<T> scaled = diff.cwiseQuotient(scale.template cast<T>());
    const typename eigen<T>::matrix r =
        scale.cwiseInverse().asDiagonal() * cov * scale.cwiseInverse().asDiagonal();

    if (diagonal) {
        vecs = eigen<T>::matrix::Identity(n, n);
        var = r.diagonal().real();
    } else {
        Eigen::SelfAdjointEigenSolver<typename eigen<T>::matrix> es(r);
        vecs = es.eigenvectors();
        var = es.eigenvalues();
        var_floor = n * (var_floor + 8 * n * eps);
    }

    diff = vecs.adjoint() * scaled;
    for (int k = 0; k != n; ++k) {
        // tolerances along the axis; that of the difference includes the
        // rounding of the rotation
        double bound_var = 0, bound_diff = 0;
        for (int i = 0; i != n; ++i) {
            const double v = std::abs(vecs(i, k));
            bound_var += v * tol_var[i] / scale[i];
            bound_diff += v * (tol_diff[i] / scale[i] + n * eps * std::abs(scaled[i]));
        }
        const double cut = bound_var * bound_var + var_floor;

        if (!(var[k] <= cut)) {
            keep.push_back(k);  // tested with its variance (also NaN)
        } else if (std::abs(diff[k]) <= bound_diff) {
            var[k] = 0;
            diff[k] = 0;
        } else {
            var[k] = cut;
        }
    }
}

/**
 * Per-component statistics for srivastava_du_test(), from the difference
 * `diff`, the variance of each mean difference `var_mean` with `nu` degrees
 * of freedom, and its covariance `cov` (null if the components are taken as
 * independent). A component is exact if `var_mean` is within `var_cut`; it is
 * then left out if the difference is within `tol_diff`, and otherwise tested
 * with `var_cut` as a known variance, as in whiten().
 */
template <typename T, typename Real>
component_stats exact_or_tested(const column<T> &diff, const column<Real> &var_mean,
                                const column<Real> &nu, const column<Real> &var_cut,
                                const column<Real> &tol_diff,
                                const typename eigen<T>::matrix *cov)
{
    component_stats out;
    out.complex_data = !std::is_same<T, make_real_type<T>>::value;
    std::vector<int> idx;
    std::vector<double> z2, dof;
    for (int i = 0; i != diff.rows(); ++i) {
        if (!(var_mean[i] <= var_cut[i])) {
            idx.push_back(i);
            z2.push_back(std::norm(diff[i]) / var_mean[i]);
            dof.push_back(nu[i]);
        } else if (!(std::abs(diff[i]) <= tol_diff[i])) {
            z2.push_back(var_cut[i] > 0 ? std::norm(diff[i]) / var_cut[i] : INFINITY);
            dof.push_back(INFINITY);
        }
    }
    out.z2 = Eigen::Map<Eigen::VectorXd>(z2.data(), z2.size());
    out.nu = Eigen::Map<Eigen::VectorXd>(dof.data(), dof.size());
    if (cov) {
        // Correlations of the components tested with their variance; the
        // exact ones (after them) are uncorrelated.
        const int p = z2.size();
        out.corr = Eigen::MatrixXd::Identity(p, p);
        for (int a = 0; a != (int)idx.size(); ++a)
            for (int b = 0; b != (int)idx.size(); ++b)
                if (a != b) {
                    const int i = idx[a], j = idx[b];
                    out.corr(a, b) = std::real((*cov)(i, j))
                            / std::sqrt(std::real((*cov)(i, i)) * std::real((*cov)(j, j)));
                }
    }
    return out;
}

/**
 * Get difference of `result` from the exact `expected`.
 *
 * Tolerances (see whiten()): a component fluctuates within atol + rtol s, and
 * its mean agrees within atol + max(rtol, accumulation_error()) s, with
 * s = sqrt(var + mean^2), in the units of the data; the variance is also
 * resolved only to its rounding, accumulation_error() s^2.
 */
template <typename Result, typename Derived,
          typename T=diff_scalar_type<Result, Derived>>
diff_stats<T> make_diff(const Result &result,
                        const Eigen::MatrixBase<Derived> &expected,
                        double rtol=default_rtol, double atol=0)
{
    if ((int)result.size() != expected.size())
        throw size_mismatch();

    typedef make_real_type<T> real;
    typedef typename eigen<T>::matrix matrix;
    const double eps = std::numeric_limits<double>::epsilon();

    diff_stats<T> diff;
    // observations() is the count consistent with var() (for autocorr_result
    // the number of bins at the binning level used), not count().
    const double nmeas = result.observations();
    diff.nmeas = nmeas;
    diff.min_obs = INFINITY;
    diff.mean = result.mean() - expected;

    // Scale per component; a component that is identically zero is scaled
    // by the expected value, so that a nonzero one rejects.
    const column<real> var = result.var();
    const column<real> exp_abs = expected.cwiseAbs().template cast<real>();
    column<real> scale = (var + result.mean().cwiseAbs2()).cwiseSqrt();
    for (int i = 0; i != scale.rows(); ++i)
        if (scale[i] == 0)
            scale[i] = exp_abs[i];

    // Tolerances of the measurement, and the rounding of the accumulation,
    // of the expected value and of the subtraction.
    const double rounding = accumulation_error(result.count());
    const column<real> tol_var = (atol + rtol * scale.array()).matrix();
    const column<real> tol_diff = (atol + std::max(rtol, rounding) * scale.array()
                                   + 8 * eps * exp_abs.array()).matrix();

    const matrix cov = get_cov(result);
    const bool diagonal = !traits<Result>::HAVE_COV;
    const column<real> var_cut =
            ((tol_var.cwiseAbs2() + rounding * scale.cwiseAbs2()) / nmeas).eval();
    diff.components = exact_or_tested<T, real>(
                diff.mean, column<real>(var / nmeas),
                column<real>::Constant(var.rows(), nmeas - 1), var_cut, tol_diff,
                diagonal ? nullptr : &cov);

    matrix vecs;
    std::vector<int> keep;
    whiten<T>(cov, diagonal, diff.mean, diff.var, tol_var, tol_diff, scale,
              rounding, vecs, keep);
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
 * with r the number of directions of S tested with their variance (see
 * whiten(), here with the scale sqrt((var_1 + mean_1^2)/N_1 + (var_2 +
 * mean_2^2)/N_2)). The result is returned as a one-sample difference with
 * nu + 1 observations, so that Hotelling's T^2 is referred to F(r, nu + 1 - r).
 * The approximation needs a well-determined covariance of each sample. If
 * the sample with the larger variance has no more observations than
 * directions, S is dominated by a singular C_k and T^2 rejects far more often
 * than alpha (47% at alpha = 0.05 for 40 components, N = 20 and 30,
 * variances 9 and 1), and just above it still too often (14% for N = 43 and
 * 200; 6% for N = 60). test_mean() therefore uses the per-component test if
 * the smaller sample has fewer than 1.5 observations per direction
 * (`min_obs`).
 * Each result has the tolerances of make_diff().
 */
template <typename Result1, typename Result2,
          typename T=joined_value_type<Result1, Result2>>
diff_stats<T> pool_var(const Result1 &r1, const Result2 &r2,
                       double rtol=default_rtol, double atol=0)
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

    // Tolerances as in make_diff(), combined for the difference of the means
    // and for the entries of S (Cauchy-Schwarz bounds the rounding of the sum
    // by the larger one).
    const double err1 = accumulation_error(r1.count());
    const double err2 = accumulation_error(r2.count());
    const column<real> t1 = (atol + rtol * s1.array()).matrix();
    const column<real> t2 = (atol + rtol * s2.array()).matrix();
    const column<real> tol_var = (t1.cwiseAbs2() / obs1 + t2.cwiseAbs2() / obs2).cwiseSqrt();
    const column<real> tol_diff = (2 * atol + std::max(rtol, err1) * s1.array()
                                   + std::max(rtol, err2) * s2.array()).matrix();
    const double rounding = std::max(err1, err2);

    const matrix c1 = get_cov(r1) / obs1;
    const matrix c2 = get_cov(r2) / obs2;
    const matrix s = c1 + c2;
    diff_stats<T> diff;
    diff.mean = r1.mean() - r2.mean();
    diff.min_obs = std::min(obs1, obs2);

    // Per component: Welch-Satterthwaite degrees of freedom
    const bool diagonal = !traits<Result1>::HAVE_COV && !traits<Result2>::HAVE_COV;
    const column<real> a1 = c1.diagonal().real(), a2 = c2.diagonal().real();
    column<real> nu(a1.rows());
    for (int i = 0; i != nu.rows(); ++i)
        nu[i] = (a1[i] + a2[i]) * (a1[i] + a2[i])
              / (a1[i] * a1[i] / (obs1 - 1) + a2[i] * a2[i] / (obs2 - 1));
    const column<real> var_cut = (tol_var.cwiseAbs2() + rounding * scale.cwiseAbs2()).eval();
    diff.components = exact_or_tested<T, real>(diff.mean, column<real>(a1 + a2), nu,
                                               var_cut, tol_diff, diagonal ? nullptr : &s);

    matrix v;
    std::vector<int> keep;
    whiten<T>(s, diagonal, diff.mean, diff.var, tol_var, tol_diff, scale,
              rounding, v, keep);

    // Whitened covariances of the two means along the axes tested with their
    // variance
    const column<real> sinv = scale.unaryExpr(
                    [](real x) { return x == 0 ? real(1) : 1 / x; });
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
    // compared with the tolerance only.
    const double nu_ky = keep.empty() ? obs1 + obs2 - 2 : (r + r * r) / denom;

    diff.nmeas = nu_ky + 1;
    // t2_test expects per-observation variances, stderror^2 * nmeas
    diff.var *= diff.nmeas;
    return diff;
}

}}} /* namespace alps::alea::internal */
