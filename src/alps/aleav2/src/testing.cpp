/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */

#include <alps/alea/testing.hpp>

namespace alps { namespace alea {

template <typename T>
t2_result t2_test(const column<T> &diff,
                  const column<typename make_real<T>::type> &var,
                  double nmeas, size_t pools, double atol)
{
    if (diff.rows() != var.rows())
        throw std::invalid_argument("Size mismatch between diff and var");
    if (pools != 1 && pools != 2)
        throw std::invalid_argument("Pools must be 1 or 2");
    if (nmeas <= 0)
        throw std::invalid_argument("Must have at least 1 measurement");

    if ((var.array() < -atol).any())
        throw std::invalid_argument("Variances must be positive");

    // remove data points where both the variance and difference are zero
    // within atol while propagating intrinsic NaNs. Exact components are
    // handled before, where their scale is known (internal::whiten()).
    double t2 = 0;
    size_t nsize = 0;
    for (int i = 0; i != diff.rows(); ++i) {
        if (std::abs(var[i]) <= atol && std::norm(diff[i]) <= atol)
            continue;

        t2 += std::norm(diff[i]) / var[i];
        ++nsize;
    }
    t2 *= nmeas;

    // nothing left to test: the results agree
    if (nsize == 0)
        return t2_result(0, 0, nmeas);

    // compute degrees of freedom.
    double dof = nmeas - nsize;
    if (dof <= pools)
        return t2_result(NAN, nsize, dof);

    // normalize T2 score to match F distribution
    double score = dof/(nsize * (nmeas - pools)) * t2;

    return t2_result(score, nsize, dof);
}

namespace internal {

t2_result srivastava_du_test(const component_stats &c)
{
    const int p = c.z2.rows();
    if (p == 0)
        return t2_result(0, 0, 1);              // nothing to test
    if (c.complex_data)
        return t2_result::normal(NAN, p);

    double t = 0, mean0 = 0, var_indep = 0, inv_nu = 0;
    for (int i = 0; i != p; ++i) {
        const double nu = c.nu[i];
        if (!(nu >= srivastava_du_min_dof))
            return t2_result::normal(NAN, p);
        t += c.z2[i];
        if (std::isinf(nu)) {                   // known variance: chi^2_1
            mean0 += 1;
            var_indep += 2;
        } else {                                // t^2 with nu dof: F(1, nu)
            mean0 += nu / (nu - 2);
            var_indep += 2 * nu * nu * (nu - 1) / ((nu - 2) * (nu - 2) * (nu - 4));
            inv_nu += 1 / nu;
        }
    }
    double var0 = var_indep;
    if (c.corr.size() != 0) {
        const double tr2 = c.corr.squaredNorm();
        // p^2 / nu with nu the harmonic mean of the nu_i
        var0 = 2 * (tr2 - p * inv_nu) * (1 + tr2 / std::pow(p, 1.5));
    }
    if (!(var0 > 0))
        return t2_result::normal(NAN, p);
    return t2_result::normal((t - mean0) / std::sqrt(var0), p);
}

} // namespace internal

template t2_result t2_test(const column<double> &diff, const column<double> &var,
                           double nmeas, size_t pools, double atol);
template t2_result t2_test(const column<std::complex<double>> &diff,
                           const column<double> &var, double nmeas, size_t pools,
                           double atol);

}} /* namespace alps::alea */

