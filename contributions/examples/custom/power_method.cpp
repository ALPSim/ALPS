// Shifted power iteration for the lowest eigenvalue of a symmetric matrix,
// with no knowledge of the onboarding harness.
#include <cmath>
#include <vector>

// Lowest eigenvalue of the symmetric n x n row-major matrix `a`; returns NaN
// if it has not converged within max_iter steps.
double lowest_eigenvalue(const std::vector<double>& a, int n, int max_iter) {
    double shift = 0.0;  // Gershgorin bound on the largest |eigenvalue|
    for (int i = 0; i < n; ++i) {
        double r = 0.0;
        for (int j = 0; j < n; ++j) r += std::abs(a[i * n + j]);
        shift = std::fmax(shift, r);
    }

    std::vector<double> v(n), w(n);
    for (int i = 0; i < n; ++i) v[i] = 1.0 + 0.1 * i;  // generic start
    double e = 0.0;
    for (int it = 0; it < max_iter; ++it) {
        double norm = 0.0;
        for (double x : v) norm += x * x;
        norm = std::sqrt(norm);
        for (double& x : v) x /= norm;

        e = 0.0;
        double r2 = 0.0;
        for (int i = 0; i < n; ++i) {
            w[i] = 0.0;
            for (int j = 0; j < n; ++j) w[i] += a[i * n + j] * v[j];
            e += v[i] * w[i];
        }
        for (int i = 0; i < n; ++i) r2 += (w[i] - e * v[i]) * (w[i] - e * v[i]);
        if (std::sqrt(r2) < 1e-10) return e;

        for (int i = 0; i < n; ++i) v[i] = shift * v[i] - w[i];
    }
    return NAN;
}
