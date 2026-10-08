// Plain eigensystem routine with no knowledge of the onboarding harness.
#include <cmath>
#include <vector>

// Eigenvalues and normalized eigenvectors of the symmetric n x n row-major
// matrix `a`, by cyclic Jacobi rotations. Eigenvector k is written to
// vectors[k*n .. k*n + n). Returned unsorted.
void jacobi_eigensystem(std::vector<double> a, int n,
                        std::vector<double>& values, std::vector<double>& vectors) {
    auto at = [&](int i, int j) -> double& { return a[i * n + j]; };

    // V accumulates the rotations; its columns become the eigenvectors.
    std::vector<double> V(static_cast<std::size_t>(n) * n, 0.0);
    for (int i = 0; i < n; ++i) V[i * n + i] = 1.0;

    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0.0;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) off += at(i, j) * at(i, j);
        if (off < 1e-30) break;

        for (int p = 0; p < n; ++p)
            for (int q = p + 1; q < n; ++q) {
                if (std::abs(at(p, q)) < 1e-300) continue;
                const double theta = (at(q, q) - at(p, p)) / (2.0 * at(p, q));
                const double t = std::copysign(1.0, theta) /
                                 (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < n; ++k) {  // A <- A J,  V <- V J
                    const double akp = at(k, p), akq = at(k, q);
                    at(k, p) = c * akp - s * akq;
                    at(k, q) = s * akp + c * akq;
                    const double vkp = V[k * n + p], vkq = V[k * n + q];
                    V[k * n + p] = c * vkp - s * vkq;
                    V[k * n + q] = s * vkp + c * vkq;
                }
                for (int k = 0; k < n; ++k) {  // A <- J^T A
                    const double apk = at(p, k), aqk = at(q, k);
                    at(p, k) = c * apk - s * aqk;
                    at(q, k) = s * apk + c * aqk;
                }
            }
    }

    values.resize(n);
    vectors.resize(static_cast<std::size_t>(n) * n);
    for (int k = 0; k < n; ++k) {
        values[k] = at(k, k);
        for (int i = 0; i < n; ++i) vectors[k * n + i] = V[i * n + k];
    }
}
