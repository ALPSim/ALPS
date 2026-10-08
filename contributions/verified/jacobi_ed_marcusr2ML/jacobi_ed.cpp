// Plain eigenvalue routine with no knowledge of the onboarding harness.
#include <cmath>
#include <vector>

// Eigenvalues of the symmetric n x n row-major matrix `a`, by cyclic Jacobi
// rotations. Returned unsorted.
std::vector<double> jacobi_eigenvalues(std::vector<double> a, int n) {
    auto at = [&](int i, int j) -> double& { return a[i * n + j]; };

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
                for (int k = 0; k < n; ++k) {  // A <- A J
                    const double akp = at(k, p), akq = at(k, q);
                    at(k, p) = c * akp - s * akq;
                    at(k, q) = s * akp + c * akq;
                }
                for (int k = 0; k < n; ++k) {  // A <- J^T A
                    const double apk = at(p, k), aqk = at(q, k);
                    at(p, k) = c * apk - s * aqk;
                    at(q, k) = s * apk + c * aqk;
                }
            }
    }

    std::vector<double> levels(n);
    for (int i = 0; i < n; ++i) levels[i] = at(i, i);
    return levels;
}
