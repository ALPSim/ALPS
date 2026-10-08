#pragma once
#include "problem.hpp"
#include "solver.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <optional>
#include <vector>

// The contributor's file: the only one that has to be written to add a method.
// Implement groundStateEnergy() against the Problem interface and leave
// makeAlgorithm's name alone; main.cpp knows nothing else about what lives here.
//
// Needs the matrix? Take a const MatrixProblem& in makeAlgorithm, as below; the
// harness then marks problems too large to offer one as N/A for this method.

namespace onboard {

// Exact diagonalization of the single-particle matrix by cyclic Jacobi
// rotations, then fill the lowest numParticles() levels.
class Algorithm final : public Solver {
public:
    explicit Algorithm(const MatrixProblem& p) : p_(p) {}

    std::optional<Estimate> groundStateEnergy() const override {
        const int n = p_.numSites();
        std::vector<double> a = p_.singleParticleMatrix();
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
        std::sort(levels.begin(), levels.end());

        double e = 0.0;
        for (int k = 0; k < p_.numParticles(); ++k) e += levels[k];
        if (!std::isfinite(e)) return std::nullopt;
        return Estimate{e};
    }

    const char* name() const override { return "JacobiED"; }

private:
    const MatrixProblem& p_;
};

// Factory the harness calls.
inline std::unique_ptr<const Solver> makeAlgorithm(const MatrixProblem& p) {
    return std::make_unique<Algorithm>(p);
}

} // namespace onboard
