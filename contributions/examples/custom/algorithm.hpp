#pragma once
#include "problem.hpp"
#include "solver.hpp"

#include <cmath>
#include <memory>
#include <optional>
#include <vector>

// A hand-written wrapper, for code that fits none of the contracts: the power
// method finds one level, so it can only answer single-particle problems. It
// says so through provides(), and declines the rest by returning nullopt.

// Defined in power_method.cpp.
double lowest_eigenvalue(const std::vector<double>& a, int n, int max_iter);

namespace onboard {

class Algorithm final : public Solver {
public:
    explicit Algorithm(const MatrixProblem& p) : p_(p) {}

    std::optional<Estimate> groundStateEnergy() const override {
        if (p_.numParticles() != 1) return std::nullopt;
        const double e = lowest_eigenvalue(p_.singleParticleMatrix(), p_.numSites(), 100000);
        if (!std::isfinite(e)) return std::nullopt;
        return Estimate{e};
    }

    const char* name() const override { return "power_method"; }

private:
    const MatrixProblem& p_;
};

inline std::unique_ptr<const Solver> makeAlgorithm(const MatrixProblem& p) {
    return std::make_unique<Algorithm>(p);
}

} // namespace onboard
