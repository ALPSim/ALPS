#pragma once
#include "problem.hpp"

#include <optional>
#include <vector>

namespace onboard {

// The abstract base class every algorithm is wrapped in. Solvers depend only on
// the Problem interface. nullopt means "no trustworthy answer". `error` is the
// method's own error bar: 0 for an exact method, a statistical error for a
// stochastic one; the judge widens the tolerance by it.

struct Estimate {
    double value;
    double error = 0.0;
};

struct MatrixEstimate {
    std::vector<double> values;  // row-major
    double error = 0.0;          // bound on any single element
};

class Solver {
public:
    virtual ~Solver() = default;

    // Which quantities this solver computes. Cases for anything else are
    // reported as N/A rather than failed.
    virtual bool provides(Quantity q) const { return q == Quantity::GroundStateEnergy; }

    virtual std::optional<Estimate> groundStateEnergy() const { return std::nullopt; }
    virtual std::optional<MatrixEstimate> correlation() const { return std::nullopt; }

    virtual const char* name() const { return "Algorithm"; }
};

} // namespace onboard
