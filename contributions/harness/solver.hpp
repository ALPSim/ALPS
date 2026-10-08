#pragma once
#include <optional>

namespace onboard {

// The abstract base class a contributor plugs into. Solvers depend only on the
// Problem interface. nullopt means "no trustworthy answer". `error` is the
// method's own error bar: 0 for an exact method, a statistical error for a
// stochastic one; the judge widens the tolerance by it.

struct Estimate {
    double value;
    double error = 0.0;
};

class Solver {
public:
    virtual ~Solver() = default;
    virtual std::optional<Estimate> groundStateEnergy() const = 0;
    virtual const char* name() const { return "Algorithm"; }
};

} // namespace onboard
