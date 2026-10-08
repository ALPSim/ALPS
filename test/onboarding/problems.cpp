#include "problem.hpp"

#include <string>
#include <utility>

namespace onboard {
namespace {

// Anonymous namespace: these names do not exist outside this translation unit,
// so no algorithm can #include, forward-declare, or special-case them.

class TightBinding final : public MatrixProblem {
public:
    TightBinding(std::string name, int n, int particles, std::vector<Hopping> hops,
                 std::vector<double> eps = {})
        : name_(std::move(name)), n_(n), particles_(particles),
          hops_(std::move(hops)), eps_(std::move(eps)) {
        eps_.resize(n_, 0.0);
    }

    const char* name() const override                    { return name_.c_str(); }
    int numSites() const override                        { return n_; }
    int numParticles() const override                    { return particles_; }
    const std::vector<Hopping>& hoppings() const override { return hops_; }
    double onsite(int site) const override               { return eps_[site]; }

    std::vector<double> singleParticleMatrix() const override {
        std::vector<double> h(static_cast<std::size_t>(n_) * n_, 0.0);
        for (int i = 0; i < n_; ++i) h[i * n_ + i] = eps_[i];
        for (const Hopping& b : hops_) {
            h[b.i * n_ + b.j] -= b.t;
            h[b.j * n_ + b.i] -= b.t;
        }
        return h;
    }

private:
    std::string name_;
    int n_, particles_;
    std::vector<Hopping> hops_;
    std::vector<double> eps_;
};

std::vector<Hopping> chain(int L, bool periodic, double t) {
    std::vector<Hopping> b;
    for (int i = 0; i + 1 < L; ++i) b.push_back({i, i + 1, t});
    if (periodic) b.push_back({L - 1, 0, t});
    return b;
}

std::vector<Hopping> openSquare(int L, double t) {
    std::vector<Hopping> b;
    for (int x = 0; x < L; ++x)
        for (int y = 0; y < L; ++y) {
            const int s = x + L * y;
            if (x + 1 < L) b.push_back({s, s + 1, t});
            if (y + 1 < L) b.push_back({s, s + L, t});
        }
    return b;
}

// The one place that knows the id -> Hamiltonian mapping.
std::unique_ptr<TightBinding> create(ProblemId id) {
    switch (id) {
        case ProblemId::OpenChain10:
            return std::make_unique<TightBinding>("open chain 10", 10, 1, chain(10, false, 1));
        case ProblemId::HalfFilledRing10:
            return std::make_unique<TightBinding>("ring 10 N=5", 10, 5, chain(10, true, 1));
        case ProblemId::BiasedDimer:
            return std::make_unique<TightBinding>("biased dimer", 2, 1, chain(2, false, 1),
                                                  std::vector<double>{1.0, -1.0});
        case ProblemId::OpenSquare4x4:
            return std::make_unique<TightBinding>("open 4x4", 16, 1, openSquare(4, 1));
    }
    return nullptr;
}

constexpr int kMaxMatrixSites = 2000;

} // namespace

std::unique_ptr<const Problem> makeProblem(ProblemId id) {
    return create(id);
}

std::unique_ptr<const MatrixProblem> makeMatrixProblem(ProblemId id) {
    auto p = create(id);
    if (p && p->numSites() > kMaxMatrixSites) return nullptr;
    return p;
}

// Hard-coded analytic ground-state energies (t = 1).
std::vector<TestCase> testCases() {
    return {
        // -2 cos(pi / 11)
        {ProblemId::OpenChain10,      -1.918985947228995, 1e-8, true},
        // -(2 + 4 cos(pi/5) + 4 cos(2 pi/5)): k = 0, +-1, +-2 filled
        {ProblemId::HalfFilledRing10, -6.472135954999579, 1e-8, true},
        // -sqrt(eps^2 + t^2) with eps = 1
        {ProblemId::BiasedDimer,      -1.414213562373095, 1e-8, true},
        // -4 cos(pi / 5)
        {ProblemId::OpenSquare4x4,    -3.236067977499790, 1e-8, false},
    };
}

} // namespace onboard
