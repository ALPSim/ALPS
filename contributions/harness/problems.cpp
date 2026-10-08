#include "problem.hpp"

#include <cmath>
#include <map>
#include <string>
#include <utility>

namespace onboard {

const char* quantityName(Quantity q) {
    switch (q) {
        case Quantity::GroundStateEnergy: return "energy";
        case Quantity::Correlation:       return "correlation";
    }
    return "?";
}

namespace {

// Anonymous namespace: these names do not exist outside this translation unit,
// so no algorithm can #include, forward-declare, or special-case them.

using Parameters = std::map<std::string, double>;

class TightBinding final : public MatrixProblem {
public:
    TightBinding(std::string name, std::string lattice, std::string geometry,
                 Parameters params, int n, int particles, std::vector<Hopping> hops,
                 std::vector<double> eps = {})
        : name_(std::move(name)), lattice_(std::move(lattice)), geometry_(std::move(geometry)),
          params_(std::move(params)), n_(n), particles_(particles),
          hops_(std::move(hops)), eps_(std::move(eps)) {
        eps_.resize(n_, 0.0);
    }

    const char* name() const override     { return name_.c_str(); }
    const char* model() const override    { return "tight-binding"; }
    const char* lattice() const override  { return lattice_.c_str(); }
    const char* geometry() const override { return geometry_.c_str(); }

    std::optional<double> parameter(const char* key) const override {
        const auto it = params_.find(key);
        if (it == params_.end()) return std::nullopt;
        return it->second;
    }

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
    std::string name_, lattice_, geometry_;
    Parameters params_;
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

std::vector<Hopping> squareTorus(int L, double t) {
    std::vector<Hopping> b;
    for (int x = 0; x < L; ++x)
        for (int y = 0; y < L; ++y) {
            const int s = x + L * y;
            b.push_back({s, (x + 1) % L + L * y, t});
            b.push_back({s, x + L * ((y + 1) % L), t});
        }
    return b;
}

// The one place that knows the id -> Hamiltonian mapping. Parameter names
// follow ALPS where it has one (L, t); V is a staggered on-site potential +-V.
std::unique_ptr<TightBinding> create(ProblemId id) {
    switch (id) {
        case ProblemId::OpenChain10:
            return std::make_unique<TightBinding>(
                "open chain 10", "open chain lattice", "1d", Parameters{{"L", 10}, {"t", 1}, {"N", 1}},
                10, 1, chain(10, false, 1));
        case ProblemId::HalfFilledRing10:
            return std::make_unique<TightBinding>(
                "ring 10 N=5", "chain lattice", "1d", Parameters{{"L", 10}, {"t", 1}, {"N", 5}},
                10, 5, chain(10, true, 1));
        case ProblemId::BiasedDimer:
            return std::make_unique<TightBinding>(
                "biased dimer", "dimer", "1d", Parameters{{"t", 1}, {"V", 1}, {"N", 1}},
                2, 1, chain(2, false, 1), std::vector<double>{1.0, -1.0});
        case ProblemId::OpenSquare3x3:
            return std::make_unique<TightBinding>(
                "open 3x3", "open square lattice", "2d", Parameters{{"L", 3}, {"t", 1}, {"N", 1}},
                9, 1, openSquare(3, 1));
        case ProblemId::HalfFilledTorus4x4:
            return std::make_unique<TightBinding>(
                "torus 4x4 N=8", "square lattice", "2d", Parameters{{"L", 4}, {"t", 1}, {"N", 8}},
                16, 8, squareTorus(4, 1));
        case ProblemId::OpenSquare4x4:
            return std::make_unique<TightBinding>(
                "open 4x4", "open square lattice", "2d", Parameters{{"L", 4}, {"t", 1}, {"N", 1}},
                16, 1, openSquare(4, 1));
    }
    return nullptr;
}

constexpr int kMaxMatrixSites = 2000;

// --- Reference correlation functions, from the analytic orbitals -----------

constexpr double kPi = 3.14159265358979323846;

// G(i,j) = sum over occupied orbitals phi(i) phi(j).
std::vector<double> fromOrbitals(int n, const std::vector<std::vector<double>>& occupied) {
    std::vector<double> g(static_cast<std::size_t>(n) * n, 0.0);
    for (const auto& phi : occupied)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) g[i * n + j] += phi[i] * phi[j];
    return g;
}

// Lowest orbital of an open chain: sqrt(2/(L+1)) sin(pi (i+1) / (L+1)).
std::vector<double> openChainGround(int L) {
    std::vector<double> phi(L);
    for (int i = 0; i < L; ++i)
        phi[i] = std::sqrt(2.0 / (L + 1)) * std::sin(kPi * (i + 1) / (L + 1));
    return phi;
}

// Lowest orbital of an open L x L square: product of chain orbitals.
std::vector<double> openSquareGround(int L) {
    const std::vector<double> c = openChainGround(L);
    std::vector<double> phi(L * L);
    for (int x = 0; x < L; ++x)
        for (int y = 0; y < L; ++y) phi[x + L * y] = c[x] * c[y];
    return phi;
}

// Ring of L sites with momenta 0, +-1, ..., +-m filled (closed shell):
// G(i,j) = (1/L) [1 + 2 sum_{k=1..m} cos(2 pi k (i-j) / L)].
std::vector<double> closedShellRing(int L, int m) {
    std::vector<double> g(static_cast<std::size_t>(L) * L);
    for (int i = 0; i < L; ++i)
        for (int j = 0; j < L; ++j) {
            double s = 1.0;
            for (int k = 1; k <= m; ++k) s += 2.0 * std::cos(2.0 * kPi * k * (i - j) / L);
            g[i * L + j] = s / L;
        }
    return g;
}

// Dimer with on-site +-V and hopping t = 1, V = 1: the ground state of
// [[1, -1], [-1, -1]] is (1, 1 + sqrt 2), normalized.
std::vector<double> biasedDimerGround() {
    const double b = 1.0 + std::sqrt(2.0), norm = std::sqrt(1.0 + b * b);
    return {1.0 / norm, b / norm};
}

TestCase energyCase(ProblemId id, double e, bool required = true) {
    return {id, Quantity::GroundStateEnergy, e, {}, 1e-8, required};
}

TestCase correlationCase(ProblemId id, std::vector<double> g, bool required = true) {
    return {id, Quantity::Correlation, 0.0, std::move(g), 1e-8, required};
}

} // namespace

std::unique_ptr<const Problem> makeProblem(ProblemId id) {
    return create(id);
}

std::unique_ptr<const MatrixProblem> makeMatrixProblem(ProblemId id) {
    auto p = create(id);
    if (p && p->numSites() > kMaxMatrixSites) return nullptr;
    return p;
}

// Analytic answers (t = 1). The half-filled torus has a degenerate ground
// state, so its correlation function is not unique and it is tested on the
// energy only.
std::vector<TestCase> testCases() {
    return {
        // -2 cos(pi / 11)
        energyCase(ProblemId::OpenChain10,        -1.918985947228995),
        // -(2 + 4 cos(pi/5) + 4 cos(2 pi/5)): k = 0, +-1, +-2 filled
        energyCase(ProblemId::HalfFilledRing10,   -6.472135954999579),
        // -sqrt(V^2 + t^2) with V = 1
        energyCase(ProblemId::BiasedDimer,        -1.414213562373095),
        // -4 cos(pi / 4)
        energyCase(ProblemId::OpenSquare3x3,      -2.828427124746190),
        // levels -2(cos kx + cos ky): -4 once, -2 four times, then zeros
        energyCase(ProblemId::HalfFilledTorus4x4, -12.0),
        // -4 cos(pi / 5)
        energyCase(ProblemId::OpenSquare4x4,      -3.236067977499790, false),

        correlationCase(ProblemId::OpenChain10,      fromOrbitals(10, {openChainGround(10)})),
        correlationCase(ProblemId::HalfFilledRing10, closedShellRing(10, 2)),
        correlationCase(ProblemId::BiasedDimer,      fromOrbitals(2, {biasedDimerGround()})),
        correlationCase(ProblemId::OpenSquare3x3,    fromOrbitals(9, {openSquareGround(3)})),
    };
}

} // namespace onboard
