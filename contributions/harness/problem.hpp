#pragma once
#include <memory>
#include <optional>
#include <vector>

namespace onboard {

// The only view a contributor gets of a test problem: spinless fermions on a
// graph,
//
//   H = - sum_bonds t (c+(i) c(j) + c+(j) c(i)) + sum_sites eps(i) n(i),
//
// filled with numParticles() particles. No concrete class is named here, so
// nothing downstream can depend on which case is being run.

struct Hopping {
    int i, j;
    double t;
};

class Problem {
public:
    virtual ~Problem() = default;
    virtual const char* name() const = 0;

    // What kind of problem this is, in ALPS vocabulary where one exists.
    virtual const char* model() const = 0;     // e.g. "tight-binding"
    virtual const char* lattice() const = 0;   // e.g. "open chain lattice"
    virtual const char* geometry() const = 0;  // "1d" or "2d"

    // A named model parameter (L, t, N, V, ...), or nullopt if this problem
    // does not define it.
    virtual std::optional<double> parameter(const char* key) const = 0;

    virtual int numSites() const = 0;
    virtual int numParticles() const = 0;
    virtual const std::vector<Hopping>& hoppings() const = 0;
    virtual double onsite(int site) const = 0;
};

// Split out so a method that needs the matrix is type-checked at compile time.
class MatrixProblem : public Problem {
public:
    // The numSites() x numSites() single-particle Hamiltonian, row-major.
    virtual std::vector<double> singleParticleMatrix() const = 0;
};

// What a test case checks.
enum class Quantity {
    GroundStateEnergy,
    Correlation,  // G(i,j) = <c+(i) c(j)> in the ground state, row-major
};

const char* quantityName(Quantity q);

enum class ProblemId {
    OpenChain10,
    HalfFilledRing10,
    BiasedDimer,
    OpenSquare3x3,
    HalfFilledTorus4x4,
    OpenSquare4x4,
};

// nullptr for an unknown id.
std::unique_ptr<const Problem> makeProblem(ProblemId id);

// nullptr if the problem is too large to offer as a matrix.
std::unique_ptr<const MatrixProblem> makeMatrixProblem(ProblemId id);

// One judged case: the problem, the quantity checked, and the answer it has to
// land on (`energy` or `correlation`, depending on `quantity`). Only `required`
// cases decide the exit code; the rest are reported for information.
struct TestCase {
    ProblemId id;
    Quantity quantity;
    double energy;
    std::vector<double> correlation;
    double tol;
    bool required;
};

std::vector<TestCase> testCases();

} // namespace onboard
