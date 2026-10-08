#pragma once
#include <memory>
#include <vector>

namespace onboard {

// The only view a contributor gets of a test problem: spinless fermions on a
// graph,
//
//   H = - sum_bonds t (c+(i) c(j) + c+(j) c(i)) + sum_sites eps(i) n(i),
//
// filled with numParticles() particles. No concrete class is named here, so
// nothing downstream can depend on which lattice or parameters a case uses.

struct Hopping {
    int i, j;
    double t;
};

class Problem {
public:
    virtual ~Problem() = default;
    virtual const char* name() const = 0;
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

// One judged case: the problem, the ground-state energy it has to land on, and
// the category it is reported under. Only `required` cases decide the exit
// code; the rest are reported for information, because no single method
// handles every Hamiltonian.
struct TestCase {
    ProblemId id;
    const char* category;
    double expected;
    double tol;
    bool required;
};

std::vector<TestCase> testCases();

} // namespace onboard
