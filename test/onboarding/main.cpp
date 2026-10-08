#include "algorithm.hpp"
#include "problem.hpp"
#include "solver.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <optional>
#include <type_traits>
#include <utility>

using namespace onboard;

namespace {

// Which view of a problem the contributed makeAlgorithm accepts. A plain
// Problem overload also accepts a MatrixProblem, so check Problem first.
template <class P, class = void>
struct Accepts : std::false_type {};
template <class P>
struct Accepts<P, std::void_t<decltype(makeAlgorithm(std::declval<const P&>()))>>
    : std::true_type {};

static_assert(Accepts<Problem>::value || Accepts<MatrixProblem>::value,
              "makeAlgorithm must take const Problem& or const MatrixProblem&");
using View = std::conditional_t<Accepts<Problem>::value, Problem, MatrixProblem>;

template <class P> std::unique_ptr<const P> view(ProblemId id);
template <> [[maybe_unused]] std::unique_ptr<const Problem> view<Problem>(ProblemId id) {
    return makeProblem(id);
}
template <> [[maybe_unused]] std::unique_ptr<const MatrixProblem> view<MatrixProblem>(ProblemId id) {
    return makeMatrixProblem(id);
}

enum class Outcome { Pass, Fail, NotApplicable };

using Clock = std::chrono::steady_clock;

// Wall time spent inside the contributed code, summed over all cases.
double totalUs = 0.0;

// Judges one answer and prints the row. Nothing here names a lattice or a
// method: the label comes from Solver::name(), the answer from the hidden
// catalog.
Outcome judge(const char* solver, const Problem& p, const TestCase& c,
              const std::optional<Estimate>& e, double us) {
    const bool ok = e && std::abs(e->value - c.expected) < c.tol + 2.0 * e->error;
    std::printf("[%s]%s %-10s %-13s ", ok ? "PASS" : "FAIL", c.required ? " " : "*",
                solver, p.name());
    if (e) std::printf("E0 = %+.10f +/- %.1e  ref %+.10f", e->value, e->error, c.expected);
    else   std::printf("no answer%34s", "");
    std::printf("  %10.1f us\n", us);
    return ok ? Outcome::Pass : Outcome::Fail;
}

// Hands the contributed method the view it asked for. A template, so that only
// the makeAlgorithm overload the contributor actually wrote gets instantiated.
template <class P>
Outcome runContributed(const TestCase& c) {
    auto p = view<P>(c.id);
    if (!p) {
        // No such view of this problem, e.g. too large to offer as a matrix.
        // Not counted against the method.
        std::printf("[N/A ]%s %-10s %s\n", c.required ? " " : "*", "-",
                    makeProblem(c.id)->name());
        return Outcome::NotApplicable;
    }

    // The clock covers construction as well as the solve: a method may do its
    // real work in the constructor.
    const auto t0 = Clock::now();
    const auto solver = makeAlgorithm(*p);
    const std::optional<Estimate> e = solver->groundStateEnergy();
    const double us = std::chrono::duration<double, std::micro>(Clock::now() - t0).count();
    totalUs += us;

    return judge(solver->name(), *p, c, e, us);
}

} // namespace

int main() {
    int failures = 0;

    std::printf("--- Contributed algorithm ---\n");
    for (const TestCase& c : testCases()) {
        if (runContributed<View>(c) == Outcome::Fail && c.required) ++failures;
    }

    std::printf("\nTotal time in contributed code: %.1f us\n", totalUs);
    std::printf("%d required case(s) failed  (* = informational, N/A = not offered)\n",
                failures);
    return failures == 0 ? 0 : 1;
}
