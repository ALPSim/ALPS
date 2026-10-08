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

// Judges one solver on one case and prints the row. Nothing here names a
// lattice or a method: the label comes from Solver::name(), the answer from
// the hidden catalog.
Outcome judge(const Solver& s, const Problem& p, const TestCase& c) {
    const auto t0 = std::chrono::steady_clock::now();
    const std::optional<Estimate> e = s.groundStateEnergy();
    const double ms = std::chrono::duration<double, std::milli>(
                          std::chrono::steady_clock::now() - t0).count();

    const bool ok = e && std::abs(e->value - c.expected) < c.tol + 2.0 * e->error;
    std::printf("[%s]%s %-10s %-13s ", ok ? "PASS" : "FAIL", c.required ? " " : "*",
                s.name(), p.name());
    if (e) std::printf("E0 = %+.10f +/- %.1e  ref %+.10f", e->value, e->error, c.expected);
    else   std::printf("no answer%34s", "");
    std::printf("  %8.1f ms\n", ms);
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
    return judge(*makeAlgorithm(*p), *p, c);
}

} // namespace

int main() {
    int failures = 0;

    std::printf("--- Contributed algorithm ---\n");
    for (const TestCase& c : testCases()) {
        if (runContributed<View>(c) == Outcome::Fail && c.required) ++failures;
    }

    std::printf("\n%d required case(s) failed  (* = informational, N/A = not offered)\n",
                failures);
    return failures == 0 ? 0 : 1;
}
