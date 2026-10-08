#include "algorithm.hpp"
#include "problem.hpp"
#include "solver.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace onboard;

// Usage: algorithm_test [report.md [gate.txt]]
//   report.md  markdown tables of every case and every category
//   gate.txt   one line "<passed> <judged>" for CI's pass-rate gate. CI reads
//              this file rather than stdout, so a contributed algorithm that
//              prints its own "[PASS]" lines cannot move the gate.

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

// One judged case, kept for the summary tables.
struct Row {
    TestCase c;
    std::string problem, solver;
    Outcome outcome;
    std::optional<Estimate> e;
    double us;
};

const char* label(Outcome o) {
    switch (o) {
        case Outcome::Pass: return "PASS";
        case Outcome::Fail: return "FAIL";
        case Outcome::NotApplicable: break;
    }
    return "N/A ";
}

void printRow(const Row& r) {
    std::printf("[%s]%s %-10s %-14s ", label(r.outcome), r.c.required ? " " : "*",
                r.solver.c_str(), r.problem.c_str());
    if (r.outcome == Outcome::NotApplicable) { std::printf("not offered\n"); return; }
    if (r.e) std::printf("E0 = %+.10f +/- %.1e  ref %+.10f", r.e->value, r.e->error, r.c.expected);
    else     std::printf("no answer%34s", "");
    std::printf("  %10.1f us\n", r.us);
}

// Judges one case. Nothing here names a lattice or a method: the label comes
// from Solver::name(), the answer from the hidden catalog. A template, so that
// only the makeAlgorithm overload the contributor actually wrote is
// instantiated.
template <class P>
Row runContributed(const TestCase& c) {
    auto p = view<P>(c.id);
    if (!p) {
        // No such view of this problem, e.g. too large to offer as a matrix.
        // Not counted against the method.
        return {c, makeProblem(c.id)->name(), "-", Outcome::NotApplicable, std::nullopt, 0.0};
    }

    // The clock covers construction as well as the solve: a method may do its
    // real work in the constructor.
    const auto t0 = Clock::now();
    const auto solver = makeAlgorithm(*p);
    const std::optional<Estimate> e = solver->groundStateEnergy();
    const double us = std::chrono::duration<double, std::micro>(Clock::now() - t0).count();

    const bool ok = e && std::abs(e->value - c.expected) < c.tol + 2.0 * e->error;
    return {c, p->name(), solver->name(), ok ? Outcome::Pass : Outcome::Fail, e, us};
}

struct Tally {
    int passed = 0, judged = 0, na = 0;
    double us = 0.0;

    void add(const Row& r) {
        if (r.outcome == Outcome::NotApplicable) { ++na; return; }
        ++judged;
        if (r.outcome == Outcome::Pass) ++passed;
        us += r.us;
    }
    double rate() const { return judged ? 100.0 * passed / judged : 0.0; }
};

void writeReport(const char* path, const std::vector<Row>& rows,
                 const std::map<std::string, Tally>& byCategory, const Tally& overall) {
    std::FILE* f = std::fopen(path, "w");
    if (!f) { std::perror(path); return; }

    std::fprintf(f, "**Overall: %d of %d judged cases passed (%.0f%%), %.1f us total.**\n\n",
                 overall.passed, overall.judged, overall.rate(), overall.us);

    std::fprintf(f, "### By category\n\n");
    std::fprintf(f, "| category | passed | rate | N/A | time (us) |\n");
    std::fprintf(f, "|---|---|---|---|---|\n");
    for (const auto& [name, t] : byCategory)
        std::fprintf(f, "| %s | %d/%d | %.0f%% | %d | %.1f |\n",
                     name.c_str(), t.passed, t.judged, t.rate(), t.na, t.us);
    std::fprintf(f, "| **overall** | **%d/%d** | **%.0f%%** | %d | %.1f |\n\n",
                 overall.passed, overall.judged, overall.rate(), overall.na, overall.us);

    std::fprintf(f, "### By case\n\n");
    std::fprintf(f, "| result | case | category | required | E0 | reference | time (us) |\n");
    std::fprintf(f, "|---|---|---|---|---|---|---|\n");
    for (const Row& r : rows) {
        std::fprintf(f, "| %s | %s | %s | %s | ", label(r.outcome), r.problem.c_str(),
                     r.c.category, r.c.required ? "yes" : "no");
        if (r.e) std::fprintf(f, "%+.10f &plusmn; %.1e", r.e->value, r.e->error);
        else     std::fprintf(f, "%s", r.outcome == Outcome::NotApplicable ? "not offered" : "no answer");
        std::fprintf(f, " | %+.10f | %.1f |\n", r.c.expected, r.us);
    }
    std::fclose(f);
}

} // namespace

int main(int argc, char** argv) {
    std::vector<Row> rows;
    std::map<std::string, Tally> byCategory;
    Tally overall;
    int failures = 0;

    std::printf("--- Contributed algorithm ---\n");
    for (const TestCase& c : testCases()) {
        rows.push_back(runContributed<View>(c));
        const Row& r = rows.back();
        printRow(r);
        byCategory[c.category].add(r);
        overall.add(r);
        if (r.outcome == Outcome::Fail && c.required) ++failures;
    }

    std::printf("\n%-10s %8s %6s %4s %12s\n", "category", "passed", "rate", "N/A", "time (us)");
    for (const auto& [name, t] : byCategory)
        std::printf("%-10s %4d/%-3d %5.0f%% %4d %12.1f\n",
                    name.c_str(), t.passed, t.judged, t.rate(), t.na, t.us);
    std::printf("%-10s %4d/%-3d %5.0f%% %4d %12.1f\n", "overall",
                overall.passed, overall.judged, overall.rate(), overall.na, overall.us);

    std::printf("\n%d required case(s) failed  (* = informational, N/A = not offered)\n",
                failures);

    if (argc > 1) writeReport(argv[1], rows, byCategory, overall);
    if (argc > 2) {
        if (std::FILE* g = std::fopen(argv[2], "w")) {
            std::fprintf(g, "%d %d\n", overall.passed, overall.judged);
            std::fclose(g);
        }
    }
    return failures == 0 ? 0 : 1;
}
