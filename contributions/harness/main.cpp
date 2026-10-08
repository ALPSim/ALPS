#include "algorithm.hpp"  // generated from the manifest's contract, or supplied (contract: custom)
#include "problem.hpp"
#include "selection.hpp"  // generated from the manifest: which cases were declared
#include "solver.hpp"

#include <algorithm>
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
//   report.md  markdown tables of every case, by geometry and by quantity
//   gate.txt   one line "<passed> <judged>" for CI's pass-rate gate. CI reads
//              this file rather than stdout, so a contributed algorithm that
//              prints its own "[PASS]" lines cannot move the gate.

namespace {

// Which view of a problem the wrapped makeAlgorithm accepts. A plain Problem
// overload also accepts a MatrixProblem, so check Problem first.
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

// One case, kept for the summary tables.
struct Row {
    TestCase c;
    std::string problem, geometry, solver;
    Outcome outcome;
    std::string value, reference, why;  // `why` explains an N/A
    double us = 0.0;
};

const char* label(Outcome o) {
    switch (o) {
        case Outcome::Pass: return "PASS";
        case Outcome::Fail: return "FAIL";
        case Outcome::NotApplicable: break;
    }
    return "N/A ";
}

std::string format(const char* fmt, double a, double b = 0.0) {
    char buf[96];
    std::snprintf(buf, sizeof buf, fmt, a, b);
    return buf;
}

void printRow(const Row& r) {
    std::printf("[%s]%s %-18s %-14s %-11s ", label(r.outcome), r.c.required ? " " : "*",
                r.solver.c_str(), r.problem.c_str(), quantityName(r.c.quantity));
    if (r.outcome == Outcome::NotApplicable) { std::printf("%s\n", r.why.c_str()); return; }
    std::printf("%-30s ref %-16s %10.1f us\n", r.value.c_str(), r.reference.c_str(), r.us);
}

// Runs and judges one case. Nothing here names a lattice or a method: the
// label comes from Solver::name(), the answer from the hidden catalog. A
// template, so that only the makeAlgorithm overload actually generated is
// instantiated.
template <class P>
Row runCase(const TestCase& c) {
    const auto meta = makeProblem(c.id);
    Row r{c, meta->name(), meta->geometry(), "-", Outcome::NotApplicable, "", "", "", 0.0};

    if (const char* why = notSelected(*meta, c.quantity)) { r.why = why; return r; }

    auto p = view<P>(c.id);
    if (!p) { r.why = "no matrix offered (too large)"; return r; }

    // The clock covers construction as well as the computation: a method may
    // do its real work in the constructor.
    const auto t0 = Clock::now();
    const auto solver = makeAlgorithm(*p);
    r.solver = solver->name();
    if (!solver->provides(c.quantity)) { r.why = "quantity not provided"; return r; }

    bool ok = false;
    if (c.quantity == Quantity::GroundStateEnergy) {
        const std::optional<Estimate> e = solver->groundStateEnergy();
        r.us = std::chrono::duration<double, std::micro>(Clock::now() - t0).count();
        ok = e && std::abs(e->value - c.energy) < c.tol + 2.0 * e->error;
        r.value = e ? format("E0 = %+.10f +/- %.1e", e->value, e->error) : "no answer";
        r.reference = format("%+.10f", c.energy);
    } else {
        const std::optional<MatrixEstimate> g = solver->correlation();
        r.us = std::chrono::duration<double, std::micro>(Clock::now() - t0).count();
        const int n = meta->numSites();
        r.reference = format("G (%.0f x %.0f)", n, n);
        if (!g) {
            r.value = "no answer";
        } else if (g->values.size() != c.correlation.size()) {
            r.value = format("wrong size %.0f, want %.0f",
                             static_cast<double>(g->values.size()),
                             static_cast<double>(c.correlation.size()));
        } else {
            double dev = 0.0;
            for (std::size_t k = 0; k < g->values.size(); ++k) {
                const double d = std::abs(g->values[k] - c.correlation[k]);
                if (!std::isfinite(d)) { dev = INFINITY; break; }  // std::max would drop a NaN
                dev = std::max(dev, d);
            }
            ok = dev < c.tol + 2.0 * g->error;
            r.value = format("max |dG| = %.1e", dev);
        }
    }
    r.outcome = ok ? Outcome::Pass : Outcome::Fail;
    return r;
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
    std::string rateText() const { return judged ? format("%.0f%%", rate()) : "-"; }
};

using Groups = std::map<std::string, Tally>;

void printGroups(const char* title, const Groups& groups) {
    std::printf("\n%-12s %8s %6s %4s %12s\n", title, "passed", "rate", "N/A", "time (us)");
    for (const auto& [name, t] : groups)
        std::printf("%-12s %4d/%-3d %6s %4d %12.1f\n",
                    name.c_str(), t.passed, t.judged, t.rateText().c_str(), t.na, t.us);
}

void reportGroups(std::FILE* f, const char* title, const Groups& groups, const Tally& overall) {
    std::fprintf(f, "### By %s\n\n", title);
    std::fprintf(f, "| %s | passed | rate | N/A | time (us) |\n", title);
    std::fprintf(f, "|---|---|---|---|---|\n");
    for (const auto& [name, t] : groups)
        std::fprintf(f, "| %s | %d/%d | %s | %d | %.1f |\n",
                     name.c_str(), t.passed, t.judged, t.rateText().c_str(), t.na, t.us);
    std::fprintf(f, "| **overall** | **%d/%d** | **%s** | %d | %.1f |\n\n",
                 overall.passed, overall.judged, overall.rateText().c_str(), overall.na, overall.us);
}

void writeReport(const char* path, const std::vector<Row>& rows, const Groups& byGeometry,
                 const Groups& byQuantity, const Tally& overall) {
    std::FILE* f = std::fopen(path, "w");
    if (!f) { std::perror(path); return; }

    std::fprintf(f, "**Overall: %d of %d judged cases passed (%.0f%%), %.1f us total.**\n\n",
                 overall.passed, overall.judged, overall.rate(), overall.us);
    reportGroups(f, "geometry", byGeometry, overall);
    reportGroups(f, "quantity", byQuantity, overall);

    std::fprintf(f, "### By case\n\n");
    std::fprintf(f, "| result | case | quantity | geometry | required | value | reference | time (us) |\n");
    std::fprintf(f, "|---|---|---|---|---|---|---|---|\n");
    for (const Row& r : rows) {
        const bool na = r.outcome == Outcome::NotApplicable;
        std::fprintf(f, "| %s | %s | %s | %s | %s | %s | %s | %s |\n", label(r.outcome),
                     r.problem.c_str(), quantityName(r.c.quantity), r.geometry.c_str(),
                     r.c.required ? "yes" : "no", na ? r.why.c_str() : r.value.c_str(),
                     na ? "" : r.reference.c_str(), na ? "" : format("%.1f", r.us).c_str());
    }
    std::fclose(f);
}

} // namespace

int main(int argc, char** argv) {
    std::vector<Row> rows;
    Groups byGeometry, byQuantity;
    Tally overall;
    int failures = 0;

    std::printf("--- Contributed algorithm ---\n");
    for (const TestCase& c : testCases()) {
        rows.push_back(runCase<View>(c));
        const Row& r = rows.back();
        printRow(r);
        byGeometry[r.geometry].add(r);
        byQuantity[quantityName(c.quantity)].add(r);
        overall.add(r);
        if (r.outcome == Outcome::Fail && c.required) ++failures;
    }

    printGroups("geometry", byGeometry);
    printGroups("quantity", byQuantity);
    std::printf("%-12s %4d/%-3d %6s %4d %12.1f\n", "overall",
                overall.passed, overall.judged, overall.rateText().c_str(), overall.na, overall.us);
    std::printf("\n%d required case(s) failed  (* = informational, N/A = not judged)\n",
                failures);

    if (argc > 1) writeReport(argv[1], rows, byGeometry, byQuantity, overall);
    if (argc > 2) {
        if (std::FILE* g = std::fopen(argv[2], "w")) {
            std::fprintf(g, "%d %d\n", overall.passed, overall.judged);
            std::fclose(g);
        }
    }
    return failures == 0 ? 0 : 1;
}
