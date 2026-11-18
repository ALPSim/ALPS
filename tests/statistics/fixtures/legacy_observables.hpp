// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// What legacy_observables.h5 holds, and how a set of observables is reduced to
// named numbers for comparison with legacy_observables.expected.
#ifndef ALPS_TESTS_STATISTICS_FIXTURES_LEGACY_OBSERVABLES_HPP
#define ALPS_TESTS_STATISTICS_FIXTURES_LEGACY_OBSERVABLES_HPP

#include "../streams.hpp"

#include <alps/alea.h>

#include <map>
#include <string>
#include <valarray>

namespace statistics_test {

// One run: a binned scalar, an unbinned scalar, a two-component vector, a sign
// with a signed observable, and an integer histogram.
inline alps::ObservableSet legacy_run(std::uint64_t seed, std::size_t n, bool with_sign = true) {
    alps::ObservableSet set;
    set << alps::RealObservable("x") << alps::SimpleRealObservable("s") << alps::RealVectorObservable("v")
        << alps::IntHistogramObservable("h", 0, 10);
    if (with_sign)
        set << alps::RealObservable("Sign") << alps::SignedObservable<alps::RealObservable>("O");
    ar1 correlated(seed, 0.7);
    normal xi(seed + 1);
    for (std::size_t t = 0; t < n; ++t) {
        const double x = correlated();
        const double sign = xi.uniform() < 0.8 ? 1. : -1.;
        set["x"] << x;
        set["s"] << x;
        set["v"] << std::valarray<double>({x, 2 * x + xi()});
        const double o = sign * (1 + 0.3 * sign + xi());
        if (with_sign) {
            set["Sign"] << sign;
            set["O"] << o;
        }
        set["h"] << static_cast<int>(10 * xi.uniform());
    }
    return set;
}

typedef std::map<std::string, double> numbers;

// The observables of a single run, as loaded from "/run".
inline numbers describe_run(const alps::ObservableSet & set) {
    numbers out;
    const alps::RealObservable & x = set.get<alps::RealObservable>("x");
    out["x/count"] = x.count();
    out["x/mean"] = x.mean();
    out["x/error"] = x.error();
    out["x/error_level0"] = x.error(0);
    out["x/error_level3"] = x.error(3);
    out["x/variance"] = x.variance();
    out["x/tau"] = x.tau();
    out["x/converged"] = x.converged_errors();
    out["x/bins"] = x.bin_number();
    out["x/binsize"] = x.bin_size();
    out["x/bin_first"] = x.bin_value(0);
    out["x/bin_last"] = x.bin_value(x.bin_number() - 1);
    const alps::SimpleRealObservable & s = set.get<alps::SimpleRealObservable>("s");
    out["s/count"] = s.count();
    out["s/mean"] = s.mean();
    out["s/error"] = s.error();
    const alps::RealVectorObservable & v = set.get<alps::RealVectorObservable>("v");
    out["v/count"] = v.count();
    for (int i = 0; i < 2; ++i) {
        out["v/mean/" + std::to_string(i)] = v.mean()[i];
        out["v/error/" + std::to_string(i)] = v.error()[i];
        out["v/tau/" + std::to_string(i)] = v.tau()[i];
    }
    out["Sign/mean"] = set.get<alps::RealObservable>("Sign").mean();
    const alps::RealObsevaluator o = set["O"];
    out["O/count"] = o.count();
    out["O/mean"] = o.mean();
    out["O/error"] = o.error();
    const alps::IntHistogramObsevaluator h = set["h"];
    out["h/count"] = h.count();
    for (int bin = 0; bin < 10; ++bin)
        out["h/" + std::to_string(bin)] = h[bin];
    return out;
}

// The result of merging two runs, as loaded from "/merged": what a task's
// result group holds.
inline numbers describe_merged(const alps::ObservableSet & set) {
    numbers out;
    const alps::RealObsevaluator x = set["x"];
    out["x/count"] = x.count();
    out["x/mean"] = x.mean();
    out["x/error"] = x.error();
    out["x/tau"] = x.tau();
    out["x/bins"] = x.bin_number();
    const alps::RealVectorObsevaluator v = set["v"];
    for (int i = 0; i < 2; ++i) {
        out["v/mean/" + std::to_string(i)] = v.mean()[i];
        out["v/error/" + std::to_string(i)] = v.error()[i];
    }
    const alps::IntHistogramObsevaluator h = set["h"];
    out["h/count"] = h.count();
    out["h/0"] = h[0];
    out["h/9"] = h[9];
    return out;
}

}

#endif
