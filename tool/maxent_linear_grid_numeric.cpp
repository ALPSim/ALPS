// Numerical regression test for the linear FREQUENCY_GRID of the maxent
// tool. ContiParameters allocates t_array_ with nfreq_+1 knots; leaving the
// last knot unwritten corrupts the final bin width to ~(OMEGA_MIN-OMEGA_MAX),
// which blows the default-model normalisation up to +-inf and makes every
// linear-grid run return an all-NaN spectrum.
//
// The test synthesises a clean fermionic G(tau) from a known two-Gaussian
// spectral function (peaks at omega = +-1.5), runs MaxEntSimulation
// in-process on the linear grid, and asserts the recovered spectrum is
// finite and satisfies the sum rule int A(omega) domega ~= 1.
#include "maxent.hpp"

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/ngs/params.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {
bool all_finite(const std::vector<double>& v) {
  if (v.empty()) return false;
  for (std::size_t i = 0; i < v.size(); ++i)
    if (!std::isfinite(v[i])) return false;
  return true;
}
double trapz(const std::vector<double>& x, const std::vector<double>& y) {
  double s = 0.0;
  for (std::size_t i = 0; i + 1 < y.size(); ++i)
    s += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i]);
  return s;
}
}

int main() {
  const double beta = 10.0;
  const int n_tau = 40;
  const int nom = 401;
  const double om_lo = -8.0, om_hi = 8.0;

  std::vector<double> tau(n_tau), om(nom), rho(nom);
  for (int i = 0; i < n_tau; ++i) tau[i] = beta * double(i) / double(n_tau - 1);
  for (int j = 0; j < nom; ++j) om[j] = om_lo + (om_hi - om_lo) * double(j) / double(nom - 1);
  const double dom = om[1] - om[0];
  for (int j = 0; j < nom; ++j)
    rho[j] = std::exp(-std::pow(om[j] - 1.5, 2) / 0.4)
           + std::exp(-std::pow(om[j] + 1.5, 2) / 0.4);
  double rho_int = 0.0;
  for (int j = 0; j + 1 < nom; ++j) rho_int += 0.5 * (rho[j] + rho[j + 1]) * dom;
  for (int j = 0; j < nom; ++j) rho[j] /= rho_int;

  std::vector<double> Gin(n_tau);
  for (int i = 0; i < n_tau; ++i) {
    double acc = 0.0;
    for (int j = 0; j < nom; ++j)
      acc += std::exp(-tau[i] * om[j]) / (1.0 + std::exp(-beta * om[j])) * rho[j];
    Gin[i] = -acc * dom;
  }
  std::vector<double> errors(n_tau, 1e-3);

  alps::params p;
  p["BETA"] = beta;
  p["NDAT"] = n_tau;
  p["NFREQ"] = 400;
  p["NORM"] = 1.0;
  p["KERNEL"] = std::string("fermionic");
  p["DATASPACE"] = std::string("time");
  p["OMEGA_MIN"] = om_lo;
  p["OMEGA_MAX"] = om_hi;
  p["N_ALPHA"] = 12;
  p["ALPHA_MIN"] = 0.1;
  p["ALPHA_MAX"] = 100.0;
  p["DEFAULT_MODEL"] = std::string("flat");
  p["FREQUENCY_GRID"] = std::string("linear");
  p["MAX_IT"] = 40;
  p["PARTICLE_HOLE_SYMMETRY"] = 0;
  p["TEXT_OUTPUT"] = 0;
  p["DATA_IN_HDF5"] = 1;
  p["MAX_TIME"] = 600;
  p["VERBOSE"] = 0;

  const char* tmpdir = std::getenv("TMPDIR");
  std::string base = std::string(tmpdir ? tmpdir : "/tmp") + "/maxent_lock_check";
  const std::string in_h5 = base + ".h5";
  const std::string out_h5 = base + ".out.h5";
  p["DATA"] = in_h5;
  {
    alps::hdf5::archive ar(in_h5, "w");
    ar << alps::make_pvp("/Data", Gin);
    ar << alps::make_pvp("/Error", errors);
  }
  std::remove(out_h5.c_str());

  {
    MaxEntSimulation sim(p, out_h5);
    sim.run(boost::function<bool()>([]() { return false; }));
  }

  std::vector<double> w, Aavg, Amax, Achi;
  {
    alps::hdf5::archive ar(out_h5, "r");
    ar >> alps::make_pvp("/spectrum/omega", w);
    ar >> alps::make_pvp("/spectrum/average", Aavg);
    ar >> alps::make_pvp("/spectrum/maximum", Amax);
    ar >> alps::make_pvp("/spectrum/chi", Achi);
  }

  int failures = 0;
  #define REQUIRE(cond, msg) do { if (!(cond)) { std::fprintf(stderr, "FAIL: %s\n", msg); ++failures; } } while (0)

  REQUIRE(!w.empty(), "omega grid empty");
  REQUIRE(all_finite(w), "omega grid not finite");
  REQUIRE(all_finite(Aavg), "A_average not all-finite (NaN regression)");
  REQUIRE(all_finite(Amax), "A_maximum not all-finite (NaN regression)");
  REQUIRE(all_finite(Achi), "A_chi2 not all-finite (NaN regression)");

  if (all_finite(w) && all_finite(Aavg)) {
    const double sumrule = trapz(w, Aavg);
    std::fprintf(stderr, "sumrule = %.4f\n", sumrule);
    REQUIRE(std::abs(sumrule - 1.0) < 0.05, "sum rule violated (expected ~1.0)");
  }

  if (failures) {
    std::fprintf(stderr, "maxent_lock_check: %d check(s) FAILED\n", failures);
    return 1;
  }
  std::fprintf(stderr, "maxent_lock_check: PASS\n");
  return 0;
}
