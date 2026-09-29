// SPDX-License-Identifier: MIT
// The same executable acts as a small external solver, so the test exercises
// real input/output handoffs and deletion of the temporary solver result.
#include "externalsolver.h"
#include <alps/hdf5.hpp>
#include <alps/utility/citation_provenance.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <stdexcept>

namespace {
struct temporary_directory {
  boost::filesystem::path path = boost::filesystem::temp_directory_path() /
                                boost::filesystem::unique_path("alps-dmft-citations-%%%%-%%%%");
  temporary_directory() { boost::filesystem::create_directory(path); }
  ~temporary_directory() { boost::system::error_code error; boost::filesystem::remove_all(path, error); }
};
void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}
}

int main(int argc, char** argv) {
  try {
    if (argc == 3) {
      alps::hdf5::archive input(argv[1], "r");
      alps::Parameters p;
      input["/parameters"] >> p;
      auto history = alps::read_citations(input);
      require(history.size() == 1 && history[0].component == "dmft", "Driver citations missing from solver input");
      alps::hdf5::archive output(argv[2], "w");
      alps::write_citations(output, history);
      alps::write_citations(output, "interaction");
      const int n = p["N"];
      itime_green_function_t tau(n + 1, 1, 2);
      for (int i = 0; i <= n; ++i) for (int f = 0; f < 2; ++f) tau(i, f) = -0.5;
      tau.write_hdf5(output, "/G_tau");
      if (p.defined("NMATSUBARA")) {
        const int nw = p["NMATSUBARA"];
        matsubara_green_function_t omega(nw, 1, 2);
        for (int i = 0; i < nw; ++i) for (int f = 0; f < 2; ++f) omega(i, f) = {0., -1.};
        omega.write_hdf5(output, "/G_omega");
      }
      return 0;
    }
    temporary_directory directory;
    alps::Parameters p;
    p["TMPNAME"] = (directory.path / "solver").string();
    p["N"] = 2; p["SITES"] = 1; p["FLAVORS"] = 2;
    ExternalSolver solver(boost::filesystem::absolute(argv[0]));
    itime_green_function_t tau(3, 1, 2);
    for (int i = 0; i < 3; ++i) for (int f = 0; f < 2; ++f) tau(i, f) = -0.5;
    require(solver.solve(tau, p)(0, 0) == -0.5, "Imaginary-time result changed");
    auto history = solver.citations();
    require(history.size() == 2, "Temporary solver output citations were lost");
    require(!boost::filesystem::exists(directory.path / "solver.out.h5"), "Temporary result was not removed");
    // A new solver exercises the frequency path independently.
    ExternalSolver frequency_solver(boost::filesystem::absolute(argv[0]));
    p["NMATSUBARA"] = 2;
    matsubara_green_function_t omega(2, 1, 2);
    for (int i = 0; i < 2; ++i) for (int f = 0; f < 2; ++f) omega(i, f) = {0., -1.};
    auto result = frequency_solver.solve_omega(omega, p);
    require(result.first(0, 0) == std::complex<double>(0., -1.), "Frequency result changed");
    alps::merge_citations(history, frequency_solver.citations());
    require(history.size() == 2, "Repeated solver use duplicated recommendations");
    alps::hdf5::archive final((directory.path / "dmft.h5").string(), "w");
    alps::write_citations(final, history);
    alps::write_citations(final, "dmft");
    require(alps::read_citations(final).size() == 2, "Final DMFT result lost solver recommendations");
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
