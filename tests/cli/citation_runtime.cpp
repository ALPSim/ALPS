// SPDX-License-Identifier: MIT
#include <alps/utility/cli.hpp>
#include <alps/parapack/parapack.h>
#include <alps/parapack/worker_factory.h>
#include <alps/scheduler.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#ifdef ALPS_HAVE_MPI
#include <boost/io/ios_state.hpp>
#include <boost/mpi/collectives/gather.hpp>
#include <boost/serialization/string.hpp>
#endif

// Exercise the real stdin scheduler with a small deterministic worker, including
// its parallel-worker path (the shipped loop executable registers serial workers).
class citation_worker : public alps::parapack::abstract_worker {
  int steps_ = 0;
public:
  explicit citation_worker(const alps::Parameters&) {}
#ifdef ALPS_HAVE_MPI
  citation_worker(const boost::mpi::communicator&, const alps::Parameters&) {}
#endif
  void init_observables(const alps::Parameters&, alps::ObservableSet& obs) override {
    obs << alps::RealObservable("Steps");
  }
  void run(alps::ObservableSet& obs) override { obs["Steps"] << double(++steps_); }
  bool is_thermalized() const override { return true; }
  double progress() const override { return steps_ / 2.; }
  void load(alps::IDump& dump) override { dump >> steps_; }
  void save(alps::ODump& dump) const override { dump << steps_; }
};

int main(int argc, char** argv) {
  try {
    if (argc < 2) throw std::invalid_argument("Expected a test mode");
    const std::string mode(argv[1]);
    --argc;
    ++argv;
    if (mode == "single") {
      alps::cli_mpi_guard mpi(argc, argv);
      alps::scheduler::SimpleMCFactory<alps::scheduler::DummyMCRun> factory;
      alps::scheduler::start_single(factory, argc, argv, "spinmc");
      alps::scheduler::stop_single(false);
      return 0;
    }
    if (mode == "single-query") {
      alps::scheduler::SimpleMCFactory<alps::scheduler::DummyMCRun> factory;
      if (alps::scheduler::start_single(factory, argc, argv, "spinmc")) return 1;
      // Existing embedded callers clean up unconditionally after --help/etc.
      alps::scheduler::stop_single();
      return 0;
    }
    if (mode == "owned-query") {
      alps::cli_mpi_guard mpi(argc, argv);
      if (!alps::handle_cli_information(argc, argv, "spinmc", [] {})) return 1;
#ifdef ALPS_HAVE_MPI
      int finalized = 0;
      MPI_Finalized(&finalized);
      if (finalized) throw std::runtime_error("Query finalized caller-owned MPI");
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      return 0;
    }
    if (mode == "parapack") {
      alps::parapack::worker_factory::instance()->register_worker<citation_worker>("citation-test");
      alps::parapack::worker_factory::instance()->set_citation_component("looper");
#ifdef ALPS_HAVE_MPI
      alps::parapack::parallel_worker_factory::instance()->register_worker<citation_worker>("citation-test");
      if (argc == 2 && std::string(argv[1]) == "--mpi") {
        alps::cli_mpi_guard mpi(argc, argv);
        boost::mpi::communicator world;
        std::ostringstream output;
        boost::io::basic_ios_rdbuf_saver<char> stdout_saver(std::cout, output.rdbuf());
        const int result = alps::parapack::start(argc, argv);
        stdout_saver.restore();
        // Preserve every rank's output without depending on launcher merge order.
        std::vector<std::string> rank_output;
        boost::mpi::gather(world, output.str(), rank_output, 0);
        if (world.rank() == 0)
          for (const auto& text : rank_output) std::cout << text;
        return result;
      }
#endif
      return alps::parapack::start(argc, argv);
    }
    throw std::invalid_argument("Unknown test mode");
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
