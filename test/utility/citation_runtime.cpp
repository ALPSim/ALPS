// SPDX-License-Identifier: MIT
#include <alps/utility/cli.hpp>
#include <alps/parapack/parapack.h>
#include <alps/parapack/worker_factory.h>
#include <alps/parapack/clone.h>
#include <alps/parapack/job.h>
#include <alps/parapack/simulation_p.h>
#include <alps/utility/citation_provenance.hpp>
#include <alps/hdf5.hpp>
#include <boost/filesystem.hpp>
#include <alps/scheduler.h>
#include <iostream>
#include <stdexcept>

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
    if (mode == "parapack" || mode == "clone-hdf5" || mode == "aggregate-hdf5") {
      alps::parapack::worker_factory::instance()->register_worker<citation_worker>("citation-test");
      alps::parapack::worker_factory::instance()->set_citation_component("looper");
#ifdef ALPS_HAVE_MPI
      alps::parapack::parallel_worker_factory::instance()->register_worker<citation_worker>("citation-test");
#endif
      if (mode == "aggregate-hdf5") {
        if (argc != 2) throw std::invalid_argument("Expected an aggregate directory");
        const boost::filesystem::path directory(argv[1]);
        alps::parapack::option opt(1, argv, true);
        alps::parapack::evaluator_factory::instance()->register_evaluator<alps::parapack::simple_evaluator>("citation-test");
        alps::Parameters p;
        p["SEED"] = 17; p["ALGORITHM"] = "citation-test"; p["NUM_CLONES"] = 1;
        const auto output = directory / "aggregate.out.h5";
        // An unrelated previous result must not become an input contribution.
        {
          alps::hdf5::archive ar(output.string(), "w");
          alps::write_citations(ar, "dmrg");
        }
        for (int pass = 0; pass < 2; ++pass) {
          alps::clone clone(directory, opt, 0, 0, p, "aggregate", true);
          while (!clone.halted()) clone.run();
          if (pass == 0) {
            alps::hdf5::archive ar((directory / clone.info().dumpfile_h5()).string(), "a");
            alps::write_citations(ar, "interaction"); // a supplied input contribution
          }
          alps::simulation_xml_writer(directory / "aggregate.out.xml", false, true,
                                     p, {}, {clone.info()});
          alps::task task(directory / "aggregate.out.xml");
          task.evaluate(opt);
          auto records = alps::read_citations(output);
          if (records.size() != (pass == 0 ? 3u : 2u))
            throw std::runtime_error("Aggregate did not follow the current clone ancestry");
          for (const auto& record : records)
            if (record.component == "dmrg") throw std::runtime_error("Aggregate inherited unrelated output citations");
          alps::hdf5::archive ar(output.string(), "r");
          double mean = 0.; ar["/simulation/results/Steps/mean/value"] >> mean;
          if (mean != 1.5) throw std::runtime_error("Aggregate numerical result changed");
        }
        return 0;
      }
      if (mode == "clone-hdf5") {
        if (argc != 2) throw std::invalid_argument("Expected a checkpoint directory");
        alps::cli_mpi_guard mpi(argc, argv);
        const boost::filesystem::path directory(argv[1]);
        alps::Parameters p;
        p["SEED"] = 17; p["ALGORITHM"] = "citation-test";
        alps::parapack::option opt(1, argv);
        int rank = 0, size = 1;
#ifdef ALPS_HAVE_MPI
        boost::mpi::communicator world;
        rank = world.rank(); size = world.size();
        alps::clone_create_msg_t message(0, 0, 0, p, "checkpoint", true);
        alps::clone_mpi clone(world, world, directory, opt, message);
#else
        alps::clone clone(directory, opt, 0, 0, p, "checkpoint", true);
#endif
        const auto file = directory / ("clone." + std::to_string(rank) + ".h5");
        {
          alps::hdf5::archive ar(file.string(), "w");
          clone.save(ar);
          alps::write_citations(ar, "looper", "analysis");
          clone.load(ar);
        }
        // Transfer the loaded history to a different checkpoint path.
        {
          alps::hdf5::archive ar((directory / ("copy." + std::to_string(rank) + ".h5")).string(), "w");
          clone.save(ar);
          if (alps::read_citations(ar).size() != 2) throw std::runtime_error("Clone lost loaded citation history");
        }
#ifdef ALPS_HAVE_MPI
        world.barrier();
#endif
        if (rank == 0) {
          for (int r = 0; r < size; ++r) {
            auto records = alps::read_citations(directory / ("copy." + std::to_string(r) + ".h5"));
            if (records.size() != 2) throw std::runtime_error("Rank checkpoint missing citations");
          }
        }
#ifdef ALPS_HAVE_MPI
        world.barrier(); // finish cross-rank reads before reusing the files
#endif
        // Fresh clones reusing an old filename do not inherit its old analysis.
#ifdef ALPS_HAVE_MPI
        alps::clone_mpi fresh(world, world, directory, opt, message);
#else
        alps::clone fresh(directory, opt, 0, 0, p, "checkpoint", true);
#endif
        {
          alps::hdf5::archive ar((directory / ("copy." + std::to_string(rank) + ".h5")).string(), "a");
          fresh.save(ar);
          auto records = alps::read_citations(ar);
          if (records.size() != 1 || records[0].activity != "calculation")
            throw std::runtime_error("Fresh clone inherited old filename ancestry");
        }
        return 0;
      }
      return alps::parapack::start(argc, argv);
    }
    throw std::invalid_argument("Unknown test mode");
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
