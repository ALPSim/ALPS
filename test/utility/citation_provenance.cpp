// SPDX-License-Identifier: MIT
#include <alps/utility/citation_provenance.hpp>
#include <alps/utility/citations.hpp>
#include <alps/scheduler.h>
#include <alps/mcbase.hpp>
#include <alps/ngs/api.hpp>
#include <alps/hdf5.hpp>
#include <boost/filesystem/operations.hpp>
#include <stdexcept>
#include <iostream>

namespace {
void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}
class tiny_mc : public alps::mcbase {
public:
  explicit tiny_mc(const alps::params& p) : alps::mcbase(p) {
    measurements << alps::accumulator::RealObservable("Value");
  }
  void update() override {}
  void measure() override { measurements["Value"] << 2.; }
  double fraction_completed() const override { return 1.; }
};
}

int main(int argc, char** argv) {
  try {
    if (argc != 2) throw std::runtime_error("Expected an output directory");
    const boost::filesystem::path dir(argv[1]);
    boost::filesystem::create_directories(dir);
    const auto calculation = alps::make_citation_snapshot("looper");
    const auto analysis = alps::make_citation_snapshot("looper", "analysis");
    require(calculation.notice == alps::citation_text("looper"), "Saved guidance differs from CLI");
    const auto file = dir / "roundtrip.h5";
    {
      alps::hdf5::archive ar(file.string(), "w");
      ar["/simulation/results/mean/value"] << 42.;
      ar["/parameters/SEED"] << 17;
      alps::write_citations(ar, alps::citation_history{calculation});
      alps::write_citations(ar, alps::citation_history{calculation});
      require(alps::read_citations(ar).size() == 1, "Repeated checkpoint duplicated metadata");
      alps::write_citations(ar, alps::citation_history{analysis});
      double value = 0;
      ar["/simulation/results/mean/value"] >> value;
      require(value == 42., "Citation writing changed numerical data");
      auto collision = calculation;
      collision.notice += "changed";
      try { alps::write_citations(ar, alps::citation_history{collision}); return 2; }
      catch (const std::runtime_error&) {}
      try { alps::replace_citations(ar, alps::citation_history{collision}); return 5; }
      catch (const std::runtime_error&) {}
      require(alps::read_citations(ar).size() == 2, "Collision changed existing records");
      // A failed append is ignored and repaired on the next write.
      const auto partial = alps::make_citation_snapshot("hirschfye");
      ar["/provenance/alps/citations/records/" + partial.id + "/component"] << partial.component;
      require(alps::read_citations(ar).size() == 2, "Incomplete append exposed a partial record");
      alps::write_citations(ar, alps::citation_history{partial});
      require(alps::read_citations(ar).size() == 3, "Incomplete record was not repaired");
    }
    auto history = alps::read_citations(file);
    require(history.size() == 3, "History lost after reopen");
    // Simulate a result copied to a new file: all old snapshots travel with it.
    {
      alps::hdf5::archive ar((dir / "copied.h5").string(), "w");
      alps::write_citations(ar, history);
    }
    require(alps::read_citations(dir / "copied.h5") == history, "Historical snapshot changed on copy");
    {
      alps::hdf5::archive ar((dir / "legacy.h5").string(), "w");
      ar["/parameters/SEED"] << 5;
      require(alps::read_citations(ar).empty(), "Legacy file acquired invented citations on read");
    }
    {
      alps::hdf5::archive ar((dir / "future.h5").string(), "w");
      ar["/provenance/alps/citations/schema_version"] << 2;
      try { alps::write_citations(ar, "looper"); return 3; }
      catch (const std::runtime_error&) {}
      try { alps::replace_citations(ar, alps::citation_history{calculation}); return 4; }
      catch (const std::runtime_error&) {}
      require(!ar.is_group("/provenance/alps/citations/records"), "Unsupported schema was mutated");
    }
    {
      alps::hdf5::archive ar((dir / "fractional-schema.h5").string(), "w");
      ar["/provenance/alps/citations/schema_version"] << 1.5;
      try { alps::write_citations(ar, "looper"); return 6; }
      catch (const std::runtime_error&) {}
      require(!ar.is_group("/provenance/alps/citations/records"), "Fractional schema was truncated to version one");
    }
    for (const int marker : {-1, 2}) {
      alps::hdf5::archive ar((dir / ("invalid-complete-" + std::to_string(marker) + ".h5")).string(), "w");
      alps::write_citations(ar, "looper");
      ar["/provenance/alps/citations/records/" + calculation.id + "/complete"] << marker;
      try { alps::read_citations(ar); return 7; }
      catch (const std::runtime_error&) {}
    }
    // Framework-only profiles exercise empty algorithm/implementation arrays.
    {
      alps::hdf5::archive ar((dir / "framework.h5").string(), "w");
      alps::write_citations(ar, "framework", "unspecified");
      auto records = alps::read_citations(ar);
      require(records.size() == 1 && records[0].algorithm.empty(), "Empty role arrays did not roundtrip");
    }
    // A real task checkpoint/reload carries the old history through .bak replacement.
    alps::Parameters old_parameters;
    old_parameters["SEED"] = 17;
    {
      alps::scheduler::MCSimulation task(alps::ProcessList(), old_parameters);
      task.set_citation_component("looper");
      task.checkpoint_hdf5(dir / "task.out.xml");
    }
    {
      alps::scheduler::MCSimulation task(alps::ProcessList(), dir / "task.out.h5");
      task.set_citation_component("looper", "analysis");
      task.checkpoint_hdf5(dir / "task.out.xml");
    }
    require(alps::read_citations(dir / "task.out.h5").size() == 2, "Task replacement lost original provenance");
    {
      alps::scheduler::Worker worker(alps::ProcessList(), old_parameters, 0);
      worker.set_citation_component("looper");
      worker.start_worker();
      alps::hdf5::archive ar((dir / "worker.h5").string(), "w");
      worker.save(ar);
    }
    {
      alps::scheduler::Worker worker(alps::ProcessList(), old_parameters, 0);
      alps::hdf5::archive ar((dir / "worker.h5").string(), "a");
      worker.load(ar);
      worker.set_citation_component("looper", "analysis");
      worker.save(ar);
      require(alps::read_citations(ar).size() == 2, "Worker checkpoint lost provenance");
    }
    // NGS checkpoint and result paths preserve explicitly inherited input histories.
    alps::params p;
    p["SEED"] = 42;
    tiny_mc mc(p);
    mc.set_citation_component("interaction");
    mc.inherit_citations(alps::citation_history{calculation});
    mc.measure(); mc.measure();
    mc.save(dir / "mc.h5");
    tiny_mc restored(p);
    restored.load(dir / "mc.h5");
    restored.set_citation_component("interaction", "analysis");
    restored.save(dir / "mc.h5");
    require(alps::read_citations(dir / "mc.h5").size() == 3, "NGS restart lost inherited citations");
    alps::save_results(restored.collect_results(), p, dir / "results.h5", "/simulation/results", restored.citations());
    require(alps::read_citations(dir / "results.h5").size() == 3, "Results lost checkpoint provenance");
    // A fresh calculation reusing a filename must not inherit unrelated old data.
    alps::save_results(mc.collect_results(), p, dir / "results.h5", "/simulation/results");
    auto fresh = alps::read_citations(dir / "results.h5");
    require(fresh.size() == 1 && fresh[0].component == "framework", "Unrelated output filename was treated as ancestry");
    // Loading unrelated/legacy state into a reused object also replaces its ancestry.
    tiny_mc unrelated(p);
    unrelated.save(dir / "legacy_checkpoint.h5");
    {
      alps::hdf5::archive ar((dir / "legacy_checkpoint.h5").string(), "a");
      ar.delete_group("/provenance/alps/citations");
    }
    restored.load(dir / "legacy_checkpoint.h5");
    require(restored.citations().size() == 1, "Loading unrelated state retained the previous object's ancestry");
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
