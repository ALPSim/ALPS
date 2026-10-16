// SPDX-License-Identifier: MIT
#include <alps/mc/api.hpp>
#include <alps/mc/mcbase.hpp>
#include <alps/mc/stop_callback.hpp>
#include <alps/hdf5/archive.hpp>
#include <cstdio>
#include <stdexcept>

namespace {
class counting_sim : public alps::mcbase {
  public:
    explicit counting_sim(parameters_type const & parameters) : alps::mcbase(parameters) {
        measurements << alps::accumulators::FullBinningAccumulator<double>("x");
    }
    void update() override { value_ = random(); }
    void measure() override { measurements["x"] << value_; ++sweeps_; }
    double fraction_completed() const override { return sweeps_ / 1024.; }

  private:
    double value_ = 0.;
    int sweeps_ = 0;
};
} // namespace

int main() {
    alps::params parameters;
    counting_sim::define_parameters(parameters);
    counting_sim simulation(parameters);
    if (!simulation.run(alps::stop_callback(60)))
        throw std::runtime_error("the simulation stopped before completing");
    const alps::results_type<counting_sim>::type results = alps::collect_results(simulation);
    simulation.save("mc_contract.h5");
    counting_sim restored(parameters);
    restored.load("mc_contract.h5");
    std::remove("mc_contract.h5");
    const alps::results_type<counting_sim>::type reloaded = alps::collect_results(restored);
    if (results["x"].count() != 1024 || reloaded["x"].count() != 1024
        || reloaded["x"].mean<double>() != results["x"].mean<double>())
        throw std::runtime_error("mc checkpoint round trip failed");
}
