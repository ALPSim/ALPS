// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Writes results_7146b9e1.h5 the way ALPSCore's Monte Carlo programs write
// their output: parameters under /parameters, results under
// /simulation/results. It then reads the file back with ALPSCore's archive
// and prints the stored values as JSON, which become expected.json.
// See README.md for how to build it against ALPSCore.
#include <alps/accumulators.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/params.hpp>

#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::string json(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

std::string json(const std::vector<double>& values) {
    std::string out = "[";
    for (std::size_t i = 0; i < values.size(); ++i)
        out += (i ? ", " : "") + json(values[i]);
    return out + "]";
}

std::string json(const std::vector<std::vector<double>>& rows) {
    std::string out = "[";
    for (std::size_t i = 0; i < rows.size(); ++i)
        out += (i ? ", " : "") + json(rows[i]);
    return out + "]";
}

// One observable's stored count, mean, error and binning analysis.
template <typename T, typename Bins>
void print(alps::hdf5::archive& archive, const std::string& name, bool last) {
    const std::string path = "/simulation/results/" + name;
    unsigned long count;
    T mean;
    archive[path + "/count"] >> count;
    archive[path + "/mean/value"] >> mean;
    std::cout << "  \"" << name << "\": {\"count\": " << count << ", \"mean\": " << json(mean);
    if (archive.is_data(path + "/mean/error")) {
        T error;
        archive[path + "/mean/error"] >> error;
        std::cout << ", \"error\": " << json(error);
    }
    if (archive.is_data(path + "/error_bins")) {
        Bins bins;
        archive[path + "/error_bins"] >> bins;
        std::cout << ", \"error_bins\": " << json(bins);
    }
    std::cout << "}" << (last ? "\n" : ",\n");
}

}

int main(int argc, char** argv) {
    namespace acc = alps::accumulators;
    // Defaults only: ALPSCore's params(argc, argv) would read argv[1] as an INI file.
    alps::params parameters;
    parameters.define<int>("L", 4, "linear size")
              .define<double>("T", 1.5, "temperature")
              .define<std::string>("MODEL", "ising", "model name")
              .define<bool>("MEASURE", true, "whether to measure");

    acc::accumulator_set measurements;
    measurements << acc::MeanAccumulator<double>("Mean")
                 << acc::NoBinningAccumulator<double>("NoBinning")
                 << acc::LogBinningAccumulator<double>("LogBinning")
                 << acc::FullBinningAccumulator<double>("FullBinning")
                 << acc::LogBinningAccumulator<std::vector<double>>("Vector LogBinning")
                 << acc::FullBinningAccumulator<std::vector<double>>("Vector FullBinning");

    // Correlated samples from mt19937's raw output, which the C++ standard fixes.
    std::mt19937 engine(42);
    double x = 0.;
    for (int i = 0; i < 4096; ++i) {
        x = 0.9 * x + engine() / 4294967296.0 - 0.5;
        for (const char* name : {"Mean", "NoBinning", "LogBinning", "FullBinning"})
            measurements[name] << x;
        measurements["Vector LogBinning"] << std::vector<double>{x, 2. * x};
        measurements["Vector FullBinning"] << std::vector<double>{x, 2. * x};
    }

    const std::string filename = argc > 1 ? argv[1] : "results_7146b9e1.h5";
    {
        alps::hdf5::archive archive(filename, "w");
        archive["/parameters"] << parameters;
        archive["/simulation/results"] << acc::result_set(measurements);
    }

    alps::hdf5::archive archive(filename, "r");
    typedef std::vector<double> vector;
    std::cout << "{\n";
    print<double, vector>(archive, "Mean", false);
    print<double, vector>(archive, "NoBinning", false);
    print<double, vector>(archive, "LogBinning", false);
    print<double, vector>(archive, "FullBinning", false);
    print<vector, std::vector<vector>>(archive, "Vector LogBinning", false);
    print<vector, std::vector<vector>>(archive, "Vector FullBinning", true);
    std::cout << "}\n";
}
