// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/ngs/params.hpp>
#include <alps/solvers.hpp>

#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using boost::property_tree::ptree;
using vector = std::vector<double>;
const double pi = boost::math::constants::pi<double>();
int failures = 0;

void require(bool condition, std::string const& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

// The solver still writes deltaOmega.dat in its working directory, even with
// TEXT_OUTPUT disabled. Give this process an isolated, automatically removed tree.
struct workspace {
    boost::filesystem::path previous = boost::filesystem::current_path();
    boost::filesystem::path path = boost::filesystem::temp_directory_path()
        / boost::filesystem::unique_path("alps-maxent-reference-%%%%%%%%");
    workspace() { boost::filesystem::create_directory(path); }
    ~workspace() {
        boost::system::error_code error;
        boost::filesystem::current_path(previous, error);
        boost::filesystem::remove_all(path, error);
    }
};

struct settings {
    double beta, alpha_min, alpha_max, omega_min, omega_max, sigma;
    int ndat, nfreq, n_alpha, max_it;
    explicit settings(ptree const& p)
        : beta(p.get<double>("beta")), alpha_min(p.get<double>("alpha_min")),
          alpha_max(p.get<double>("alpha_max")), omega_min(p.get<double>("omega_min")),
          omega_max(p.get<double>("omega_max")), sigma(p.get<double>("sigma")),
          ndat(p.get<int>("ndat")), nfreq(p.get<int>("nfreq")),
          n_alpha(p.get<int>("n_alpha")), max_it(p.get<int>("max_it")) {}
};

struct reference_case {
    std::string name, dataspace, kernel, grid;
    double norm;
    bool covariance, omit_last, singular_last;
    std::vector<std::pair<double, double>> poles;
    explicit reference_case(ptree const& p)
        : name(p.get<std::string>("name")), dataspace(p.get<std::string>("dataspace")),
          kernel(p.get<std::string>("kernel")), grid(p.get<std::string>("grid")),
          norm(p.get<double>("norm")), covariance(p.get<bool>("covariance")),
          omit_last(p.get<bool>("omit_last")), singular_last(p.get<bool>("singular_last")) {
        for (auto const& entry : p.get_child("poles")) {
            auto item = entry.second.begin();
            if (entry.second.size() != 2) throw std::runtime_error("A pole needs energy and weight");
            double energy = item++->second.get_value<double>();
            poles.emplace_back(energy, item->second.get_value<double>());
        }
    }
};

// Analytic pole Green functions generate the observations independently of the
// solver's grid and kernel. Frequency cases contain one real datum per Matsubara
// frequency and use particle-hole symmetry, not interleaved complex input.
double observation(reference_case const& c, settings const& s, int i) {
    double value = 0;
    for (auto const& pole : c.poles) {
        double energy = pole.first, weight = pole.second;
        if (c.dataspace == "time") {
            double tau = s.beta * i / (s.ndat - 1);
            value -= weight * std::exp(-tau * energy) / (1 + std::exp(-s.beta * energy));
        } else if (c.kernel == "fermionic") {
            double nu = (2 * i + 1) * pi / s.beta;
            value -= weight * nu / (nu * nu + energy * energy);
        } else {
            double nu = 2 * i * pi / s.beta;
            value -= weight * energy * energy / (nu * nu + energy * energy);
        }
    }
    return c.norm * value;
}

struct output {
    vector omega, average, maximum, chi, variance, alpha, probability;
    vector bin_widths;
    vector bosonic_average, bosonic_maximum;
    double norm;
};

output run_case(reference_case const& c, settings const& s,
                boost::filesystem::path const& directory) {
    boost::filesystem::create_directory(directory);
    boost::filesystem::current_path(directory);
    int count = s.ndat - (c.omit_last ? 1 : 0);
    vector data(count), error(count), covariance(count * count, 0);
    alps::params p;
    p["BETA"] = s.beta;
    p["NDAT"] = count;
    p["NFREQ"] = s.nfreq;
    p["N_ALPHA"] = s.n_alpha;
    p["ALPHA_MIN"] = s.alpha_min;
    p["ALPHA_MAX"] = s.alpha_max;
    p["MAX_IT"] = s.max_it;
    p["MAX_TIME"] = 600;
    p["OMEGA_MIN"] = s.omega_min;
    p["OMEGA_MAX"] = s.omega_max;
    p["NORM"] = c.norm;
    p["DEFAULT_MODEL"] = std::string("flat");
    p["DATASPACE"] = c.dataspace;
    p["KERNEL"] = c.kernel;
    p["FREQUENCY_GRID"] = c.grid;
    p["CUT"] = 0.01;
    p["PARTICLE_HOLE_SYMMETRY"] = c.dataspace == "frequency";
    p["TEXT_OUTPUT"] = false;
    p["VERBOSE"] = false;
    p["DATA_IN_HDF5"] = true;
    p["DATA"] = std::string("input.h5");
    if (c.covariance) p["COVARIANCE_MATRIX"] = std::string("input.h5");
    for (int i = 0; i < count; ++i) {
        data[i] = observation(c, s, i);
        error[i] = c.norm * s.sigma * (1.0 + double(i) / s.ndat);
        covariance[i * count + i] = (c.singular_last && i == count - 1) ? 0 : error[i] * error[i];
        p["TAU_" + std::to_string(i)] = s.beta * i / (s.ndat - 1);
    }
    {
        alps::hdf5::archive input("input.h5", "w");
        input["/Data"] << data;
        input["/Error"] << error;
        if (c.covariance) input["/Covariance"] << covariance;
    }
    alps::solvers::maxent(p, "result.out.h5");
    output result;
    result.norm = c.norm;
    {
        // This existing diagnostic is the only public-run observation of the
        // actual quadrature widths. Its text uses six significant digits.
        std::ifstream widths("deltaOmega.dat");
        int index;
        double width, default_weight;
        while (widths >> index >> width >> default_weight) {
            if (index != static_cast<int>(result.bin_widths.size()))
                throw std::runtime_error("Unexpected bin index in deltaOmega.dat");
            result.bin_widths.push_back(width);
        }
    }
    alps::hdf5::archive archive("result.out.h5", "r");
    archive["/spectrum/omega"] >> result.omega;
    archive["/spectrum/average"] >> result.average;
    archive["/spectrum/maximum"] >> result.maximum;
    archive["/spectrum/chi"] >> result.chi;
    archive["/spectrum/variance"] >> result.variance;
    archive["/alpha/values"] >> result.alpha;
    archive["/alpha/probability"] >> result.probability;
    if (c.kernel == "bosonic") {
        archive["/spectrum/bosonic/average"] >> result.bosonic_average;
        archive["/spectrum/bosonic/maximum"] >> result.bosonic_maximum;
    }
    return result;
}

bool finite(vector const& values) {
    return std::all_of(values.begin(), values.end(), [](double x) { return std::isfinite(x); });
}

// Compare arrays with a normwise tolerance: vanishing tails or variances should
// not turn harmless roundoff into a large pointwise relative error.
void compare(vector const& a, vector const& b, double a_scale, double b_scale,
             std::string const& label, double relative = 2e-5, double absolute = 1e-8) {
    if (a.size() != b.size() || !finite(a) || !finite(b)) {
        require(false, label + " has inconsistent shape or nonfinite values");
        return;
    }
    double difference = 0, magnitude = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        difference = std::max(difference, std::abs(a[i] / a_scale - b[i] / b_scale));
        magnitude = std::max({magnitude, std::abs(a[i] / a_scale), std::abs(b[i] / b_scale)});
    }
    double tolerance = absolute + relative * magnitude;
    std::cout << "PAIR " << label << " max_difference=" << difference << " tolerance=" << tolerance << '\n';
    require(difference <= tolerance, label + " equivalence");
}

void validate(reference_case const& c, settings const& s, output const& out) {
    for (auto const* values : {&out.omega, &out.average, &out.maximum, &out.chi, &out.variance}) {
        if (values->size() != static_cast<std::size_t>(s.nfreq) || !finite(*values))
            throw std::runtime_error(c.name + ": spectral output has invalid size or nonfinite values");
    }
    if (out.alpha.size() != static_cast<std::size_t>(s.n_alpha)
        || out.probability.size() != out.alpha.size() || !finite(out.alpha) || !finite(out.probability))
        throw std::runtime_error(c.name + ": invalid alpha output");
    if (out.bin_widths.size() != static_cast<std::size_t>(s.nfreq) || !finite(out.bin_widths))
        throw std::runtime_error(c.name + ": invalid bin-width diagnostic");

    vector edges(s.nfreq + 1), widths(s.nfreq), expected(s.nfreq);
    const double lower = std::tan(pi * (0.01 - 0.5));
    const double upper = std::tan(pi * (0.99 - 0.5));
    for (int j = 0; j <= s.nfreq; ++j) {
        double t = double(j) / s.nfreq;
        if (c.grid == "Lorentzian") t = (std::tan(pi * (0.98 * t + 0.01 - 0.5)) - lower) / (upper - lower);
        edges[j] = s.omega_min + (s.omega_max - s.omega_min) * t;
    }
    for (int j = 0; j < s.nfreq; ++j) {
        widths[j] = edges[j + 1] - edges[j];
        expected[j] = (edges[j + 1] + edges[j]) / 2;
        require(widths[j] > 0, c.name + ": positive bin width");
        require(out.bin_widths[j] > 0
                    && std::abs(out.bin_widths[j] - widths[j]) <= 5.1e-6 * widths[j] + 1e-12,
                c.name + ": quadrature width matches grid within diagnostic text precision");
        require(out.variance[j] >= 0, c.name + ": nonnegative variance");
    }
    compare(out.omega, expected, 1, 1, c.name + "/grid", 1e-13, 1e-13);

    int count = s.ndat - ((c.omit_last || c.singular_last) ? 1 : 0);
    for (auto const& named : std::vector<std::pair<std::string, vector const*>>{
             {"average", &out.average}, {"maximum", &out.maximum}, {"chi", &out.chi}}) {
        vector const& spectrum = *named.second;
        double mass = 0, residual_sum = 0;
        for (int j = 0; j < s.nfreq; ++j) {
            require(spectrum[j] >= 0, c.name + "/" + named.first + ": nonnegative density");
            mass += spectrum[j] * widths[j] / c.norm;
        }
        for (int i = 0; i < count; ++i) {
            double prediction = 0;
            for (int j = 0; j < s.nfreq; ++j) {
                double w = out.omega[j], kernel;
                if (c.dataspace == "time") {
                    double tau = s.beta * i / (s.ndat - 1);
                    kernel = -std::exp(-tau * w) / (1 + std::exp(-s.beta * w));
                } else if (c.kernel == "fermionic") {
                    double nu = (2 * i + 1) * pi / s.beta;
                    kernel = -nu / (nu * nu + w * w);
                } else {
                    double nu = 2 * i * pi / s.beta;
                    kernel = -w * w / (nu * nu + w * w);
                }
                prediction += kernel * spectrum[j] * widths[j];
            }
            double residual = (prediction - observation(c, s, i))
                / (c.norm * s.sigma * (1.0 + double(i) / s.ndat));
            residual_sum += residual * residual;
        }
        double rms = std::sqrt(residual_sum / count);
        std::cout << "CASE " << c.name << '/' << named.first << " normalized_mass=" << mass
                  << " forward_rms_sigma=" << rms << '\n';
        require(std::abs(mass - 1) <= 0.05, c.name + "/" + named.first + ": spectral mass within 5%");
        require(rms <= 5, c.name + "/" + named.first + ": forward RMS within five input sigmas");
    }

    double probability_integral = 0;
    for (int i = 0; i < s.n_alpha; ++i) {
        double expected_alpha = s.alpha_max * std::pow(s.alpha_min / s.alpha_max, double(i) / (s.n_alpha - 1));
        require(std::abs(out.alpha[i] - expected_alpha) <= 1e-12 * expected_alpha,
                c.name + ": geometric alpha grid");
        require(out.probability[i] >= 0, c.name + ": nonnegative alpha probability");
        if (i) {
            require(out.alpha[i - 1] > out.alpha[i], c.name + ": descending alpha grid");
            probability_integral += 0.5 * (out.probability[i - 1] + out.probability[i])
                * (out.alpha[i - 1] - out.alpha[i]);
        }
    }
    std::cout << "CASE " << c.name << " alpha_probability_integral=" << probability_integral << '\n';
    require(std::abs(probability_integral - 1) <= 1e-12, c.name + ": normalized alpha probability density");
    if (c.kernel == "bosonic") {
        vector average(s.nfreq), maximum(s.nfreq);
        for (int j = 0; j < s.nfreq; ++j) {
            average[j] = pi * out.omega[j] * out.average[j];
            maximum[j] = pi * out.omega[j] * out.maximum[j];
        }
        compare(out.bosonic_average, average, 1, 1, c.name + "/bosonic-average", 1e-13, 1e-13);
        compare(out.bosonic_maximum, maximum, 1, 1, c.name + "/bosonic-maximum", 1e-13, 1e-13);
    }
}

void compare_outputs(output const& a, output const& b, std::string const& label) {
    compare(a.omega, b.omega, 1, 1, label + "/omega", 1e-13, 1e-13);
    compare(a.average, b.average, a.norm, b.norm, label + "/average");
    compare(a.maximum, b.maximum, a.norm, b.norm, label + "/maximum");
    compare(a.chi, b.chi, a.norm, b.norm, label + "/chi");
    compare(a.variance, b.variance, a.norm * a.norm, b.norm * b.norm, label + "/variance");
    compare(a.alpha, b.alpha, 1, 1, label + "/alpha", 1e-13, 1e-13);
    compare(a.probability, b.probability, 1, 1, label + "/probability");
}
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <reference_cases.json>\n";
        return 2;
    }
    try {
        ptree specification;
        boost::property_tree::read_json(argv[1], specification);
        settings defaults(specification.get_child("defaults"));
        workspace temporary;
        std::map<std::string, output> results;
        std::cout << std::setprecision(12);
        for (auto const& entry : specification.get_child("cases")) {
            reference_case current(entry.second);
            try {
                auto result = run_case(current, defaults, temporary.path / current.name);
                validate(current, defaults, result);
                results.emplace(current.name, std::move(result));
            } catch (std::exception const& error) {
                require(false, current.name + ": " + error.what());
            }
        }
        for (auto const& pair : std::vector<std::pair<std::string, std::string>>{
                 {"time_linear", "time_scaled"}, {"time_linear", "time_covariance"},
                 {"time_covariance_singular", "time_omitted"}}) {
            if (results.count(pair.first) && results.count(pair.second))
                compare_outputs(results.at(pair.first), results.at(pair.second), pair.first + "=" + pair.second);
            else
                require(false, "Missing output for equivalence pair " + pair.first + "=" + pair.second);
        }
    } catch (std::exception const& error) {
        require(false, error.what());
    }
    if (failures) std::cerr << "MaxEnt reference checks failed: " << failures << '\n';
    return failures ? 1 : 0;
}
