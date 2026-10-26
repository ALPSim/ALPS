// SPDX-License-Identifier: MIT
#include <alps/params.hpp>
#include <alps/params/text_values.hpp>
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

void require(bool condition) {
    if (!condition) throw std::runtime_error("Typed params SDK contract failed");
}

int main() {
    alps::params parameters;
    parameters.supply("count", "3").supply("values", "1,2").supply("unused", "x");
    parameters.define<int>("count", "a count")
              .define<std::vector<double>>("values", "some values")
              .define<std::string>("label", "sample", "a label");
    require(parameters.ok() && parameters["count"].as<int>() == 3);
    require(parameters["values"].as<std::vector<double>>() == std::vector<double>({1., 2.}));
    require(parameters.defaulted("label") && parameters.supplied("count"));
    std::ostringstream unused;
    require(parameters.has_unused(unused) && unused.str().find("unused") != std::string::npos);

    alps::params malformed;
    malformed.supply("count", "2.5");
    malformed.define<int>("count", 1, "a count");
    std::ostringstream missing;
    require(malformed.has_missing(missing) && missing.str().find("count") != std::string::npos);

    const char* filename = "params-component-contract.h5";
    {
        alps::hdf5::archive archive(filename, "w");
        archive["/parameters"] << parameters;
        archive["/plain/count"] << 3;
        archive["/plain/beta"] << 2.0;
    }
    {
        alps::hdf5::archive archive(filename, "r");
        alps::params restored;
        archive["/parameters"] >> restored;
        require(restored == parameters);
        alps::params plain;
        archive["/plain"] >> plain;
        plain.define<int>("count", "a count").define<double>("beta", "inverse temperature");
        require(plain.ok() && plain["count"].as<int>() == 3 && plain["beta"].as<double>() == 2.);
    }
    std::remove(filename);
    require(alps::params_ns::to_text(24.0) == "24" && alps::params_ns::to_text(0.1) == "0.1");
    std::ostringstream output;
    output << parameters;
    require(output.str().find("label") != std::string::npos);
}
