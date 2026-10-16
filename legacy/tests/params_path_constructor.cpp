// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// ALPS 2.x code constructs parameters straight from a text file. The deprecated
// constructor must keep doing that, with the same result as params_from_file,
// rather than implicitly opening the file as HDF5. Its deprecation warnings are
// intentional: tests/cmake/test_sdk_contracts.py requires them.
#include <alps/ngs/params.hpp>
#include <alps/ngs/params_from_file.hpp>

#include <boost/filesystem.hpp>

#include <fstream>
#include <stdexcept>
#include <type_traits>

static_assert(std::is_constructible<alps::params, boost::filesystem::path>::value,
              "the legacy layer must declare alps::params(path)");

void require(bool condition) {
    if (!condition) throw std::runtime_error("Legacy params(path) contract failed");
}

int main() {
    const boost::filesystem::path text("legacy-params-path.parm");
    {
        std::ofstream output(text.string());
        output << "LATTICE=\"chain lattice\"; L=10; T=2.25;\n";
    }
    const alps::params current = alps::params_from_file(text);
    const alps::params from_path{text};
    const alps::params from_string{text.string()};
    const alps::params from_c_string{text.string().c_str()};
    alps::params assigned;
    assigned = alps::params(text);
    const alps::params* restored_parameters[] = {&from_path, &from_string, &from_c_string, &assigned};
    for (const auto* restored : restored_parameters) {
        require(restored->size() == current.size());
        require(restored->begin()->first == "LATTICE");
        require((*restored)["LATTICE"].cast<std::string>() == "chain lattice");
        require((*restored)["L"].cast<int>() == 10);
        require((*restored)["T"].cast<double>() == 2.25);
        require(restored->find("L")->which() == current.find("L")->which());
    }
    {
        std::ofstream output(text.string());
        output << "L=10; garbage @\n";
    }
    bool caught = false;
    try { alps::params invalid{text}; }
    catch (const std::runtime_error&) { caught = true; }
    require(caught);
    boost::filesystem::remove(text);
}
