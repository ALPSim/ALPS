// SPDX-License-Identifier: MIT
#include <alps/ngs/params.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <sstream>
#include <stdexcept>

void require(bool condition) {
    if (!condition) throw std::runtime_error("Typed params SDK contract failed");
}

int main() {
    alps::params parameters;
    parameters["count"] = 3;
    parameters["values"] = std::vector<double>{1., 2.};
    parameters["label"] = "sample";
    require(parameters.size() == 3 && !parameters.defined("absent"));
    require((parameters["absent"] | 7) == 7 && parameters.find("absent") == nullptr);
    require(parameters["values"] + std::vector<double>{3., 4.}
            == std::vector<double>({4., 6.}));
    const alps::params copy(parameters);
    require(copy["count"].cast<int>() == 3 && copy.begin()->first == "count");
    bool caught = false;
    try { (void)copy["absent"].cast<int>(); }
    catch (const std::runtime_error&) { caught = true; }
    require(caught);

    const char* filename = "params-component-contract.h5";
    {
        alps::hdf5::archive archive(filename, "w");
        archive["/parameters"] << parameters;
    }
    {
        alps::hdf5::archive archive(filename, "r");
        alps::params restored(archive);
        require(restored["values"].cast<std::vector<double>>() == std::vector<double>({1., 2.}));
        require(restored["label"].cast<std::string>() == "sample");
    }
    boost::filesystem::remove(filename);
    std::stringstream buffer;
    { boost::archive::text_oarchive archive(buffer); archive << copy; }
    alps::params restored;
    { boost::archive::text_iarchive archive(buffer); archive >> restored; }
    require(restored["count"].cast<int>() == 3);
    restored.erase("count");
    require(restored.size() == 2 && !restored.defined("count"));
    std::ostringstream output;
    output << restored;
    require(output.str().find("label = sample") != std::string::npos);
}
