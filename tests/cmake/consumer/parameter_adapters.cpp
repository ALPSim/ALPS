// SPDX-License-Identifier: MIT
#include <alps/params.hpp>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/ngs/make_parameters_from_file.hpp>
#include <alps/ngs/make_parameters_from_xml.hpp>
#include <boost/filesystem/operations.hpp>
#include <fstream>
#include <map>
#include <stdexcept>

void require(bool condition) {
    if (!condition) throw std::runtime_error("Legacy parameter adapter contract failed");
}

int main() {
    const boost::filesystem::path text("legacy-params-contract.parm");
    {
        std::ofstream output(text.string());
        output << "// Existing ALPS parameter grammar\n"
               << "LATTICE=\"chain lattice\"; L=10; T=2.25;\n";
    }
    // File input supplies text; the program's definitions give the types.
    auto parameters = alps::make_parameters_from_file(text);
    const std::map<std::string, std::string> supplied = {
        {"L", "10"}, {"LATTICE", "chain lattice"}, {"T", "2.25"}};
    require(parameters.supplied_text() == supplied);
    parameters.define<std::string>("LATTICE", "lattice").define<int>("L", "length");
    require(parameters.ok() && parameters["LATTICE"].as<std::string>() == "chain lattice");
    require(parameters["L"].as<int>() == 10);
    // The bridge passes on supplied text, including parameters not defined.
    auto legacy = alps::make_deprecated_parameters(parameters);
    require(legacy.defined("L") && std::string(legacy["L"]) == "10");
    require(std::string(legacy["LATTICE"]) == "chain lattice" && std::string(legacy["T"]) == "2.25");
    parameters.define<int>("SITES", 1, "sites");
    require(!alps::make_deprecated_parameters(parameters).defined("SITES"));
    {
        std::ofstream output(text.string());
        output << "L=10; garbage @\n";
    }
    bool caught = false;
    try { (void)alps::make_parameters_from_file(text); }
    catch (const std::runtime_error&) { caught = true; }
    require(caught);
    boost::filesystem::remove(text);

    const boost::filesystem::path xml("legacy-params-contract.xml");
    for (bool explicit_seed : {false, true}) {
        {
            std::ofstream output(xml.string());
            output << "<SIMULATION><IGNORED/><PARAMETERS>"
                   << "<PARAMETER name=\"L\">12</PARAMETER>";
            if (explicit_seed) output << "<PARAMETER name=\"SEED\">23</PARAMETER>";
            output << "</PARAMETERS></SIMULATION>\n";
        }
        auto from_xml = alps::make_parameters_from_xml(xml);
        from_xml.define<int>("L", "length").define<int>("SEED", "seed");
        require(from_xml.ok() && from_xml["L"].as<int>() == 12);
        require(from_xml["SEED"].as<int>() == (explicit_seed ? 23 : 0));
    }
    boost::filesystem::remove(xml);
}
