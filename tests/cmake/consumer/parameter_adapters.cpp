// SPDX-License-Identifier: MIT
#include <alps/ngs/params_from_file.hpp>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/ngs/make_parameters_from_xml.hpp>
#include <fstream>
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
    auto parameters = alps::params_from_file(text);
    require(parameters.size() == 3 && parameters.begin()->first == "LATTICE");
    require(parameters["LATTICE"].cast<std::string>() == "chain lattice");
    require(parameters["L"].cast<int>() == 10 && parameters["T"].cast<double>() == 2.25);
    // File input stores text values; conversion preserves their existing form.
    require(parameters.find("L")->which() == alps::detail::paramvalue_index<std::string>::value);
    auto legacy = alps::make_deprecated_parameters(parameters);
    require(legacy.defined("L") && std::string(legacy["L"]) == "10");
    require(std::string(legacy["LATTICE"]) == "chain lattice");
    {
        std::ofstream output(text.string());
        output << "L=10; garbage @\n";
    }
    bool caught = false;
    try { alps::params_from_file(text); }
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
        const auto from_xml = alps::make_parameters_from_xml(xml);
        require(from_xml["L"].cast<int>() == 12);
        require(from_xml["SEED"].cast<int>() == (explicit_seed ? 23 : 0));
    }
    boost::filesystem::remove(xml);
}
