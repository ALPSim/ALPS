// SPDX-License-Identifier: MIT
#include <alps/ngs/params_from_file.hpp>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/ngs/make_parameters_from_xml.hpp>
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include <fstream>
#include <stdexcept>

TEST(ParamsLegacyAdapters, TextGrammarAndLegacyConversion) {
    alps::testing::TemporaryDirectory temporary;
    const boost::filesystem::path text((temporary.path() / "legacy-params-contract.parm").string());
    {
        std::ofstream output(text.string());
        output << "// Existing ALPS parameter grammar\n"
               << "LATTICE=\"chain lattice\"; L=10; T=2.25;\n";
    }
    auto parameters = alps::params_from_file(text);
    EXPECT_TRUE(parameters.size() == 3 && parameters.begin()->first == "LATTICE");
    EXPECT_EQ(parameters["LATTICE"].cast<std::string>(), "chain lattice");
    EXPECT_TRUE(parameters["L"].cast<int>() == 10 && parameters["T"].cast<double>() == 2.25);
    // File input stores text values; conversion preserves their existing form.
    EXPECT_EQ(parameters.find("L")->which(), alps::detail::paramvalue_index<std::string>::value);
    auto legacy = alps::make_deprecated_parameters(parameters);
    EXPECT_TRUE(legacy.defined("L"));
    EXPECT_EQ(std::string(legacy["L"]), "10");
    EXPECT_EQ(std::string(legacy["LATTICE"]), "chain lattice");
    {
        std::ofstream output(text.string());
        output << "L=10; garbage @\n";
    }
    EXPECT_THROW({ alps::params_from_file(text); }, std::runtime_error);
    boost::filesystem::remove(text);

}

TEST(ParamsLegacyAdapters, XmlSeedDefaultsAndExplicitValues) {
    alps::testing::TemporaryDirectory temporary;
    const boost::filesystem::path xml((temporary.path() / "legacy-params-contract.xml").string());
    for (bool explicit_seed : {false, true}) {
        {
            std::ofstream output(xml.string());
            output << "<SIMULATION><IGNORED/><PARAMETERS>"
                   << "<PARAMETER name=\"L\">12</PARAMETER>";
            if (explicit_seed) output << "<PARAMETER name=\"SEED\">23</PARAMETER>";
            output << "</PARAMETERS></SIMULATION>\n";
        }
        const auto from_xml = alps::make_parameters_from_xml(xml);
        EXPECT_EQ(from_xml["L"].cast<int>(), 12);
        EXPECT_EQ(from_xml["SEED"].cast<int>(), (explicit_seed ? 23 : 0));
    }
    boost::filesystem::remove(xml);
}
