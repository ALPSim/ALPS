// SPDX-License-Identifier: MIT
#include <alps/params.hpp>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/ngs/make_parameters_from_file.hpp>
#include <alps/ngs/make_parameters_from_xml.hpp>
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include <boost/filesystem/operations.hpp>
#include <fstream>
#include <map>
#include <stdexcept>

TEST(ParamsLegacyAdapters, TextGrammarAndLegacyConversion) {
    alps::testing::TemporaryDirectory temporary;
    const boost::filesystem::path text((temporary.path() / "legacy-params-contract.parm").string());
    {
        std::ofstream output(text.string());
        output << "// Existing ALPS parameter grammar\n"
               << "LATTICE=\"chain lattice\"; L=10; T=2.25;\n";
    }
    // File input supplies text; the program's definitions give the types.
    auto parameters = alps::make_parameters_from_file(text);
    const std::map<std::string, std::string> supplied = {
        {"L", "10"}, {"LATTICE", "chain lattice"}, {"T", "2.25"}};
    EXPECT_EQ(parameters.supplied_text(), supplied);
    parameters.define<std::string>("LATTICE", "lattice")
              .define<int>("L", "length")
              .define<double>("T", "temperature");
    ASSERT_TRUE(parameters.ok());
    EXPECT_EQ(parameters["LATTICE"].as<std::string>(), "chain lattice");
    EXPECT_TRUE(parameters["L"].as<int>() == 10 && parameters["T"].as<double>() == 2.25);
    auto legacy = alps::make_deprecated_parameters(parameters);
    EXPECT_TRUE(legacy.defined("L"));
    EXPECT_EQ(std::string(legacy["L"]), "10");
    EXPECT_EQ(std::string(legacy["LATTICE"]), "chain lattice");
    {
        std::ofstream output(text.string());
        output << "L=10; garbage @\n";
    }
    EXPECT_THROW({ (void)alps::make_parameters_from_file(text); }, std::runtime_error);
    boost::filesystem::remove(text);
}

TEST(ParamsLegacyAdapters, TextExpressionsAreEvaluated) {
    alps::testing::TemporaryDirectory temporary;
    const boost::filesystem::path text((temporary.path() / "legacy-params-expressions.parm").string());
    {
        std::ofstream output(text.string());
        output << "LATTICE=\"chain lattice\"; J=2; T=J/4; L=2*5; H=0.123456789123;\n"
               << "MODEL=\"spin\"; SIZES=\"1,2,3\"\n";
    }
    auto parameters = alps::make_parameters_from_file(text);
    const std::map<std::string, std::string> supplied = {
        {"H", "0.123456789123"}, {"J", "2"}, {"L", "10"}, {"LATTICE", "chain lattice"},
        {"MODEL", "spin"}, {"SIZES", "1,2,3"}, {"T", "0.5"}};
    EXPECT_EQ(parameters.supplied_text(), supplied);
    parameters.define<int>("L", "length").define<double>("T", "temperature");
    ASSERT_TRUE(parameters.ok());
    EXPECT_EQ(parameters["L"].as<int>(), 10);
    EXPECT_EQ(parameters["T"].as<double>(), 0.5);
    boost::filesystem::remove(text);
}

TEST(ParamsLegacyAdapters, TextBlocks) {
    alps::testing::TemporaryDirectory temporary;
    const boost::filesystem::path text((temporary.path() / "legacy-params-blocks.parm").string());
    {
        std::ofstream output(text.string());
        output << "LATTICE=\"chain lattice\"\nL=10\nJ=2\n{ T = J/4 }\n";
    }
    // One block is the parameter set, with the keys set before it.
    auto parameters = alps::make_parameters_from_file(text);
    const std::map<std::string, std::string> supplied = {
        {"J", "2"}, {"L", "10"}, {"LATTICE", "chain lattice"}, {"T", "0.5"}};
    EXPECT_EQ(parameters.supplied_text(), supplied);
    {
        std::ofstream output(text.string());
        output << "L=10\n{ T = 0.5 }\n{ T = 1.0 }\n";
    }
    try {
        (void)alps::make_parameters_from_file(text);
        FAIL() << "two blocks were accepted";
    } catch (const std::runtime_error& error) {
        EXPECT_NE(std::string(error.what()).find("2 parameter sets"), std::string::npos)
            << error.what();
    }
    boost::filesystem::remove(text);
}

// The bridge passes on what was supplied, as text, and leaves defaults to the
// code that reads the alps::Parameters.
TEST(ParamsLegacyAdapters, BridgePassesSuppliedTextOnly) {
    alps::params parameters;
    parameters.supply("U", "2").supply("H", "0.123456789").supply("FLAVORS", "4");
    parameters.define<double>("U", "interaction").define<int>("SITES", 1, "sites");
    const auto legacy = alps::make_deprecated_parameters(parameters);
    EXPECT_EQ(std::string(legacy["U"]), "2");
    EXPECT_EQ(std::string(legacy["H"]), "0.123456789");
    EXPECT_EQ(std::string(legacy["FLAVORS"]), "4");
    EXPECT_FALSE(legacy.defined("SITES"));
    EXPECT_EQ(legacy.size(), 3u);
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
        auto from_xml = alps::make_parameters_from_xml(xml);
        from_xml.define<int>("L", "length").define<int>("SEED", "seed");
        ASSERT_TRUE(from_xml.ok());
        EXPECT_EQ(from_xml["L"].as<int>(), 12);
        EXPECT_EQ(from_xml["SEED"].as<int>(), (explicit_seed ? 23 : 0));
    }
    boost::filesystem::remove(xml);
}

TEST(ParamsLegacyAdapters, XmlExpressionsAreEvaluated) {
    alps::testing::TemporaryDirectory temporary;
    const boost::filesystem::path xml((temporary.path() / "legacy-params-expressions.xml").string());
    {
        std::ofstream output(xml.string());
        output << "<SIMULATION><PARAMETERS>"
               << "<PARAMETER name=\"LATTICE\">chain lattice</PARAMETER>"
               << "<PARAMETER name=\"J\">2</PARAMETER>"
               << "<PARAMETER name=\"T\">J/4</PARAMETER>"
               << "</PARAMETERS></SIMULATION>\n";
    }
    auto from_xml = alps::make_parameters_from_xml(xml);
    const std::map<std::string, std::string> supplied = {
        {"J", "2"}, {"LATTICE", "chain lattice"}, {"SEED", "0"}, {"T", "0.5"}};
    EXPECT_EQ(from_xml.supplied_text(), supplied);
    boost::filesystem::remove(xml);
}
