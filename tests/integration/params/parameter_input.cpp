// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/parameter_input.hpp>
#include <alps/expression.h>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include <hdf5.h>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using path = std::filesystem::path;
using value_map = std::map<std::string, std::string>;

void write_text(path const & file, std::string const & text) {
    std::ofstream(file, std::ios::binary) << text;
}

std::string read_text(path const & file) {
    std::ifstream in(file, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

value_map values(alps::Parameters const & parameters) {
    value_map result;
    for (alps::Parameter const & p : parameters)
        result[p.key()] = p.value().c_str();
    return result;
}

alps::Parameters make_parameters(value_map const & entries) {
    alps::Parameters parameters;
    for (auto const & entry : entries)
        parameters[entry.first] = entry.second;
    return parameters;
}

alps::Parameters read_text_file(path const & file) {
    return alps::read_parameter_file(file.string());
}

// Text values that are not numbers and must survive any round trip
// character for character: solver names, file paths, lattice names,
// symbolic expressions and measurement specifications.
value_map const text_values{
    {"SOLVER", "Hirsch-Fye"},
    {"EXTERNAL_SOLVER", "hybridization"},
    {"LATTICE", "chain lattice"},
    {"U_MATRIX", "u-matrix.dat"},
    {"G0TAU_INPUT", "../g0/tau-input.dat"},
    {"LATTICE_LIBRARY", "/opt/alps/lib/xml/lattices.xml"},
    {"SHIFT", "U/2+JH"},
    {"PHASE", "x*cos(2*pi)"},
    {"MEASURE_CORRELATIONS[Spin]", "Sz:Sz"},
    {"FLAVOR_LIST", "up,down"},
    {"NOTE", "it's"},
    {"QUOTE", "say \"hi\""},
    {"EMPTY", ""},
};

// Numbers as users write them; the text format keeps them as written.
value_map const number_values{
    {"L", "16"},
    {"NEGATIVE", "-3"},
    {"BETA", "6"},
    {"t", "0.707106781186547"},
    {"TOLERANCE", "1e-3"},
    {"BIG", "-2.5E+2"},
    {"THIRD", "1/3"},
    {"U", "3"},
    {"MU", "U/2"},
    {"N", "2*L"},
};

}

// --------------------------------------------------------------------------
// Reading ALPS parameter text

TEST(ParameterInputText, AcceptsEverySeparator) {
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "parm", "A = 1; B = 2, C = 3\nD = 4\n\n\nE = 5;\n");
    EXPECT_EQ(values(read_text_file(dir.path() / "parm")),
              (value_map{{"A", "1"}, {"B", "2"}, {"C", "3"}, {"D", "4"}, {"E", "5"}}));
}

TEST(ParameterInputText, AcceptsWindowsLineEndings) {
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "parm", "A = 1\r\nB = \"two words\"\r\n");
    EXPECT_EQ(values(read_text_file(dir.path() / "parm")),
              (value_map{{"A", "1"}, {"B", "two words"}}));
}

TEST(ParameterInputText, IgnoresComments) {
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "parm",
               "// DMFT input\n"
               "L = 4; // trailing comment\n"
               "/* block\n   comment */ BETA = 2;\n");
    EXPECT_EQ(values(read_text_file(dir.path() / "parm")),
              (value_map{{"L", "4"}, {"BETA", "2"}}));
}

TEST(ParameterInputText, KeepsQuotedValuesVerbatim) {
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "parm",
               "LATTICE = \"chain lattice\"\n"
               "MODEL = 'spin'\n"
               "PADDED = ' two spaces  '\n"
               "NOTE = \"it's\"\n"
               "QUOTE = 'say \"hi\"'\n"
               "EMPTY = \"\"\n");
    EXPECT_EQ(values(read_text_file(dir.path() / "parm")),
              (value_map{{"LATTICE", "chain lattice"}, {"MODEL", "spin"},
                         {"PADDED", " two spaces  "}, {"NOTE", "it's"},
                         {"QUOTE", "say \"hi\""}, {"EMPTY", ""}}));
}

TEST(ParameterInputText, KeepsNumbersAsWritten) {
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "parm",
               "L = 16; NEGATIVE = -3; t = 0.707106781186547\n"
               "TOLERANCE = 1e-3; BIG = -2.5E+2\n");
    alps::Parameters parameters = read_text_file(dir.path() / "parm");
    EXPECT_EQ(values(parameters),
              (value_map{{"L", "16"}, {"NEGATIVE", "-3"}, {"t", "0.707106781186547"},
                         {"TOLERANCE", "1e-3"}, {"BIG", "-2.5E+2"}}));
    EXPECT_EQ(static_cast<int>(parameters["L"]), 16);
    EXPECT_EQ(static_cast<int>(parameters["NEGATIVE"]), -3);
    EXPECT_EQ(static_cast<double>(parameters["t"]), 0.707106781186547);
    EXPECT_EQ(static_cast<double>(parameters["TOLERANCE"]), 1e-3);
    EXPECT_EQ(static_cast<double>(parameters["BIG"]), -250.);
}

TEST(ParameterInputText, KeepsExpressionsUnevaluated) {
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "parm", "U = 3; MU = U/2; N = 2*L; L = 8\n");
    alps::Parameters parameters = read_text_file(dir.path() / "parm");
    EXPECT_EQ(std::string(parameters["MU"].c_str()), "U/2");
    EXPECT_EQ(alps::evaluate<double>(parameters["MU"], parameters), 1.5);
    EXPECT_EQ(alps::evaluate<double>(parameters["N"], parameters), 16.);
}

TEST(ParameterInputText, AcceptsKeysWithPrimesHashesAndBrackets) {
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "parm",
               "J' = 1; T#1 = 2\n"
               "MEASURE_AVERAGE[Local magnetization] = Sz\n"
               "U_MATRIX = u-matrix.dat; G0TAU_INPUT = ../g0/tau.dat\n");
    EXPECT_EQ(values(read_text_file(dir.path() / "parm")),
              (value_map{{"J'", "1"}, {"T#1", "2"},
                         {"MEASURE_AVERAGE[Local magnetization]", "Sz"},
                         {"U_MATRIX", "u-matrix.dat"}, {"G0TAU_INPUT", "../g0/tau.dat"}}));
}

TEST(ParameterInputText, LastDuplicateWins) {
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "parm", "BETA = 1\nL = 4\nBETA = 2\n");
    alps::Parameters parameters = read_text_file(dir.path() / "parm");
    EXPECT_EQ(parameters.size(), 2u);
    EXPECT_EQ(std::string(parameters["BETA"].c_str()), "2");
}

TEST(ParameterInputText, ReplacesEnvironmentVariables) {
    alps::testing::TemporaryDirectory dir;
    ASSERT_EQ(setenv("ALPS_PARAMETER_INPUT_TEST_DIR", "/data/run1", 1), 0);
    write_text(dir.path() / "parm", "G0OMEGA_INPUT = ${ALPS_PARAMETER_INPUT_TEST_DIR}/g0.dat\n");
    EXPECT_EQ(std::string(read_text_file(dir.path() / "parm")["G0OMEGA_INPUT"].c_str()),
              "/data/run1/g0.dat");
}

// The text grammar requires at least one parameter.
TEST(ParameterInputText, RejectsFileWithoutParameters) {
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "empty", "");
    write_text(dir.path() / "blank", "\n\n   \n");
    write_text(dir.path() / "comments", "// nothing\n/* here */\n");
    for (char const * name : {"empty", "blank", "comments"}) {
        EXPECT_THROW(read_text_file(dir.path() / name), std::runtime_error) << name;
        EXPECT_THROW(alps::read_params_file((dir.path() / name).string()), std::runtime_error) << name;
    }
}

TEST(ParameterInputText, AcceptsEmptyUnquotedValue) {
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "parm", "L = 4\nOUTPUT_PREFIX =\nBETA = 2\n");
    EXPECT_EQ(values(read_text_file(dir.path() / "parm")),
              (value_map{{"L", "4"}, {"OUTPUT_PREFIX", ""}, {"BETA", "2"}}));
}

TEST(ParameterInputText, RejectsMalformedText) {
    alps::testing::TemporaryDirectory dir;
    for (std::string const text : {"= 3", "L 3", "L = two words", "L = \"open", "1L = 3", "L = 4; ; BETA = 2"}) {
        write_text(dir.path() / "parm", text);
        EXPECT_THROW(read_text_file(dir.path() / "parm"), std::runtime_error) << text;
        EXPECT_THROW(alps::read_params_file((dir.path() / "parm").string()), std::runtime_error) << text;
    }
}

// The format follows the contents, so a text file named like HDF5 is text.
TEST(ParameterInputText, TextNamedH5IsReadAsText) {
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "parm.h5", "L = 4\n");
    EXPECT_EQ(values(read_text_file(dir.path() / "parm.h5")), (value_map{{"L", "4"}}));
}

TEST(ParameterInputText, TypedReaderConvertsOnRequest) {
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "parm",
               "L = 16; t = 0.707106781186547; TOLERANCE = 1e-3; SOLVER = \"Hirsch-Fye\"\n");
    alps::params parameters = alps::read_params_file((dir.path() / "parm").string());
    EXPECT_EQ(parameters["L"].cast<int>(), 16);
    EXPECT_EQ(parameters["t"].cast<double>(), 0.707106781186547);
    EXPECT_EQ(parameters["TOLERANCE"].cast<double>(), 1e-3);
    EXPECT_EQ(parameters["SOLVER"].cast<std::string>(), "Hirsch-Fye");
}

// --------------------------------------------------------------------------
// Reading HDF5

TEST(ParameterInputHdf5, TypedReaderKeepsStoredTypes) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "solver.in.h5";
    {
        alps::hdf5::archive ar(file.string(), "w");
        ar["/parameters/L"] << 16;
        ar["/parameters/NEGATIVE"] << -3;
        ar["/parameters/t"] << 0.707106781186547;
        ar["/parameters/TOLERANCE"] << 1e-3;
        ar["/parameters/SOLVER"] << std::string("Hirsch-Fye");
        ar["/parameters/U_MATRIX"] << std::string("u-matrix.dat");
        ar["/parameters/WEIGHTS"] << std::vector<double>{0.25, 0.75};
    }
    alps::params parameters = alps::read_params_file(file.string());
    EXPECT_EQ(parameters.size(), 7u);
    EXPECT_EQ(parameters["L"].cast<int>(), 16);
    EXPECT_EQ(parameters["NEGATIVE"].cast<int>(), -3);
    EXPECT_EQ(parameters["t"].cast<double>(), 0.707106781186547);
    EXPECT_EQ(parameters["TOLERANCE"].cast<double>(), 1e-3);
    EXPECT_EQ(parameters["SOLVER"].cast<std::string>(), "Hirsch-Fye");
    EXPECT_EQ(parameters["U_MATRIX"].cast<std::string>(), "u-matrix.dat");
    EXPECT_EQ(parameters["WEIGHTS"].cast<std::vector<double>>(), (std::vector<double>{0.25, 0.75}));
}

TEST(ParameterInputHdf5, LegacyReaderKeepsScalarValues) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "solver.in.h5";
    {
        alps::hdf5::archive ar(file.string(), "w");
        ar["/parameters/L"] << 16;
        ar["/parameters/t"] << 0.707106781186547;
        ar["/parameters/TOLERANCE"] << 1e-3;
        ar["/parameters/SOLVER"] << std::string("Hirsch-Fye");
    }
    alps::Parameters parameters = alps::read_parameter_file(file.string());
    EXPECT_EQ(parameters.size(), 4u);
    EXPECT_EQ(static_cast<int>(parameters["L"]), 16);
    EXPECT_EQ(static_cast<double>(parameters["t"]), 0.707106781186547);
    EXPECT_EQ(static_cast<double>(parameters["TOLERANCE"]), 1e-3);
    EXPECT_EQ(std::string(parameters["SOLVER"].c_str()), "Hirsch-Fye");
}

// Solver input files also carry Green's functions; only /parameters is read.
TEST(ParameterInputHdf5, ReadsOnlyTheParametersGroup) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "solver.in.h5";
    {
        alps::hdf5::archive ar(file.string(), "w");
        ar["/G0"] << std::vector<double>{1., 2., 3.};
        ar["/simulation/results/sign"] << 1.;
        ar["/parameters/L"] << 4;
    }
    EXPECT_EQ(values(alps::read_parameter_file(file.string())), (value_map{{"L", "4"}}));
    alps::params parameters = alps::read_params_file(file.string());
    EXPECT_EQ(parameters.size(), 1u);
    EXPECT_EQ(parameters["L"].cast<int>(), 4);
}

TEST(ParameterInputHdf5, EmptyParametersGroupHasNoParameters) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "solver.in.h5";
    {
        alps::hdf5::archive ar(file.string(), "w");
        ar.create_group("/parameters");
    }
    EXPECT_EQ(alps::read_parameter_file(file.string()).size(), 0u);
    EXPECT_EQ(alps::read_params_file(file.string()).size(), 0u);
}

TEST(ParameterInputHdf5, RejectsFileWithoutParametersGroup) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "results.h5";
    {
        alps::hdf5::archive ar(file.string(), "w");
        ar["/G0"] << std::vector<double>{1., 2.};
    }
    EXPECT_THROW(alps::read_parameter_file(file.string()), std::runtime_error);
    EXPECT_THROW(alps::read_params_file(file.string()), std::runtime_error);
}

// HDF5 is recognized by content, so the name does not have to end in .h5.
TEST(ParameterInputHdf5, RecognizedWhateverTheName) {
    alps::testing::TemporaryDirectory dir;
    for (char const * name : {"solver.in.h5", "solver.input", "parm", "parm.txt", "parm.toml"}) {
        auto const file = dir.path() / name;
        {
            alps::hdf5::archive ar(file.string(), "w");
            ar["/parameters/L"] << 16;
        }
        EXPECT_EQ(values(alps::read_parameter_file(file.string())), (value_map{{"L", "16"}})) << name;
        EXPECT_EQ(alps::read_params_file(file.string())["L"].cast<int>(), 16) << name;
    }
}

// A user block moves the HDF5 signature from offset 0 to 512, 1024, ...
TEST(ParameterInputHdf5, RecognizedWithUserBlock) {
    alps::testing::TemporaryDirectory dir;
    for (hsize_t block : {512, 4096}) {
        auto const file = dir.path() / ("userblock" + std::to_string(block));
        hid_t properties = H5Pcreate(H5P_FILE_CREATE);
        ASSERT_GE(H5Pset_userblock(properties, block), 0);
        hid_t id = H5Fcreate(file.string().c_str(), H5F_ACC_TRUNC, properties, H5P_DEFAULT);
        ASSERT_GE(id, 0);
        H5Fclose(id);
        H5Pclose(properties);
        {
            alps::hdf5::archive ar(file.string(), "a");
            ar["/parameters/L"] << 16;
        }
        ASSERT_NE(read_text(file).substr(0, 4), "\211HDF");
        EXPECT_EQ(values(alps::read_parameter_file(file.string())), (value_map{{"L", "16"}})) << block;
        EXPECT_EQ(alps::read_params_file(file.string())["L"].cast<int>(), 16) << block;
    }
}

// --------------------------------------------------------------------------
// Missing and unreadable input

TEST(ParameterInputErrors, MissingFileThrows) {
    alps::testing::TemporaryDirectory dir;
    auto const file = (dir.path() / "missing").string();
    EXPECT_THROW(alps::read_parameter_file(file), std::runtime_error);
    EXPECT_THROW(alps::read_params_file(file), std::runtime_error);
}

TEST(ParameterInputErrors, DirectoryThrows) {
    alps::testing::TemporaryDirectory dir;
    EXPECT_THROW(alps::read_parameter_file(dir.path().string()), std::runtime_error);
    EXPECT_THROW(alps::read_params_file(dir.path().string()), std::runtime_error);
}

TEST(ParameterInputErrors, ErrorNamesTheFile) {
    alps::testing::TemporaryDirectory dir;
    auto const file = (dir.path() / "missing-parm").string();
    try {
        alps::read_parameter_file(file);
        FAIL() << "expected an exception";
    } catch (std::runtime_error const & error) {
        EXPECT_NE(std::string(error.what()).find("missing-parm"), std::string::npos) << error.what();
    }
}

TEST(ParameterInputErrors, UnwritableLocationThrows) {
    alps::testing::TemporaryDirectory dir;
    alps::Parameters parameters = make_parameters({{"L", "4"}});
    EXPECT_THROW(alps::write_parameter_file((dir.path() / "no-such-dir" / "parm").string(), parameters),
                 std::runtime_error);
    EXPECT_ANY_THROW(alps::write_parameter_file((dir.path() / "no-such-dir" / "parm.h5").string(), parameters));
}

// --------------------------------------------------------------------------
// Writing

TEST(ParameterInputWrite, TextKeepsSimpleValuesUnquoted) {
    alps::testing::TemporaryDirectory dir;
    alps::Parameters parameters;
    parameters["L"] = 16;
    parameters["SOLVER"] = "Hirsch-Fye";
    parameters["LATTICE"] = "chain lattice";
    parameters["MU"] = "U/2";
    alps::write_parameter_file((dir.path() / "parm").string(), parameters);
    EXPECT_EQ(read_text(dir.path() / "parm"),
              "L = 16;\nSOLVER = Hirsch-Fye;\nLATTICE = \"chain lattice\";\nMU = U/2;\n");
}

TEST(ParameterInputWrite, TextRoundTripKeepsEveryValue) {
    alps::testing::TemporaryDirectory dir;
    for (value_map const & entries : {text_values, number_values}) {
        auto const file = dir.path() / "parm";
        alps::write_parameter_file(file.string(), make_parameters(entries));
        EXPECT_EQ(values(read_text_file(file)), entries) << read_text(file);
    }
}

TEST(ParameterInputWrite, TextReplacesPreviousContents) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "parm";
    write_text(file, "L = 4; BETA = 2; SOLVER = hybridization; MAX_TIME = 300\n");
    alps::write_parameter_file(file.string(), make_parameters({{"L", "8"}}));
    EXPECT_EQ(values(read_text_file(file)), (value_map{{"L", "8"}}));
}

TEST(ParameterInputWrite, Hdf5KeepsTextValuesAsWritten) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "parm.h5";
    alps::write_parameter_file(file.string(), make_parameters(text_values));
    EXPECT_EQ(values(alps::read_parameter_file(file.string())), text_values);
    alps::hdf5::archive ar(file.string(), "r");
    for (auto const & entry : text_values) {
        std::string stored;
        ar["/parameters/" + entry.first] >> stored;
        EXPECT_EQ(stored, entry.second) << entry.first;
    }
}

TEST(ParameterInputWrite, Hdf5StoresNumbersWithTheirType) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "parm.h5";
    alps::write_parameter_file(file.string(), make_parameters(number_values));
    {
        alps::hdf5::archive ar(file.string(), "r");
        for (char const * key : {"L", "NEGATIVE", "BETA", "BIG", "U", "N"})
            EXPECT_TRUE(ar.is_datatype<int>(std::string("/parameters/") + key)) << key;
        for (char const * key : {"t", "TOLERANCE", "THIRD", "MU"})
            EXPECT_TRUE(ar.is_datatype<double>(std::string("/parameters/") + key)) << key;
    }
    alps::params parameters = alps::read_params_file(file.string());
    EXPECT_EQ(parameters["L"].cast<int>(), 16);
    EXPECT_EQ(parameters["NEGATIVE"].cast<int>(), -3);
    EXPECT_EQ(parameters["BIG"].cast<int>(), -250);
    EXPECT_EQ(parameters["N"].cast<int>(), 32);
    EXPECT_EQ(parameters["t"].cast<double>(), 0.707106781186547);
    EXPECT_EQ(parameters["TOLERANCE"].cast<double>(), 1e-3);
    EXPECT_EQ(parameters["THIRD"].cast<double>(), 1. / 3.);
    EXPECT_EQ(parameters["MU"].cast<double>(), 1.5);
}

// dmft writes its parameters into the results file next to the input.
TEST(ParameterInputWrite, Hdf5KeepsOtherGroups) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "parm.h5";
    {
        alps::hdf5::archive ar(file.string(), "w");
        ar["/simulation/results/G_tau"] << std::vector<double>{0.5, 0.25};
        ar["/parameters/L"] << 4;
    }
    alps::write_parameter_file(file.string(), make_parameters({{"L", "8"}, {"BETA", "2"}}));
    alps::hdf5::archive ar(file.string(), "r");
    std::vector<double> green;
    ar["/simulation/results/G_tau"] >> green;
    EXPECT_EQ(green, (std::vector<double>{0.5, 0.25}));
    int length = 0;
    ar["/parameters/L"] >> length;
    EXPECT_EQ(length, 8);
    EXPECT_TRUE(ar.is_data("/parameters/BETA"));
}

// The output format follows an existing file, so HDF5 stays HDF5.
TEST(ParameterInputWrite, ExistingHdf5StaysHdf5WhateverTheName) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "parm";
    {
        alps::hdf5::archive ar(file.string(), "w");
        ar["/parameters/L"] << 4;
    }
    alps::write_parameter_file(file.string(), make_parameters({{"L", "8"}}));
    EXPECT_EQ(read_text(file).substr(0, 4), "\211HDF");
    EXPECT_EQ(alps::read_params_file(file.string())["L"].cast<int>(), 8);
}

TEST(ParameterInputWrite, NewFileFormatFollowsTheName) {
    alps::testing::TemporaryDirectory dir;
    alps::Parameters parameters = make_parameters({{"L", "8"}});
    alps::write_parameter_file((dir.path() / "parm.h5").string(), parameters);
    alps::write_parameter_file((dir.path() / "parm.txt").string(), parameters);
    EXPECT_EQ(read_text(dir.path() / "parm.h5").substr(0, 4), "\211HDF");
    EXPECT_EQ(read_text(dir.path() / "parm.txt"), "L = 8;\n");
}

// --------------------------------------------------------------------------
// The DMFT pipeline: dmft reads text, writes text and HDF5, and an
// external solver reads the HDF5 copy with typed parameters.

TEST(ParameterInputPipeline, TextToHdf5ToTypedSolverInput) {
    alps::testing::TemporaryDirectory dir;
    auto const input = dir.path() / "parm";
    write_text(input,
               "SOLVER = \"hybridization\"; U_MATRIX = u-matrix.dat\n"
               "BETA = 6; t = 0.707106781186547; U = 3; MU = U/2\n"
               "N = 250; NMATSUBARA = N; FLAVORS = 2\n");
    alps::Parameters parameters = alps::read_parameter_file(input.string());
    parameters["BASENAME"] = input.string();

    alps::write_parameter_file(input.string(), parameters);
    alps::write_parameter_file(input.string() + ".h5", parameters);

    EXPECT_EQ(values(alps::read_parameter_file(input.string())), values(parameters));

    alps::params solver = alps::read_params_file(input.string() + ".h5");
    EXPECT_EQ(solver["SOLVER"].cast<std::string>(), "hybridization");
    EXPECT_EQ(solver["U_MATRIX"].cast<std::string>(), "u-matrix.dat");
    EXPECT_EQ(solver["BASENAME"].cast<std::string>(), input.string());
    EXPECT_EQ(solver["BETA"].cast<int>(), 6);
    EXPECT_EQ(solver["t"].cast<double>(), 0.707106781186547);
    EXPECT_EQ(solver["MU"].cast<double>(), 1.5);
    EXPECT_EQ(solver["NMATSUBARA"].cast<int>(), 250);
}
