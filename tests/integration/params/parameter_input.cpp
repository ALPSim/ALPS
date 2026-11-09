// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/parameter_input.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/utilities/vectorio.hpp>
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include <hdf5.h>
#include <sys/stat.h>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

// Values read from a file are supplied as text; a program defines each key
// with its type before reading it, as here.
template <typename T>
T typed(alps::params parameters, std::string const & key) {
    parameters.define<T>(key, "");
    if (!parameters.ok())
        throw std::runtime_error("cannot read " + key + " as the requested type");
    return parameters[key].template as<T>();
}

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
    {"BLANK", "  "},
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

// The format follows the contents, so a text file named like HDF5 is text.
TEST(ParameterInputText, TextNamedH5IsReadAsText) {
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "parm.h5", "L = 4\n");
    EXPECT_EQ(values(read_text_file(dir.path() / "parm.h5")), (value_map{{"L", "4"}}));
}

// A solver reading text gets the same typed values it would get from the
// HDF5 file dmft writes for it.
TEST(ParameterInputText, TypedReaderEvaluatesLikeHdf5) {
    alps::testing::TemporaryDirectory dir;
    auto const text = dir.path() / "parm";
    write_text(text,
               "U = 3; MU = U/2; L = 8; N = 2*L; THIRD = 1/3\n"
               "SOLVER = \"Hirsch-Fye\"; SHIFT = \"U/2+JH\"; EMPTY = \"\"\n");
    alps::write_parameter_file((dir.path() / "parm.h5").string(), read_text_file(text));

    alps::params from_text = alps::read_params_file(text.string());
    alps::params from_hdf5 = alps::read_params_file((dir.path() / "parm.h5").string());
    for (char const * key : {"U", "L", "N"})
        EXPECT_EQ(typed<int>(from_text, key), typed<int>(from_hdf5, key)) << key;
    for (char const * key : {"MU", "THIRD"})
        EXPECT_EQ(typed<double>(from_text, key), typed<double>(from_hdf5, key)) << key;
    for (char const * key : {"SOLVER", "SHIFT", "EMPTY"})
        EXPECT_EQ(typed<std::string>(from_text, key), typed<std::string>(from_hdf5, key)) << key;

    EXPECT_EQ(typed<double>(from_text, "MU"), 1.5);
    EXPECT_EQ(typed<int>(from_text, "N"), 16);
    EXPECT_EQ(typed<double>(from_text, "THIRD"), 1. / 3.);
    EXPECT_EQ(typed<std::string>(from_text, "SHIFT"), "U/2+JH");
}

// Input can come from a pipe, e.g. dmft <(generate-parameters).
TEST(ParameterInputText, ReadsFromPipe) {
    alps::testing::TemporaryDirectory dir;
    auto const fifo = dir.path() / "parm.fifo";
    ASSERT_EQ(mkfifo(fifo.c_str(), 0600), 0);
    for (int reader = 0; reader < 2; ++reader) {
        std::thread writer([&] { std::ofstream(fifo) << "L = 4; SOLVER = hybridization\n"; });
        if (reader == 0) {
            EXPECT_EQ(values(alps::read_parameter_file(fifo.string())),
                      (value_map{{"L", "4"}, {"SOLVER", "hybridization"}}));
        } else {
            EXPECT_EQ(typed<int>(alps::read_params_file(fifo.string()), "L"), 4);
        }
        writer.join();
    }
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
    EXPECT_EQ(parameters.supplied_text().size(), 7u);
    EXPECT_EQ(typed<int>(parameters, "L"), 16);
    EXPECT_EQ(typed<int>(parameters, "NEGATIVE"), -3);
    EXPECT_EQ(typed<double>(parameters, "t"), 0.707106781186547);
    EXPECT_EQ(typed<double>(parameters, "TOLERANCE"), 1e-3);
    EXPECT_EQ(typed<std::string>(parameters, "SOLVER"), "Hirsch-Fye");
    EXPECT_EQ(typed<std::string>(parameters, "U_MATRIX"), "u-matrix.dat");
    EXPECT_EQ(typed<std::vector<double>>(parameters, "WEIGHTS"), (std::vector<double>{0.25, 0.75}));
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
    EXPECT_EQ(parameters.supplied_text().size(), 1u);
    EXPECT_EQ(typed<int>(parameters, "L"), 4);
}

// Legacy Parameters hold text, so vectors arrive as the text read_vector reads.
TEST(ParameterInputHdf5, LegacyReaderKeepsVectorsAsText) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "solver.in.h5";
    {
        alps::hdf5::archive ar(file.string(), "w");
        ar["/parameters/L"] << 4;
        ar["/parameters/WEIGHTS"] << std::vector<double>{0.25, 0.75};
        ar["/parameters/SITES"] << std::vector<int>{1, 2, 3};
        ar["/parameters/NAMES"] << std::vector<std::string>{"up", "down"};
    }
    alps::Parameters parameters = alps::read_parameter_file(file.string());
    EXPECT_EQ(std::string(parameters["WEIGHTS"].c_str()), "0.25 0.75");
    EXPECT_EQ(std::string(parameters["SITES"].c_str()), "1 2 3");
    EXPECT_EQ(std::string(parameters["NAMES"].c_str()), "up down");
    EXPECT_EQ(alps::read_vector<std::vector<double>>(std::string(parameters["WEIGHTS"].c_str())),
              (std::vector<double>{0.25, 0.75}));
}

// HDF5 names cannot hold '/', so keys are encoded on disk and decoded on read.
TEST(ParameterInputHdf5, KeysWithSlashesRoundTrip) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "parm.h5";
    alps::write_parameter_file(file.string(),
        make_parameters({{"MEASURE_LOCAL[a/b]", "Sz"}, {"T[x&y]", "1"}, {"L", "4"}}));
    EXPECT_EQ(values(alps::read_parameter_file(file.string())),
              (value_map{{"MEASURE_LOCAL[a/b]", "Sz"}, {"T[x&y]", "1"}, {"L", "4"}}));
    alps::params parameters = alps::read_params_file(file.string());
    EXPECT_EQ(typed<std::string>(parameters, "MEASURE_LOCAL[a/b]"), "Sz");
    EXPECT_EQ(typed<int>(parameters, "T[x&y]"), 1);
    alps::hdf5::archive ar(file.string(), "r");
    EXPECT_EQ(ar.list_children("/parameters").size(), 3u);
}

// Files from other tools may hold a raw '&' in a name.
TEST(ParameterInputHdf5, ReadsRawAmpersandNames) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "parm.h5";
    {
        alps::hdf5::archive ar(file.string(), "w");
        ar["/parameters/A&B"] << 2;
    }
    EXPECT_EQ(values(alps::read_parameter_file(file.string())), (value_map{{"A&B", "2"}}));
    EXPECT_EQ(typed<int>(alps::read_params_file(file.string()), "A&B"), 2);
}

TEST(ParameterInputHdf5, EmptyParametersGroupHasNoParameters) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "solver.in.h5";
    {
        alps::hdf5::archive ar(file.string(), "w");
        ar.create_group("/parameters");
    }
    EXPECT_EQ(alps::read_parameter_file(file.string()).size(), 0u);
    EXPECT_EQ(alps::read_params_file(file.string()).supplied_text().size(), 0u);
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
        EXPECT_EQ(typed<int>(alps::read_params_file(file.string()), "L"), 16) << name;
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
        EXPECT_EQ(typed<int>(alps::read_params_file(file.string()), "L"), 16) << block;
    }
}

// --------------------------------------------------------------------------
// Missing and unreadable input

TEST(ParameterInputErrors, MissingFileThrowsNamingTheFile) {
    alps::testing::TemporaryDirectory dir;
    auto const file = (dir.path() / "missing-parm").string();
    auto const expect_named = [&](auto const & read, char const * reader) {
        try {
            read();
            ADD_FAILURE() << reader << ": expected an exception";
        } catch (std::runtime_error const & error) {
            EXPECT_NE(std::string(error.what()).find("missing-parm"), std::string::npos) << reader << ": " << error.what();
        }
    };
    expect_named([&] { alps::read_parameter_file(file); }, "read_parameter_file");
    expect_named([&] { alps::read_params_file(file); }, "read_params_file");
    expect_named([&] { alps::read_parameter_list(file); }, "read_parameter_list");
}

TEST(ParameterInputErrors, DirectoryThrows) {
    alps::testing::TemporaryDirectory dir;
    EXPECT_THROW(alps::read_parameter_file(dir.path().string()), std::runtime_error);
    EXPECT_THROW(alps::read_params_file(dir.path().string()), std::runtime_error);
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
    EXPECT_EQ(typed<int>(alps::read_params_file(file.string()), "L"), 8);
}

// An existing text file is never replaced by HDF5, whatever its name.
TEST(ParameterInputWrite, ExistingTextStaysTextWhateverTheName) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "parm.h5";
    write_text(file, "L = 4\n");
    alps::write_parameter_file(file.string(), make_parameters({{"L", "8"}}));
    EXPECT_EQ(read_text(file), "L = 8;\n");
}

TEST(ParameterInputWrite, Hdf5RewriteDropsRemovedParameters) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "parm.h5";
    alps::write_parameter_file(file.string(), make_parameters({{"L", "4"}, {"OLD", "1"}}));
    alps::write_parameter_file(file.string(), make_parameters({{"L", "8"}}));
    EXPECT_EQ(values(alps::read_parameter_file(file.string())), (value_map{{"L", "8"}}));
    EXPECT_EQ(alps::read_params_file(file.string()).supplied_text().size(), 1u);
}

// Values from XML or HDF5 input may hold characters the text grammar
// cannot read; writing them fails before the existing file is touched.
TEST(ParameterInputWrite, TextRefusesValuesItCannotReadBack) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "parm";
    for (std::string const value : {"line1\nline2", "a\tb", "it's \"quoted\""}) {
        write_text(file, "L = 4\n");
        alps::Parameters parameters = make_parameters({{"L", "8"}, {"BAD", value}});
        try {
            alps::write_parameter_file(file.string(), parameters);
            ADD_FAILURE() << "expected an exception for " << value;
        } catch (std::runtime_error const & error) {
            EXPECT_NE(std::string(error.what()).find("BAD"), std::string::npos) << error.what();
        }
        EXPECT_EQ(read_text(file), "L = 4\n") << value;
    }
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
    EXPECT_EQ(typed<std::string>(solver, "SOLVER"), "hybridization");
    EXPECT_EQ(typed<std::string>(solver, "U_MATRIX"), "u-matrix.dat");
    EXPECT_EQ(typed<std::string>(solver, "BASENAME"), input.string());
    EXPECT_EQ(typed<int>(solver, "BETA"), 6);
    EXPECT_EQ(typed<double>(solver, "t"), 0.707106781186547);
    EXPECT_EQ(typed<double>(solver, "MU"), 1.5);
    EXPECT_EQ(typed<int>(solver, "NMATSUBARA"), 250);
}

// --------------------------------------------------------------------------
// Several runs and standard input

namespace {

std::vector<value_map> runs(alps::ParameterList const & list) {
    std::vector<value_map> result;
    for (alps::Parameters const & p : list)
        result.push_back(values(p));
    return result;
}

}

TEST(ParameterInputList, EachBlockIsARunWithTheParametersBeforeIt) {
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "parm", "LATTICE = \"chain lattice\"\nL = 8\n{ T = 0.5 }\n{ T = 1.0 }\nL = 16\n{ T = 2.0 }\n");
    EXPECT_EQ(runs(alps::read_parameter_list((dir.path() / "parm").string())), (std::vector<value_map>{
        {{"LATTICE", "chain lattice"}, {"L", "8"}, {"T", "0.5"}},
        {{"LATTICE", "chain lattice"}, {"L", "8"}, {"T", "1.0"}},
        {{"LATTICE", "chain lattice"}, {"L", "16"}, {"T", "2.0"}}}));
}

TEST(ParameterInputList, TextWithoutBlocksIsOneRun) {
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "parm", "L = 8\nT = 0.5\n");
    EXPECT_EQ(runs(alps::read_parameter_list((dir.path() / "parm").string())),
              (std::vector<value_map>{{{"L", "8"}, {"T", "0.5"}}}));
}

// As for a single parameter set, text without parameters is rejected.
TEST(ParameterInputList, RejectsTextWithoutParameters) {
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "parm", "// nothing here\n");
    EXPECT_THROW(alps::read_parameter_list((dir.path() / "parm").string()), std::runtime_error);
}

TEST(ParameterInputList, Hdf5AndXmlHoldOneRun) {
    alps::testing::TemporaryDirectory dir;
    alps::Parameters p = make_parameters({{"L", "8"}, {"SOLVER", "Hirsch-Fye"}});
    alps::write_parameter_file((dir.path() / "parm.h5").string(), p);
    EXPECT_EQ(runs(alps::read_parameter_list((dir.path() / "parm.h5").string())),
              (std::vector<value_map>{values(p)}));
    // The .xml extension is matched case-insensitively.
    write_text(dir.path() / "task.in.XML",
               "<SIMULATION><PARAMETERS><PARAMETER name=\"L\">8</PARAMETER>"
               "<PARAMETER name=\"SOLVER\">Hirsch-Fye</PARAMETER></PARAMETERS></SIMULATION>\n");
    EXPECT_EQ(runs(alps::read_parameter_list((dir.path() / "task.in.XML").string())),
              (std::vector<value_map>{values(p)}));
}

TEST(ParameterInputList, MalformedTextThrows) {
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "parm", "{ T = 0.5 \n");
    EXPECT_THROW(alps::read_parameter_list((dir.path() / "parm").string()), std::runtime_error);
}

TEST(ParameterInputStream, ReadsRunsLikeAFile) {
    std::string const text = "L = 8\n{ T = 0.5 }\n{ T = 1.0 }\n";
    alps::testing::TemporaryDirectory dir;
    write_text(dir.path() / "parm", text);
    std::istringstream in(text);
    EXPECT_EQ(runs(alps::read_parameter_list_text(in)),
              runs(alps::read_parameter_list((dir.path() / "parm").string())));
    std::istringstream flat("L = 8\nT = 0.5\n");
    EXPECT_EQ(runs(alps::read_parameter_list_text(flat)), (std::vector<value_map>{{{"L", "8"}, {"T", "0.5"}}}));
}
