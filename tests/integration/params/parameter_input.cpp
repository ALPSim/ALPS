// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/parameter_input.hpp>
#include <alps/hdf5/archive.hpp>
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {

void write_text(std::filesystem::path const & path, std::string const & text) {
    std::ofstream(path) << text;
}

void write_hdf5(std::filesystem::path const & path) {
    alps::hdf5::archive ar(path.string(), "w");
    ar["/parameters/L"] << 16;
    ar["/parameters/BETA"] << 2.5;
    ar["/parameters/SOLVER"] << std::string("hybridization");
}

}

TEST(ParameterInput, ReadsText) {
    alps::testing::TemporaryDirectory dir;
    auto const file = dir.path() / "parm";
    write_text(file, "L = 16;\nBETA = 2.5\nSOLVER = \"hybridization\"\n");

    alps::Parameters legacy = alps::read_parameter_file(file.string());
    EXPECT_EQ(static_cast<int>(legacy["L"]), 16);
    EXPECT_EQ(static_cast<double>(legacy["BETA"]), 2.5);
    EXPECT_EQ(static_cast<std::string>(legacy["SOLVER"]), "hybridization");

    alps::params typed = alps::read_params_file(file.string());
    EXPECT_EQ(typed["L"].cast<int>(), 16);
    EXPECT_EQ(typed["BETA"].cast<double>(), 2.5);
    EXPECT_EQ(typed["SOLVER"].cast<std::string>(), "hybridization");
}

// HDF5 is recognized by content, so the name does not have to end in .h5.
TEST(ParameterInput, ReadsHdf5WhateverTheName) {
    alps::testing::TemporaryDirectory dir;
    for (char const * name : {"solver.in.h5", "solver.input"}) {
        auto const file = dir.path() / name;
        write_hdf5(file);

        alps::Parameters legacy = alps::read_parameter_file(file.string());
        EXPECT_EQ(static_cast<int>(legacy["L"]), 16) << name;
        EXPECT_EQ(static_cast<std::string>(legacy["SOLVER"]), "hybridization") << name;

        alps::params typed = alps::read_params_file(file.string());
        EXPECT_EQ(typed["L"].cast<int>(), 16) << name;
        EXPECT_EQ(typed["BETA"].cast<double>(), 2.5) << name;
    }
}

TEST(ParameterInput, MissingFileThrows) {
    alps::testing::TemporaryDirectory dir;
    auto const file = (dir.path() / "missing").string();
    EXPECT_THROW(alps::read_parameter_file(file), std::runtime_error);
    EXPECT_THROW(alps::read_params_file(file), std::runtime_error);
}

TEST(ParameterInput, WritesTextAndHdf5) {
    alps::testing::TemporaryDirectory dir;
    alps::Parameters parameters;
    parameters["L"] = 16;
    parameters["SOLVER"] = "hybridization";

    auto const text = dir.path() / "parm";
    auto const hdf5 = dir.path() / "parm.h5";
    alps::write_parameter_file(text.string(), parameters);
    alps::write_parameter_file(hdf5.string(), parameters);

    std::ifstream in(text);
    std::string first_line;
    std::getline(in, first_line);
    EXPECT_EQ(first_line, "L = 16;");

    for (auto const & file : {text, hdf5}) {
        alps::Parameters read = alps::read_parameter_file(file.string());
        EXPECT_EQ(static_cast<int>(read["L"]), 16) << file;
        EXPECT_EQ(static_cast<std::string>(read["SOLVER"]), "hybridization") << file;
    }
}
