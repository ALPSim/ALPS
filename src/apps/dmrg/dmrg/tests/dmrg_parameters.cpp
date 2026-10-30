// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include "dmrg.h"

#include <alps/parameter.h>
#include <gtest/gtest.h>

#include <sstream>
#include <stdexcept>
#include <string>

namespace {

typedef DMRGTask<double> task_type;

alps::params read(std::string const & text) {
    std::istringstream in(text);
    alps::Parameters parameters;
    in >> parameters;
    return task_type::read_parameters(parameters);
}

TEST(DmrgParameters, Defaults) {
    alps::params const p = read("NUMSTATES = 10");
    EXPECT_FALSE(p["COMPLEX"].as<bool>());
    EXPECT_EQ(p["NUMBER_EIGENVALUES"].as<int>(), 1);
    EXPECT_EQ(p["SWEEPS"].as<int>(), 4);
    EXPECT_EQ(p["START_SWEEP"].as<int>(), 0);
    EXPECT_EQ(p["START_DIR"].as<int>(), 0);
    EXPECT_EQ(p["START_ITER"].as<int>(), 1);
    EXPECT_EQ(p["VERBOSE"].as<int>(), 0);
    EXPECT_FALSE(p["SAVE_MPS"].as<bool>());
    EXPECT_EQ(p["NUM_WARMUP_STATES"].as<int>(), 20);
    EXPECT_EQ(p["SEED"].as<long>(), 0);
    EXPECT_EQ(p["CONSERVED_QUANTUMNUMBERS"].as<std::string>(), "");
    EXPECT_FALSE(p.exists("TRUNCATION_ERROR"));
    EXPECT_FALSE(p.exists("LANCZOS_TOLERANCE"));
    EXPECT_FALSE(p.exists("TEMP_DIRECTORY"));
    EXPECT_EQ(p["NUMSTATES"].as<int>(), 10);
}

TEST(DmrgParameters, ReadsEveryParameter) {
    alps::params const p = read(
        "COMPLEX = true\nNUMBER_EIGENVALUES = 2\nSWEEPS = 3\nSTART_SWEEP = 1\nSTART_DIR = 1\n"
        "START_ITER = 5\nVERBOSE = 1\nSAVE_MPS = 1\nNUM_WARMUP_STATES = 40\nSEED = 1234567\n"
        "CONSERVED_QUANTUMNUMBERS = \"N,Sz\"\nTRUNCATION_ERROR = 1e-7\nLANCZOS_TOLERANCE = 1e-9\n"
        "MAXSTATES = 50\nTEMP_DIRECTORY = \"/tmp/x\"");
    EXPECT_TRUE(p["COMPLEX"].as<bool>());
    EXPECT_EQ(p["NUMBER_EIGENVALUES"].as<int>(), 2);
    EXPECT_EQ(p["SWEEPS"].as<int>(), 3);
    EXPECT_EQ(p["START_SWEEP"].as<int>(), 1);
    EXPECT_EQ(p["START_DIR"].as<int>(), 1);
    EXPECT_EQ(p["START_ITER"].as<int>(), 5);
    EXPECT_EQ(p["VERBOSE"].as<int>(), 1);
    EXPECT_TRUE(p["SAVE_MPS"].as<bool>());
    EXPECT_EQ(p["NUM_WARMUP_STATES"].as<int>(), 40);
    EXPECT_EQ(p["SEED"].as<long>(), 1234567);
    EXPECT_EQ(p["CONSERVED_QUANTUMNUMBERS"].as<std::string>(), "N,Sz");
    EXPECT_EQ(p["TRUNCATION_ERROR"].as<double>(), 1e-7);
    EXPECT_EQ(p["LANCZOS_TOLERANCE"].as<double>(), 1e-9);
    EXPECT_EQ(p["MAXSTATES"].as<int>(), 50);
    EXPECT_EQ(p["TEMP_DIRECTORY"].as<std::string>(), "/tmp/x");
}

TEST(DmrgParameters, StatesBeforeMaxstatesBeforeNumstates) {
    // only the one that is read is defined, so the others need not be numbers
    alps::params const states = read("STATES = \"5,6\"\nMAXSTATES = abc\nNUMSTATES = abc");
    EXPECT_EQ(states["STATES"].as<std::string>(), "5,6");
    EXPECT_FALSE(states.exists("MAXSTATES") || states.exists("NUMSTATES"));
    alps::params const maxstates = read("MAXSTATES = 100\nNUMSTATES = abc");
    EXPECT_EQ(maxstates["MAXSTATES"].as<int>(), 100);
    EXPECT_FALSE(maxstates.exists("NUMSTATES"));
}

TEST(DmrgParameters, NonIntegerNamesParameter) {
    for (std::string const name : {"NUMBER_EIGENVALUES", "SWEEPS", "MAXSTATES", "SEED"}) {
        for (std::string const value : {"2.5", "L/2", "\"abc\""}) {
            try {
                read(name + " = " + value);
                ADD_FAILURE() << name << " = " << value << " was accepted";
            } catch (std::invalid_argument const & e) {
                EXPECT_NE(std::string(e.what()).find("parameter '" + name + "'"), std::string::npos) << e.what();
            }
        }
    }
}

}
