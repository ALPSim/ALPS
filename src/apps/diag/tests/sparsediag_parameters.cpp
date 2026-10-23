// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include "../sparsediag/sparsediag.h"

#include <alps/parameter.h>
#include <gtest/gtest.h>

#include <sstream>
#include <stdexcept>
#include <string>

namespace {

typedef SparseDiagMatrix<double> matrix_type;

alps::params read(std::string const & text) {
    std::istringstream in(text);
    alps::Parameters parameters;
    in >> parameters;
    return matrix_type::read_parameters(parameters, &matrix_type::define_parameters);
}

TEST(SparsediagParameters, Defaults) {
    alps::params const p = read("L = 4");
    EXPECT_EQ(p["CONSERVED_QUANTUMNUMBERS"].as<std::string>(), "");
    EXPECT_TRUE(p["TRANSLATION_SYMMETRY"].as<bool>());
    EXPECT_FALSE(p["PRINT_EIGENVECTORS"].as<bool>());
    EXPECT_FALSE(p["SAVE_EIGENVECTORS"].as<bool>());
    EXPECT_TRUE(p["MEASURE_ENERGY"].as<bool>());
    EXPECT_FALSE(p.exists("COMPLEX"));
    EXPECT_FALSE(p.exists("NUMBER_EIGENVALUES"));
    EXPECT_FALSE(p.exists("MAX_ITERATIONS"));
    EXPECT_FALSE(p.supplied("TOTAL_MOMENTUM"));
}

TEST(SparsediagParameters, ReadsEveryParameter) {
    alps::params const p = read(
        "CONSERVED_QUANTUMNUMBERS = \"N,Sz\"\nTRANSLATION_SYMMETRY = False\nTOTAL_MOMENTUM = \"0 0\"\n"
        "COMPLEX = 0\nPRINT_EIGENVECTORS = 1\nSAVE_EIGENVECTORS = true\nMEASURE_ENERGY = false\n"
        "NUMBER_EIGENVALUES = 3\nMAX_ITERATIONS = 200");
    EXPECT_EQ(p["CONSERVED_QUANTUMNUMBERS"].as<std::string>(), "N,Sz");
    EXPECT_FALSE(p["TRANSLATION_SYMMETRY"].as<bool>());
    EXPECT_TRUE(p.supplied("TOTAL_MOMENTUM"));
    EXPECT_FALSE(p["COMPLEX"].as<bool>());
    EXPECT_TRUE(p["PRINT_EIGENVECTORS"].as<bool>());
    EXPECT_TRUE(p["SAVE_EIGENVECTORS"].as<bool>());
    EXPECT_FALSE(p["MEASURE_ENERGY"].as<bool>());
    EXPECT_EQ(p["NUMBER_EIGENVALUES"].as<int>(), 3);
    EXPECT_EQ(p["MAX_ITERATIONS"].as<int>(), 200);
}

TEST(SparsediagParameters, NonIntegerNamesParameter) {
    for (std::string const name : {"NUMBER_EIGENVALUES", "MAX_ITERATIONS"}) {
        for (std::string const value : {"2.5", "J*3", "\"abc\""}) {
            try {
                read(name + " = " + value);
                ADD_FAILURE() << name << " = " << value << " was accepted";
            } catch (std::invalid_argument const & e) {
                EXPECT_NE(std::string(e.what()).find("parameter '" + name + "'"), std::string::npos) << e.what();
            }
        }
    }
}

TEST(SparsediagParameters, ComplexUnlessGraphOrNoMomenta) {
    EXPECT_TRUE(matrix_type::use_complex(read("L = 4"), false));
    EXPECT_FALSE(matrix_type::use_complex(read("L = 4"), true));
    EXPECT_FALSE(matrix_type::use_complex(read("TRANSLATION_SYMMETRY = false"), false));
    EXPECT_TRUE(matrix_type::use_complex(read("TRANSLATION_SYMMETRY = false\nTOTAL_MOMENTUM = 0"), false));
    EXPECT_FALSE(matrix_type::use_complex(read("COMPLEX = false"), false));
    EXPECT_TRUE(matrix_type::use_complex(read("COMPLEX = true"), true));
}

}
