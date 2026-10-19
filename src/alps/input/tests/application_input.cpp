// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/input/dmft.hpp>
#include <alps/input/hirschfye.hpp>
#include <alps/parameter_input.hpp>
#include <gtest/gtest.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

alps::Parameters text(std::string const & parameters) {
    std::istringstream in(parameters);
    return alps::read_parameter_text(in);
}

TEST(DmftInput, RunDefaults) {
    alps::input::dmft_run_parameters const p = alps::input::read_dmft_run_parameters(text("SOLVER = \"hybridization\""));
    EXPECT_FALSE(p.cluster_loop);
    EXPECT_FALSE(p.omega_loop);
    EXPECT_EQ(p.solver, "hybridization");
    EXPECT_FALSE(p.sc_write_delta);
    EXPECT_FALSE(p.antiferromagnet);
    EXPECT_FALSE(p.symmetrization);
    EXPECT_FALSE(p.semicircle_hilbert);
    EXPECT_EQ(p.sizes.flavors, 2);
    EXPECT_EQ(p.sizes.sites, 1);
}

TEST(DmftInput, RunNeedsSolverUnlessClusterLoop) {
    EXPECT_THROW(alps::input::read_dmft_run_parameters(text("OMEGA_LOOP = 1")), std::runtime_error);
    EXPECT_TRUE(alps::input::read_dmft_run_parameters(text("CLUSTER_LOOP = 1")).cluster_loop);
}

TEST(DmftInput, LoopReadsWhatEachLoopNeeds) {
    std::string const common = "N = 16\nBETA = 10\nCONVERGED = 0.001\nSYMMETRIZATION = 1\nBASENAME = \"run\"\n";
    alps::input::dmft_loop_parameters const dressed =
        alps::input::read_dmft_loop_parameters(text(common), alps::input::dmft_loop::dressed_imaginary_time);
    EXPECT_EQ(dressed.n_tau, 16);
    EXPECT_EQ(dressed.beta, 10.);
    EXPECT_EQ(dressed.converged, 0.001);
    EXPECT_TRUE(dressed.symmetrization);
    EXPECT_EQ(dressed.max_it, 1000);
    EXPECT_EQ(dressed.h, 0.);
    EXPECT_EQ(dressed.basename, "run");
    EXPECT_EQ(dressed.g0_tau_output, "G0tau_output");
    EXPECT_FALSE(dressed.h_init);
    // MU is needed by the other loops, NMATSUBARA by the Matsubara loop.
    EXPECT_THROW(alps::input::read_dmft_loop_parameters(text(common), alps::input::dmft_loop::imaginary_time), std::runtime_error);
    EXPECT_THROW(alps::input::read_dmft_loop_parameters(text(common + "MU = 0.5"), alps::input::dmft_loop::matsubara), std::runtime_error);
    alps::input::dmft_loop_parameters const omega = alps::input::read_dmft_loop_parameters(
        text(common + "MU = 0.5\nNMATSUBARA = 64\nSITES = 1\nRELAX_RATE = 0.5\nH = 0.1\nH_INIT = 0.2\nG0TAU_input = \"g0\""),
        alps::input::dmft_loop::matsubara);
    EXPECT_EQ(omega.n_matsubara, 64);
    EXPECT_EQ(omega.mu, 0.5);
    EXPECT_EQ(omega.relax_rate, 0.5);
    EXPECT_EQ(omega.h, 0.1);
    EXPECT_EQ(omega.h_init, std::string("0.2"));
    EXPECT_EQ(omega.initial_g0_tau_output, std::string("g0"));
    EXPECT_EQ(omega.g0_omega_output, "G0omega_output");
}

TEST(DmftInput, InitialG0IgnoresEmptyInputFile) {
    alps::input::dmft_initial_g0_parameters const p =
        alps::input::read_dmft_initial_g0_tau(text("N = 8\nG0TAU_INPUT = \"\"\nG0TAU_input = \"out\""));
    EXPECT_EQ(p.n_tau, 8);
    EXPECT_FALSE(p.input);
    EXPECT_EQ(p.output, std::string("out"));
    alps::input::dmft_initial_g0_parameters const q = alps::input::read_dmft_initial_g0_omega(
        text("N = 8\nNMATSUBARA = 32\nBETA = 5\nMU = 1\nINSULATING = 1\nG0OMEGA_INPUT = \"g0\""));
    EXPECT_EQ(q.n_matsubara, 32);
    EXPECT_TRUE(q.insulating);
    EXPECT_EQ(q.input, std::string("g0"));
}

TEST(DmftInput, Bandstructures) {
    EXPECT_EQ(alps::input::read_dmft_bandstructure(text("t = 1")).kind, alps::input::dmft_bandstructure::semicircle);
    EXPECT_EQ(alps::input::read_dmft_bandstructure(text("TWODBS = square")).kind, alps::input::dmft_bandstructure::square);
    EXPECT_EQ(alps::input::read_dmft_bandstructure(text("TWODBS = hexagonal")).kind, alps::input::dmft_bandstructure::hexagonal);
    alps::input::dmft_bandstructure_choice const dos = alps::input::read_dmft_bandstructure(text("DOSFILE = \"dos.dat\"\nTWODBS = square"));
    EXPECT_EQ(dos.kind, alps::input::dmft_bandstructure::dos);
    EXPECT_EQ(dos.dosfile, "dos.dat");
    EXPECT_EQ(alps::input::read_semicircle_hoppings(text("FLAVORS = 4\nt = 1\nt1 = 0.5")), (std::vector<double>{1., 0.5}));
    EXPECT_THROW(alps::input::read_semicircle_hoppings(text("FLAVORS = 3\nt = 1")), std::logic_error);
    alps::input::dmft_lattice_hoppings const h = alps::input::read_lattice_hoppings(text("tprime = 0.2"));
    EXPECT_EQ(h.t, 1.);
    EXPECT_EQ(h.tprime, 0.2);
    EXPECT_EQ(h.L, 128);
}

TEST(DmftInput, FourierDefaultsTheMoments) {
    alps::input::dmft_fourier_parameters const p =
        alps::input::read_dmft_fourier_parameters(text("BETA = 2\nMU = 0.5\nU = 3\nEPS_1 = 0.25"), true, true);
    EXPECT_EQ(p.beta, 2.);
    EXPECT_EQ(p.mu, 0.5);
    EXPECT_EQ(p.U, 3.);
    EXPECT_EQ(p.eps, (std::vector<double>{0., 0.25}));
    EXPECT_EQ(p.epssq, (std::vector<double>{1., 1.}));
    alps::input::dmft_band_moments const m = alps::input::read_band_moments(text("EPSSQ_0 = 2"), 2);
    EXPECT_FALSE(m.eps[0]);
    EXPECT_EQ(m.epssq[0], 2.);
    EXPECT_FALSE(m.epssq[1]);
}

TEST(DmftInput, InteractionMatrix) {
    alps::input::dmft_interaction_parameters const two = alps::input::read_dmft_interaction(text("U = 4"));
    EXPECT_EQ(two.U, 4.);
    EXPECT_FALSE(two.u_matrix);
    alps::input::dmft_interaction_parameters const four = alps::input::read_dmft_interaction(text("FLAVORS = 4\nU = 4\nJ = 0.5"));
    EXPECT_EQ(four.Uprime, 3.);
    EXPECT_EQ(alps::input::read_dmft_interaction(text("FLAVORS = 4\nU = 4\nJ = 0.5\nU' = 2")).Uprime, 2.);
    EXPECT_THROW(alps::input::read_dmft_interaction(text("FLAVORS = 4\nU = 4")), std::runtime_error);
    EXPECT_EQ(alps::input::read_dmft_interaction(text("U_MATRIX = \"u.dat\"")).u_matrix, std::string("u.dat"));
}

TEST(DmftInput, ExternalSolverReadsHybridizationOnlyWhenWritten) {
    alps::input::dmft_external_solver_parameters const plain =
        alps::input::read_dmft_external_solver(text("N = 8\nTMPNAME = \"tmp\""), alps::input::dmft_loop::imaginary_time);
    EXPECT_EQ(plain.n_tau, 8);
    EXPECT_EQ(plain.tmpname, std::string("tmp"));
    EXPECT_FALSE(plain.sc_write_delta);
    alps::input::dmft_external_solver_parameters const delta = alps::input::read_dmft_external_solver(
        text("N = 8\nNMATSUBARA = 16\nSC_WRITE_DELTA = 1\nBETA = 4\nMU = 1\nH = 0.1"), alps::input::dmft_loop::matsubara);
    EXPECT_EQ(delta.n_matsubara, 16);
    EXPECT_EQ(delta.beta, 4.);
    EXPECT_EQ(delta.mu, 1.);
    EXPECT_EQ(delta.h, 0.1);
    EXPECT_EQ(alps::input::read_dmft_max_time(text("MAX_TIME = 30"), 60), 30u);
    EXPECT_EQ(alps::input::read_dmft_max_time(text("N = 1"), 60), 60u);
}

TEST(HirschFyeInput, ReadsSolverParameters) {
    std::string const required = "THERMALIZATION = 10\nSWEEPS = 100\nBETA = 8\nN = 16\nU = 2\nTOLERANCE = 1e-6\nNMATSUBARA = 32\n";
    alps::input::hirschfye_parameters const p = alps::input::read_hirschfye_parameters(text(required));
    EXPECT_EQ(p.thermalization, 10);
    EXPECT_EQ(p.sweeps, 100);
    EXPECT_EQ(p.beta, 8.);
    EXPECT_EQ(p.n_tau, 16);
    EXPECT_EQ(p.sites, 1);
    EXPECT_EQ(p.U, 2.);
    EXPECT_EQ(p.fourpoint_interval, 1000);
    EXPECT_FALSE(p.measure_fourpoint);
    EXPECT_EQ(p.tolerance, 1e-6);
    EXPECT_EQ(p.max_time, 86400);
    EXPECT_EQ(p.n_matsubara, 32);
    EXPECT_THROW(alps::input::read_hirschfye_parameters(text("BETA = 8")), std::runtime_error);
}

TEST(HirschFyeInput, ResultNeedsSites) {
    EXPECT_THROW(alps::input::read_hirschfye_result_parameters(text("N = 16\nNMATSUBARA = 32")), std::runtime_error);
    alps::input::hirschfye_result_parameters const p = alps::input::read_hirschfye_result_parameters(
        text("N = 16\nNMATSUBARA = 32\nSITES = 1\nOUTFILE = \"out.h5\""));
    EXPECT_EQ(p.sites, 1);
    EXPECT_FALSE(p.checkpoint);
    EXPECT_EQ(p.outfile, std::string("out.h5"));
    alps::input::hirschfye_input_sizes const sizes = alps::input::read_hirschfye_input_sizes(text("NMATSUBARA = 32"));
    EXPECT_EQ(sizes.n_matsubara, 32);
    EXPECT_EQ(sizes.sites, 1);
    EXPECT_EQ(sizes.flavors, 2);
}

}
