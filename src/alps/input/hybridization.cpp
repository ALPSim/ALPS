// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/input/hybridization.hpp>

#include <stdexcept>

namespace alps { namespace input {

    namespace {

        // The checks of parameters that depend on each other. Passing them
        // does not guarantee that all parameters are meaningful.
        void check(hybridization_parameters const & p) {
            if (p.n_tau == 0) throw std::invalid_argument("define parameter N_TAU, the number of hybridization time slices!");
            if (p.measure_freq && p.n_matsubara == 0) throw std::invalid_argument("please specify parameter N_MATSUBARA for # of Matsubara frequencies to be measured");
            if (p.measure_legendre && p.n_legendre == 0) throw std::invalid_argument("please specify parameter N_LEGENDRE for # of Legendre coefficients to be measured");
            if (p.measure_legendre && p.n_matsubara == 0) throw std::invalid_argument("please specify parameter N_MATSUBARA for # of Matsubara frequencies");
            if (p.measure_nnt && p.n_nn == 0) throw std::invalid_argument("please specify the parameter N_nn for # of imaginary time points for the density-density correlator");
            if (p.measure_nnw && p.n_W == 0) throw std::invalid_argument("please specify the parameter N_W for # of bosonic frequencies for the density-density correlator");
            if (p.measure_g2w || p.measure_h2w) {
                if (p.n_w2 == 0) throw std::invalid_argument("please specify the parameter N_w2 for # of fermionic Matsubara frequencies for two-particle functions");
                if (p.n_W == 0) throw std::invalid_argument("please specify the parameter N_W for # of bosonic Matsubara frequencies for two-particle functions");
                if (p.n_w2 % 2 != 0) throw std::invalid_argument("parameter N_w2 must be even");
            }
            if (p.compute_vertex) {
                if (!p.measure_freq) throw std::invalid_argument("frequency measurement is required for computing the vertex, please set MEASURE_freq=1");
                if (!(p.measure_g2w || p.measure_h2w)) throw std::invalid_argument("at least one two-particle quantity is required for computing the vertex, set MEASURE_g2w=1 or MEASURE_h2w=1");
                if (p.n_matsubara + 1 < p.n_w2 / 2 + p.n_W) throw std::invalid_argument("for computing the vertex, N_MATSUBARA must be at least N_w2/2+N_W-1");
            }
        }

    }

    void define_hybridization_parameters(params & parms) {
        // simulation
        parms.define<long>("SWEEPS", "number of sweeps after thermalization");
        parms.define<long>("THERMALIZATION", "number of thermalization sweeps");
        parms.define<long>("MAX_TIME", "time limit in seconds");
        parms.define<int>("N_MEAS", "number of updates between measurements");
        parms.define<int>("OUTPUT_PERIOD", 100000, "number of sweeps between progress reports");
        parms.define<int>("N_HISTOGRAM_ORDERS", 50, "number of orders measured for the order histogram");
        parms.define<bool>("TIMESERIES", false, "use fixed binning instead of detailed binning");
        parms.define<int>("NUM_BINS", 0, "number of bins for the measurements");
        parms.define<bool>("VERBOSE", false, "print extra information");
        // updates
        parms.define<bool>("SPINFLIP", false, "perform local spin-flip updates");
        parms.define<bool>("GLOBALFLIP", false, "perform global spin-flip updates");
        // measurements
        parms.define<bool>("MEASURE_time", true, "measure the Green's function in imaginary time");
        parms.define<bool>("MEASURE_freq", false, "measure the Green's function on Matsubara frequencies");
        parms.define<bool>("MEASURE_legendre", false, "measure the Green's function in Legendre polynomials");
        parms.define<bool>("MEASURE_nn", false, "measure equal-time density-density correlations");
        parms.define<bool>("MEASURE_nnt", false, "measure the density-density correlation function in imaginary time");
        parms.define<bool>("MEASURE_nnw", false, "measure the density-density correlation function on Matsubara frequencies");
        parms.define<bool>("MEASURE_g2w", false, "measure the two-particle Green's function");
        parms.define<bool>("MEASURE_h2w", false, "measure the higher-order two-particle correlator");
        parms.define<bool>("MEASURE_sector_statistics", false, "measure sector statistics");
        parms.define<bool>("COMPUTE_VERTEX", false, "compute the vertex function");
        parms.define<bool>("ACCURATE_COVARIANCE", true, "compute the covariance of the imaginary-time Green's function accurately");
        parms.define<int>("N_MATSUBARA", 0, "number of Matsubara frequencies; required with MEASURE_freq or MEASURE_legendre");
        parms.define<int>("N_LEGENDRE", 0, "number of Legendre coefficients; required with MEASURE_legendre");
        parms.define<int>("N_nn", 0, "number of imaginary-time points for the density-density correlator; required with MEASURE_nnt");
        parms.define<int>("N_W", 0, "number of bosonic Matsubara frequencies; required with MEASURE_nnw, MEASURE_g2w or MEASURE_h2w");
        parms.define<int>("N_w2", 0, "number of fermionic Matsubara frequencies for two-particle functions; required with MEASURE_g2w or MEASURE_h2w");
        // input and output
        parms.define<std::string>("BASEPATH", "", "path in the output file under which results are stored");
        parms.define<bool>("TEXT_OUTPUT", false, "also write results as text files");
        parms.define<bool>("DELTA_IN_HDF5", false, "read the hybridization function DELTA from an HDF5 file");
        parms.define<bool>("DMFT_FRAMEWORK", false, "read HDF5 input in the layout of the DMFT framework");
        parms.define<std::string>("U_MATRIX", "", "file with the interaction matrix; replaces U, U' and J");
        parms.define<bool>("UMATRIX_IN_HDF5", false, "read U_MATRIX from an HDF5 file");
        parms.define<std::string>("MU_VECTOR", "", "file with the orbital chemical potentials; replaces MU");
        parms.define<bool>("MU_IN_HDF5", false, "read MU_VECTOR from an HDF5 file");
        parms.define<std::string>("RET_INT_K", "", "file with the retarded interaction K(tau); enables the retarded interaction");
        parms.define<bool>("K_IN_HDF5", false, "read RET_INT_K from an HDF5 file");
        // physics
        parms.define<double>("BETA", "inverse temperature");
        parms.define<int>("N_TAU", "number of imaginary-time slices of the hybridization function");
        parms.define<int>("N_ORBITALS", "number of orbitals");
        parms.define<std::string>("DELTA", "file with the hybridization function");
        parms.define<double>("MU", 0., "chemical potential; MU_VECTOR replaces it");
        // interaction: U, U' and J, or a matrix from the file U_MATRIX
        parms.define<double>("J", 0., "Hund's coupling");
        if (parms["U_MATRIX"].as<std::string>().empty()) {
            parms.define<double>("U", "interaction");
            if (parms.exists("U"))
                parms.define<double>("U'", parms["U"].as<double>() - 2 * parms["J"].as<double>(), "interaction between different orbitals; U-2J by default");
        }
    }

    hybridization_parameters read_hybridization_parameters(params const & parms) {
        hybridization_parameters p;
        p.sweeps = parms["SWEEPS"].as<std::uint64_t>();
        p.thermalization = parms["THERMALIZATION"].as<std::uint64_t>();
        p.max_time = parms["MAX_TIME"].as<long>();
        p.n_meas = parms["N_MEAS"].as<std::size_t>();
        p.output_period = parms["OUTPUT_PERIOD"].as<int>();
        p.n_histogram_orders = parms["N_HISTOGRAM_ORDERS"].as<std::size_t>();
        p.timeseries = parms["TIMESERIES"].as<bool>();
        p.num_bins = parms["NUM_BINS"].as<std::size_t>();
        p.verbose = parms["VERBOSE"].as<bool>();
        p.spinflip = parms["SPINFLIP"].as<bool>();
        p.globalflip = parms["GLOBALFLIP"].as<bool>();
        p.measure_time = parms["MEASURE_time"].as<bool>();
        p.measure_freq = parms["MEASURE_freq"].as<bool>();
        p.measure_legendre = parms["MEASURE_legendre"].as<bool>();
        p.measure_nn = parms["MEASURE_nn"].as<bool>();
        p.measure_nnt = parms["MEASURE_nnt"].as<bool>();
        p.measure_nnw = parms["MEASURE_nnw"].as<bool>();
        p.measure_g2w = parms["MEASURE_g2w"].as<bool>();
        p.measure_h2w = parms["MEASURE_h2w"].as<bool>();
        p.measure_sector_statistics = parms["MEASURE_sector_statistics"].as<bool>();
        p.compute_vertex = parms["COMPUTE_VERTEX"].as<bool>();
        p.accurate_covariance = parms["ACCURATE_COVARIANCE"].as<bool>();
        p.n_matsubara = parms["N_MATSUBARA"].as<std::size_t>();
        p.n_legendre = parms["N_LEGENDRE"].as<std::size_t>();
        p.n_nn = parms["N_nn"].as<std::size_t>();
        p.n_W = parms["N_W"].as<std::size_t>();
        p.n_w2 = parms["N_w2"].as<std::size_t>();
        p.basepath = parms["BASEPATH"].as<std::string>();
        p.text_output = parms["TEXT_OUTPUT"].as<bool>();
        p.delta_in_hdf5 = parms["DELTA_IN_HDF5"].as<bool>();
        p.dmft_framework = parms["DMFT_FRAMEWORK"].as<bool>();
        p.u_matrix = parms["U_MATRIX"].as<std::string>();
        p.umatrix_in_hdf5 = parms["UMATRIX_IN_HDF5"].as<bool>();
        p.mu_vector = parms["MU_VECTOR"].as<std::string>();
        p.mu_in_hdf5 = parms["MU_IN_HDF5"].as<bool>();
        p.ret_int_k = parms["RET_INT_K"].as<std::string>();
        p.k_in_hdf5 = parms["K_IN_HDF5"].as<bool>();
        p.beta = parms["BETA"].as<double>();
        p.n_tau = parms["N_TAU"].as<std::size_t>();
        p.n_orbitals = parms["N_ORBITALS"].as<std::size_t>();
        p.delta = parms["DELTA"].as<std::string>();
        p.mu = parms["MU"].as<double>();
        p.J = parms["J"].as<double>();
        p.U = p.u_matrix.empty() ? parms["U"].as<double>() : 0.;
        p.U_prime = p.u_matrix.empty() ? parms["U'"].as<double>() : 0.;
        p.u_ignored = !p.u_matrix.empty() && parms.supplied("U");
        p.mu_ignored = !p.mu_vector.empty() && parms.supplied("MU");
        check(p);
        return p;
    }

} }
