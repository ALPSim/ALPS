// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/parameter_input.hpp>

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/make_parameters_from_file.hpp>

#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/throw_exception.hpp>

#include <algorithm>
#include <stdexcept>

namespace alps {

    namespace {

        // The HDF5 signature sits at offset 0, or at 512, 1024, 2048, ...
        // when the file has a user block.
        bool is_hdf5_file(boost::filesystem::path const & filename) {
            static char const signature[8] = {'\211', 'H', 'D', 'F', '\r', '\n', '\032', '\n'};
            boost::filesystem::ifstream in(filename, std::ios::binary);
            for (std::streamoff offset = 0; in; offset = offset ? 2 * offset : 512) {
                char buffer[8];
                in.seekg(offset);
                if (in.read(buffer, sizeof buffer) && std::equal(buffer, buffer + sizeof buffer, signature))
                    return true;
            }
            return false;
        }

        void require_readable(boost::filesystem::path const & filename) {
            if (!boost::filesystem::is_regular_file(filename))
                boost::throw_exception(std::runtime_error("Cannot open parameter file: " + filename.string()));
        }

    }

    Parameters read_parameter_file(boost::filesystem::path const & filename) {
        require_readable(filename);
        Parameters parameters;
        if (is_hdf5_file(filename)) {
            hdf5::archive ar(filename, "r");
            ar["/parameters"] >> parameters;
        } else {
            boost::filesystem::ifstream in(filename);
            if (!in.is_open())
                boost::throw_exception(std::runtime_error("Cannot open parameter file: " + filename.string()));
            in >> parameters;
        }
        return parameters;
    }

    params read_params_file(boost::filesystem::path const & filename) {
        require_readable(filename);
        if (!is_hdf5_file(filename))
            return make_parameters_from_file(filename);
        params parameters;
        hdf5::archive ar(filename.string(), hdf5::archive::READ);
        ar["/parameters"] >> parameters;
        return parameters;
    }

    void write_parameter_file(boost::filesystem::path const & filename, Parameters const & parameters) {
        if (filename.extension() == ".h5" || is_hdf5_file(filename)) {
            hdf5::archive ar(filename, "a");
            ar << make_pvp("/parameters", parameters);
        } else {
            boost::filesystem::ofstream out(filename);
            if (!out.is_open())
                boost::throw_exception(std::runtime_error("Cannot write parameter file: " + filename.string()));
            out << parameters;
        }
    }

}
