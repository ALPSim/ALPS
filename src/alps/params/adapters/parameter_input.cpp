// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/parameter_input.hpp>

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/make_parameters_from_file.hpp>

#include <boost/throw_exception.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <cctype>
#include <iterator>
#include <sstream>
#include <stdexcept>

namespace alps {

    namespace {

        // The HDF5 signature sits at offset 0, or at 512, 1024, 2048, ...
        // when the file has a user block. Only regular files are inspected,
        // so a pipe is never read before it is parsed.
        bool is_hdf5_file(std::filesystem::path const & filename) {
            if (!std::filesystem::is_regular_file(filename))
                return false;
            static char const signature[8] = {'\211', 'H', 'D', 'F', '\r', '\n', '\032', '\n'};
            std::ifstream in(filename, std::ios::binary);
            for (std::streamoff offset = 0; in; offset = offset ? 2 * offset : 512) {
                char buffer[8];
                in.seekg(offset);
                if (in.read(buffer, sizeof buffer) && std::equal(buffer, buffer + sizeof buffer, signature))
                    return true;
            }
            return false;
        }

        void require_readable(std::filesystem::path const & filename) {
            if (!std::filesystem::exists(filename) || std::filesystem::is_directory(filename))
                boost::throw_exception(std::runtime_error("Cannot open parameter file: " + filename.string()));
        }

        std::string read_file(std::filesystem::path const & filename) {
            std::ifstream in(filename, std::ios::binary);
            if (!in.is_open())
                boost::throw_exception(std::runtime_error("Cannot open parameter file: " + filename.string()));
            return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
        }

        Parameters read_text(std::filesystem::path const & filename) {
            std::ifstream in(filename);
            if (!in.is_open())
                boost::throw_exception(std::runtime_error("Cannot open parameter file: " + filename.string()));
            Parameters parameters;
            in >> parameters;
            return parameters;
        }

        Parameters read_xml(std::filesystem::path const & filename) {
            std::ifstream in(filename);
            if (!in.is_open())
                boost::throw_exception(std::runtime_error("Cannot open parameter file: " + filename.string()));
            Parameters parameters;
            parameters.extract_from_xml(in);
            return parameters;
        }

        bool is_xml_name(std::filesystem::path const & filename) {
            std::string extension = filename.extension().string();
            std::transform(extension.begin(), extension.end(), extension.begin(),
                           [](unsigned char c) { return std::tolower(c); });
            return extension == ".xml";
        }

        // Text without { } blocks is one run of its parameters.
        ParameterList list_from_text(std::string const & text) {
            ParameterList list;
            {
                std::istringstream in(text);
                in >> list;
            }
            if (list.empty()) {
                std::istringstream in(text);
                Parameters parameters;
                in >> parameters;
                if (!parameters.empty())
                    list.push_back(parameters);
            }
            return list;
        }

        // The text grammar reads printable characters only, inside quotes
        // that the value itself does not contain.
        void require_text_representable(Parameters const & parameters, std::filesystem::path const & filename) {
            for (Parameter const & p : parameters) {
                std::string const value = p.value().c_str();
                bool const printable = std::all_of(value.begin(), value.end(),
                    [](unsigned char c) { return std::isprint(c); });
                if (!printable || (value.find('"') != std::string::npos && value.find('\'') != std::string::npos))
                    boost::throw_exception(std::runtime_error("Parameter " + p.key()
                        + " cannot be written as parameter text to " + filename.string()));
            }
        }

    }

    Parameters read_parameter_file(std::filesystem::path const & filename) {
        require_readable(filename);
        if (!is_hdf5_file(filename))
            return is_xml_name(filename) ? read_xml(filename) : read_text(filename);
        Parameters parameters;
        hdf5::archive ar(filename, "r");
        ar["/parameters"] >> parameters;
        return parameters;
    }

    params read_params_file(std::filesystem::path const & filename) {
        require_readable(filename);
        if (!is_hdf5_file(filename))
            return make_parameters_from_file(filename);
        params parameters;
        hdf5::archive ar(filename.string(), hdf5::archive::READ);
        ar["/parameters"] >> parameters;
        return parameters;
    }

    void write_parameter_file(std::filesystem::path const & filename, Parameters const & parameters) {
        // An existing file keeps its format; a new one is HDF5 if named .h5.
        bool const hdf5 = std::filesystem::exists(filename) ? is_hdf5_file(filename)
                                                               : filename.extension() == ".h5";
        if (hdf5) {
            hdf5::archive ar(filename, "a");
            if (ar.is_group("/parameters"))
                ar.delete_group("/parameters");
            ar << make_pvp("/parameters", parameters);
        } else {
            require_text_representable(parameters, filename);
            std::ofstream out(filename);
            if (!out.is_open())
                boost::throw_exception(std::runtime_error("Cannot write parameter file: " + filename.string()));
            out << parameters;
        }
    }

    ParameterList read_parameter_list(std::filesystem::path const & filename) {
        require_readable(filename);
        if (!is_hdf5_file(filename) && !is_xml_name(filename))
            return list_from_text(read_file(filename));
        ParameterList list;
        list.push_back(read_parameter_file(filename));
        return list;
    }

    Parameters read_parameter_text(std::istream & in) {
        Parameters parameters;
        in >> parameters;
        return parameters;
    }

    ParameterList read_parameter_list_text(std::istream & in) {
        return list_from_text(std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()));
    }

}
