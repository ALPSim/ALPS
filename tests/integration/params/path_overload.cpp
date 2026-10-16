// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// hdf5::archive is implicitly constructible from a path. A params(path) call
// must therefore either be rejected at compile time or parse the text file;
// it must never open a text file as HDF5 at run time.
#include <alps/ngs/params.hpp>

#include <boost/filesystem.hpp>

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include <fstream>
#include <stdexcept>
#include <type_traits>

#ifndef ALPS_PARAMS_LEGACY_PATH_CONSTRUCTOR
static_assert(!std::is_constructible<alps::params, boost::filesystem::path>::value);
static_assert(!std::is_constructible<alps::params, std::string>::value);
static_assert(!std::is_constructible<alps::params, const char*>::value);
#endif
static_assert(std::is_copy_constructible<alps::params>::value);
static_assert(std::is_constructible<alps::params, alps::hdf5::archive>::value);

template <typename Params>
void check_text_constructor(boost::filesystem::path const & text, std::true_type) {
#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
    Params parameters(text);
#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif
    if (parameters["L"].template cast<int>() != 10)
        throw std::runtime_error("params(path) did not parse the text grammar");
}

template <typename Params>
void check_text_constructor(boost::filesystem::path const &, std::false_type) {}

TEST(ParamsPathOverload, RejectsImplicitHdf5OrReadsText) {
    alps::testing::TemporaryDirectory temporary;
    const boost::filesystem::path text((temporary.path() / "params.parm").string());
    {
        std::ofstream output(text.string());
        output << "L=10;\n";
    }
    check_text_constructor<alps::params>(
        text, std::is_constructible<alps::params, boost::filesystem::path>());
    boost::filesystem::remove(text);
}
