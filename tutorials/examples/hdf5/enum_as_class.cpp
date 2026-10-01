/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/hdf5/archive.hpp>
 
#include <boost/filesystem.hpp>

#include <string>
#include <iostream>

typedef enum { PLUS, MINUS } enum_type;

void save(
      alps::hdf5::archive & ar
    , std::string const & path
    , enum_type const & value
    , std::vector<std::size_t> size = std::vector<std::size_t>()
    , std::vector<std::size_t> chunk = std::vector<std::size_t>()
    , std::vector<std::size_t> offset = std::vector<std::size_t>()
) {
    switch (value) {
        case PLUS: ar << alps::make_pvp(path, std::string("plus")); break;
        case MINUS: ar << alps::make_pvp(path, std::string("minus")); break;
    }
}

void load(
      alps::hdf5::archive & ar
    , std::string const & path
    , enum_type & value
    , std::vector<std::size_t> chunk = std::vector<std::size_t>()
    , std::vector<std::size_t> offset = std::vector<std::size_t>()
) {
    std::string s;
    ar >> alps::make_pvp(path, s);
    value = (s == "plus" ? PLUS : MINUS);
}

int main() {

    std::string const filename = "example_enum_as_class.h5";

    if (boost::filesystem::exists(boost::filesystem::path(filename)))
        boost::filesystem::remove(boost::filesystem::path(filename));

    enum_type read, write = PLUS;

    {
        alps::hdf5::archive ar(filename, "a");
        ar << alps::make_pvp("/enum", write);
    }

    {
        alps::hdf5::archive ar(filename);
        ar >> alps::make_pvp("/enum", read);
    }

    boost::filesystem::remove(boost::filesystem::path(filename));

    bool match = (read == write);
    std::cout << (match ? "SUCCESS" : "FAILURE") << std::endl;
    return match ? EXIT_SUCCESS : EXIT_FAILURE;
}
