// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
//
// Writes or checks archive_7146b9e1.h5 with whichever alps::hdf5::archive it
// is compiled against. The committed fixture was written by ALPSCore 7146b9e1;
// see README.md for the commands. Not part of the ALPS build.
//
//   writer write <file>   write every value in items.hpp
//   writer read <file>    read them back and report differences
#include "items.hpp"

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/boost_optional.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/map.hpp>
#include <alps/hdf5/multi_array.hpp>
#include <alps/hdf5/pair.hpp>
#include <alps/hdf5/valarray.hpp>
#include <alps/hdf5/vector.hpp>

#include <boost/multi_array.hpp>
#include <boost/optional.hpp>

#include <algorithm>
#include <exception>
#include <iostream>

namespace {
    int failures = 0;
    typedef boost::multi_array<double, 2> array_type;

    template <class T> void check(alps::hdf5::archive & ar, std::string const & path, T const & expected) {
        try {
            T value;
            ar[path] >> value;
            if (!alps_fixture::same(value, expected)) { std::cout << "DIFF " << path << "\n"; ++failures; }
        } catch (std::exception const & e) { std::cout << "FAIL " << path << ": " << e.what() << "\n"; ++failures; }
    }
}

int main(int argc, char ** argv) {
    if (argc != 3) { std::cerr << "usage: writer write|read <file>\n"; return 2; }
    std::string const mode = argv[1], file = argv[2];
    array_type array(boost::extents[2][3]);
    for (std::size_t i = 0; i != array.num_elements(); ++i)
        array.data()[i] = alps_fixture::multi_array_element(i);
    if (mode == "write") {
        alps::hdf5::archive ar(file, "w");
#define ALPS_FIXTURE_WRITE(name, type, value) { typedef ALPS_FIXTURE_ARG type T; ar["/data/" #name] << T(value); }
        ALPS_FIXTURE_ITEMS(ALPS_FIXTURE_WRITE)
        ar["/data/multi_array"] << array;
        ar["/data/optional_full"] << boost::optional<int>(17);
        ar["/data/optional_empty"] << boost::optional<int>();
        // ALPSCore's attribute spelling: no slash before the "@".
        ar["/data/int_scalar@unit"] << std::string("meV");
        ar["/data@version"] << 3;
        return 0;
    }
    alps::hdf5::archive ar(file, "r");
#define ALPS_FIXTURE_READ(name, type, value) { typedef ALPS_FIXTURE_ARG type T; check<T>(ar, "/data/" #name, T(value)); }
    ALPS_FIXTURE_ITEMS(ALPS_FIXTURE_READ)
    array_type stored;
    ar["/data/multi_array"] >> stored;
    if (!std::equal(array.data(), array.data() + array.num_elements(), stored.data())) { std::cout << "DIFF multi_array\n"; ++failures; }
    check(ar, "/data/int_scalar@unit", std::string("meV"));
    check(ar, "/data@version", 3);
    std::cout << "failures " << failures << "\n";
    return failures ? 1 : 0;
}
