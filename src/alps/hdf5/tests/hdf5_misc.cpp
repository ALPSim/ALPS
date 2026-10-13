/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include <alps/hdf5.hpp>

#include <boost/filesystem.hpp>

class my_class {
    public:
        my_class(double v = 0): d(v) {}
        void save(alps::hdf5::archive & ar) const {
            using alps::make_pvp;
            ar << make_pvp("value", d);
        }
        void load(alps::hdf5::archive & ar) {
            using alps::make_pvp;
            ar >> make_pvp("value", d);
        }
        double value() const { return d; }
    private:
        double d;
};

TEST(Hdf5, Misc) {
    alps::testing::TemporaryDirectory temporary;

    if (boost::filesystem::exists(boost::filesystem::path((temporary.path() / "data.h5").string())))
        boost::filesystem::remove(boost::filesystem::path((temporary.path() / "data.h5").string()));

    {
        alps::hdf5::archive ar((temporary.path() / "data.h5").string(), "w");
        ar << alps::make_pvp("/value", 42);
    }

    {
        alps::hdf5::archive ar((temporary.path() / "data.h5").string());
        int i;
        ar >> alps::make_pvp("/value", i);
        EXPECT_EQ(i, 42);
    }

    {
        alps::hdf5::archive ar((temporary.path() / "data.h5").string());
        std::string s;
        ar >> alps::make_pvp("/value", s);
        EXPECT_EQ(s, "42");
    }

    {
        alps::hdf5::archive ar((temporary.path() / "data.h5").string(), "w");
        std::vector<double> vec(5, 42);
        ar << alps::make_pvp("/path/2/vec", vec);
    }

    {
        std::vector<double> vec;
        // fill the vector
        alps::hdf5::archive ar((temporary.path() / "data.h5").string());
        ar >> alps::make_pvp("/path/2/vec", vec);
        EXPECT_EQ(vec, std::vector<double>(5, 42.));
    }

    {
        std::string str("foobar");
        alps::hdf5::archive ar((temporary.path() / "data.h5").string(), "w");
        ar << alps::make_pvp("/foo/bar", str);
    }

    {
        alps::hdf5::archive ar((temporary.path() / "data.h5").string());
        std::string str;
        ar >> alps::make_pvp("/foo/bar", str);
        EXPECT_EQ(str, "foobar");
    }

    {
        long *d = new long[17];
        for (long i = 0; i < 17; ++i) d[i] = i * i;
        alps::hdf5::archive ar((temporary.path() / "data.h5").string(), "w");
        ar << alps::make_pvp("/c/array", d, 17);
        delete[] d;
    }

    {
        alps::hdf5::archive ar((temporary.path() / "data.h5").string());
        std::size_t size = ar.extent("/c/array")[0];
        ASSERT_EQ(size, 17u);
        long *d = new long[size];
        ar >> alps::make_pvp("/c/array", d, size);
        for (std::size_t i = 0; i < size; ++i) EXPECT_EQ(d[i], static_cast<long>(i * i));
        delete[] d;
    }

    {
        {
                my_class c(42);
                alps::hdf5::archive ar((temporary.path() / "data.h5").string(), "w");
                ar << alps::make_pvp("/my/class", c);
        }
        {
                my_class c;
                alps::hdf5::archive ar((temporary.path() / "data.h5").string());
                ar >> alps::make_pvp("/my/class", c);
                EXPECT_EQ(c.value(), 42.);
        }
    }

    {
        alps::hdf5::archive ar((temporary.path() / "data.h5").string(), "w");
        // the parent of an attribute must exist
        ar.create_group("/foo");
        ar << alps::make_pvp("/foo/@bar", std::string("hello"));
    }

    {
        alps::hdf5::archive ar((temporary.path() / "data.h5").string());
        std::string str;
        ar >> alps::make_pvp("/foo/@bar", str);
        EXPECT_EQ(str, "hello");
    }

    boost::filesystem::remove(boost::filesystem::path((temporary.path() / "data.h5").string()));

}
