/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010        by Andreas Hehn <hehn@phys.ethz.ch>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/hdf5/matrix.hpp>
#include <alps/hdf5/numeric_vector.hpp>
#include "matrix_unit_tests.hpp"

using alps::numeric::matrix;

BOOST_AUTO_TEST_CASE_TEMPLATE( hdf5, T, test_types )
{
    std::string const filename = "alps_matrix_test.h5";
    if (boost::filesystem::exists(boost::filesystem::path(filename)))
        boost::filesystem::remove(boost::filesystem::path(filename));

    {
        matrix<T> a(10,20);
        resize(a,40,40);
        fill_matrix_with_numbers(a);
        resize(a,10,20);
        matrix<T> b(a);
        b.shrink_to_fit();


        BOOST_CHECK_EQUAL(a.capacity() > b.capacity(), true); // maybe this should be an assert instead

        {
            alps::hdf5::archive ar(filename, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
            ar["/matrix"] << a;
        }

        BOOST_CHECK_EQUAL(a,b);
        matrix<T> c;
        alps::hdf5::archive ar2(filename);
        ar2["/matrix"] >> c;

        BOOST_CHECK_EQUAL(a,c);
        BOOST_CHECK_EQUAL(b,c);
    }

    boost::filesystem::remove(boost::filesystem::path(filename));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( hdf5_matrix_matrix, T, test_types )
{
    typedef typename matrix<matrix<T> >::col_element_iterator col_iterator;

    std::string const filename = "alps_matrix_test.h5";
    if (boost::filesystem::exists(boost::filesystem::path(filename)))
        boost::filesystem::remove(boost::filesystem::path(filename));

    {
        matrix<matrix<T> > a(10,20);
        a.resize(40,40);
        T iota = 0;
        for(std::size_t j=0; j < num_cols(a); ++j) {
            for(std::pair<col_iterator,col_iterator> r = col(a,j); r.first != r.second; ++r.first) {
                resize(*r.first,3,5);
                for(std::size_t j=0; j < num_cols(*r.first); ++j) {
                    for(std::pair<typename matrix<T>::col_element_iterator,typename matrix<T>::col_element_iterator> rm = col(*r.first,j); rm.first != rm.second; ++rm.first) {
                        *rm.first = iota;
                        iota += 1;
                    }
                }
            }
        }
        a.resize(10,20);
        matrix<matrix<T> > b(a);
        b.shrink_to_fit();

        BOOST_CHECK_EQUAL(a.capacity() > b.capacity(), true); // maybe this should be an assert instead

        {
            alps::hdf5::archive ar(filename, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
            ar["/matrix"] << a;
        }

        BOOST_CHECK_EQUAL(a,b);
        matrix<matrix<T> > c;
        alps::hdf5::archive ar2(filename);
        ar2["/matrix"] >> c;

        BOOST_CHECK_EQUAL(a,c);
        BOOST_CHECK_EQUAL(b,c);
    }

    boost::filesystem::remove(boost::filesystem::path(filename));
}

BOOST_AUTO_TEST_CASE(explicit_adapter_storage_layout)
{
    std::string const filename = "alps_numeric_adapter_test.h5";
    matrix<double> dense(2, 3);
    dense.reserve(7, 5);
    for (std::size_t row = 0; row != 2; ++row)
        for (std::size_t column = 0; column != 3; ++column)
            dense(row, column) = 10 * row + column + 0.5;
    std::vector<std::complex<double>> const values{{1., 2.}, {-3., 4.}};
    alps::numeric::vector<std::complex<double>> vector(values.begin(), values.end());
    {
        alps::hdf5::archive archive(filename, "w");
        archive["/dense"] << dense;
        archive["/vector"] << vector;
        archive["/empty"] << matrix<double>();
    }
    {
        alps::hdf5::archive archive(filename, "r");
        // The on-disk axes stay [columns, rows], independent of reserved storage.
        BOOST_CHECK(archive.extent("/dense") == std::vector<std::size_t>({3, 2}));
        std::vector<double> raw(6);
        archive.read("/dense", raw.data(), {3, 2});
        BOOST_CHECK(raw == std::vector<double>({0.5, 10.5, 1.5, 11.5, 2.5, 12.5}));
        matrix<double> restored;
        archive["/dense"] >> restored;
        BOOST_CHECK_EQUAL(restored, dense);

        BOOST_CHECK(archive.is_complex("/vector"));
        BOOST_CHECK(archive.extent("/vector") == std::vector<std::size_t>({2, 2}));
        alps::numeric::vector<std::complex<double>> restored_vector;
        archive["/vector"] >> restored_vector;
        BOOST_CHECK(std::equal(values.begin(), values.end(), restored_vector.begin(), restored_vector.end()));
        matrix<double> empty(2, 3);
        archive["/empty"] >> empty;
        BOOST_CHECK_EQUAL(empty.num_rows(), 0);
        BOOST_CHECK_EQUAL(empty.num_cols(), 0);
    }
    boost::filesystem::remove(filename);
}
