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
#include <alps/testing/temporary_directory.hpp>

using alps::numeric::matrix;

template<class T> class MatrixHdf5 : public ::testing::Test {};
TYPED_TEST_SUITE(MatrixHdf5, test_types);

TYPED_TEST(MatrixHdf5, hdf5)
{
    using T = TypeParam;
    alps::testing::TemporaryDirectory temporary;
    const auto filename = (temporary.path() / "alps_matrix_test.h5").string();

    {
        matrix<T> a(10,20);
        resize(a,40,40);
        fill_matrix_with_numbers(a);
        resize(a,10,20);
        matrix<T> b(a);
        b.shrink_to_fit();


        EXPECT_EQ(a.capacity() > b.capacity(), true); // maybe this should be an assert instead

        {
            alps::hdf5::archive ar(filename, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
            ar["/matrix"] << a;
        }

        EXPECT_EQ(a,b);
        matrix<T> c;
        alps::hdf5::archive ar2(filename);
        ar2["/matrix"] >> c;

        EXPECT_EQ(a,c);
        EXPECT_EQ(b,c);
    }

}

TYPED_TEST(MatrixHdf5, hdf5_matrix_matrix)
{
    using T = TypeParam;
    typedef typename matrix<matrix<T> >::col_element_iterator col_iterator;

    alps::testing::TemporaryDirectory temporary;
    const auto filename = (temporary.path() / "alps_matrix_test.h5").string();

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

        EXPECT_EQ(a.capacity() > b.capacity(), true); // maybe this should be an assert instead

        {
            alps::hdf5::archive ar(filename, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
            ar["/matrix"] << a;
        }

        EXPECT_EQ(a,b);
        matrix<matrix<T> > c;
        alps::hdf5::archive ar2(filename);
        ar2["/matrix"] >> c;

        EXPECT_EQ(a,c);
        EXPECT_EQ(b,c);
    }

}

TEST(MatrixHdf5Storage, explicit_adapter_storage_layout)
{
    alps::testing::TemporaryDirectory temporary;
    const auto filename = (temporary.path() / "alps_numeric_adapter_test.h5").string();
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
        EXPECT_TRUE(archive.extent("/dense") == std::vector<std::size_t>({3, 2}));
        std::vector<double> raw(6);
        archive.read("/dense", raw.data(), {3, 2});
        EXPECT_TRUE(raw == std::vector<double>({0.5, 10.5, 1.5, 11.5, 2.5, 12.5}));
        matrix<double> restored;
        archive["/dense"] >> restored;
        EXPECT_EQ(restored, dense);

        EXPECT_TRUE(archive.is_complex("/vector"));
        EXPECT_TRUE(archive.extent("/vector") == std::vector<std::size_t>({2, 2}));
        alps::numeric::vector<std::complex<double>> restored_vector;
        archive["/vector"] >> restored_vector;
        EXPECT_TRUE(std::equal(values.begin(), values.end(), restored_vector.begin(), restored_vector.end()));
        matrix<double> empty(2, 3);
        archive["/empty"] >> empty;
        EXPECT_EQ(empty.num_rows(), 0);
        EXPECT_EQ(empty.num_cols(), 0);
    }
}
