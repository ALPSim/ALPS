/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */
#include <alps/alea.hpp>
#include <alps/alea/util/model.hpp>
#include <alps/testing/near.hpp>

#include <iostream>
#include <random>
#include <limits>
#include <Eigen/Dense>

#include "gtest/gtest.h"

using namespace alps;

// ALPS regression coverage for the helper imported from ALPSCore 7146b9e1.
TEST(near, matrix_absolute_error)
{
    const auto near = [](const auto &a, const auto &b, double tolerance) {
        return alps::testing::NearPredFormat("a", "b", "tolerance", a, b, tolerance);
    };

    Eigen::Matrix2d a = Eigen::Matrix2d::Constant(1e10);
    Eigen::Matrix2d b = a;
    b(0, 0) += 1.0;
    EXPECT_FALSE(near(a, b, 0.25)); // Relative comparison would hide this error.
    EXPECT_TRUE(near(a, b, 1.0));  // The absolute bound is inclusive.

    a.setZero();
    b.setConstant(0.125);
    EXPECT_TRUE(near(a, b, 0.125)); // Each entry, not the matrix norm, is bounded.
    EXPECT_FALSE(near(a, b, 0.0625));

    Eigen::MatrixXd empty(0, 2), same_empty(0, 2), other_empty(2, 0);
    EXPECT_TRUE(near(empty, same_empty, 0.0));
    EXPECT_FALSE(near(empty, other_empty, 0.0));
    EXPECT_FALSE(near(a, Eigen::MatrixXd::Zero(3, 2), 1.0));

    b.setZero();
    b(0, 0) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(near(a, b, 1.0));
    b(0, 0) = std::numeric_limits<double>::infinity();
    EXPECT_FALSE(near(a, b, 1.0));

    Eigen::MatrixXcd complex_a = Eigen::MatrixXcd::Zero(1, 1);
    Eigen::MatrixXcd complex_b = complex_a;
    complex_b(0, 0) = std::complex<double>(0.3, 0.4);
    EXPECT_TRUE(near(complex_a, complex_b, 0.5));
    EXPECT_FALSE(near(complex_a, complex_b, 0.4));
}

template <typename Acc>
class eq_case
    : public ::testing::Test
{
private:
    alea::util::var1_model<double> model_;

public:
    eq_case()
    {
        Eigen::VectorXd phi0(2), veps(2);
        Eigen::MatrixXd phi1(2,2);
        phi0 << 2, 3;
        phi1 << .90, 0, 0, .30;
        veps << 1.0, 0.25;
        model_ = alea::util::var1_model<double>(phi0, phi1, veps);
    }

    void fill(Acc &acc, size_t init, size_t n)
    {
        alea::util::var1_run<double> run = model_.start();
        std::mt19937 rng(init);
        for (; run.t() < 100; run.step(rng))
            ;
        for (; run.t() < 100 + n; run.step(rng))
            acc << run.xt();
    }

    void test_empty()
    {
        Acc acc1(2), acc2;
        acc2 = acc1;
        EXPECT_EQ(acc1.result(), acc2.result());
    }

    void test_one()
    {
        Acc acc1(2), acc2(2);

        fill(acc1, 0, 1);
        EXPECT_NE(acc1.result(), acc2.result());

        acc2 = acc1;
        EXPECT_EQ(acc1.result(), acc2.result());
    }

    void test_more()
    {
        Acc acc1(2), acc2(2);
        fill(acc1, 0, 100);
        acc2 = acc1;
        EXPECT_EQ(acc1.result(), acc2.result());
    }

};

typedef ::testing::Types<
      alps::alea::mean_acc<double>
    , alps::alea::var_acc<double>
    , alps::alea::cov_acc<double>
    , alps::alea::autocorr_acc<double>
    , alps::alea::batch_acc<double>
    > has_eq;

TYPED_TEST_CASE(eq_case, has_eq);

TYPED_TEST(eq_case, test_empty) { this->test_empty(); }

TYPED_TEST(eq_case, test_one) { this->test_one(); }

TYPED_TEST(eq_case, test_more) { this->test_more(); }
