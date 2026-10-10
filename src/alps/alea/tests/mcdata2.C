/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Lukas Gamper <gamperl@gmail.com>,
*                            Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: nobinning.h 3520 2009-12-11 16:49:53Z gamperl $ */

#include <alps/alea/mcdata.hpp>
#include "observable_checks.hpp"
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>

TEST(AleaMcdata2, HistoricalNumericalScenarios)
{
  using namespace alps::alea;

  mcdata<double> a(0.81,0.1);
  alps_test::expect_estimate(a, {{0.81, 0.1, 5.00001e-11, 5.00001e-11}}, "a");

  mcdata<double> b = mcdata<double>(1.21,0.2);
  alps_test::expect_estimate(b, {{1.21, 0.2, 5.00001e-10, 5.00001e-11}}, "b");

  mcdata<double> c = mcdata<double>(4.5,0.3);
  alps_test::expect_estimate(c, {{4.5, 0.3, 5.00001e-10, 5.00001e-11}}, "c");

  alps_test::expect_estimate(+a, {{0.81, 0.1, 5.00001e-11, 5.00001e-11}}, "+a");
  alps_test::expect_estimate(-a, {{-0.81, 0.1, 5.00001e-11, 5.00001e-11}}, "-a");
  alps_test::expect_estimate(abs(c), {{4.5, 0.3, 5.00001e-10, 5.00001e-11}}, "abs(c)");

  mcdata<double> e = a;

  a = e;
  a += b;
  alps_test::expect_estimate(a, {{2.02, 0.2236067977, 5.00001e-10, 5.00001e-11}}, "a");

  a = e;
  a -= b;
  alps_test::expect_estimate(a, {{-0.4, 0.2236067977, 5.00001e-11, 5.00001e-11}}, "a");

  a = e;
  a *= b;
  alps_test::expect_estimate(a, {{0.9801, 0.2022003956, 5.00001e-11, 5.00001e-11}}, "a");

  a = e;
  a /= b;
  alps_test::expect_estimate(a, {{0.6694214876, 0.1381055909, 5.00001e-11, 5.00001e-11}}, "a");

  a = e;
  a += 2.;
  alps_test::expect_estimate(a, {{2.81, 0.1, 5.00001e-10, 5.00001e-11}}, "a");

  a = e;
  a -= 2.;
  alps_test::expect_estimate(a, {{-1.19, 0.1, 5.00001e-10, 5.00001e-11}}, "a");

  a = e;
  a *= 2.;
  alps_test::expect_estimate(a, {{1.62, 0.2, 5.00001e-10, 5.00001e-11}}, "a");

  a = e;
  a /= 2.;
  alps_test::expect_estimate(a, {{0.405, 0.05, 5.00001e-11, 5.00001e-12}}, "a");

  a = e;

  alps_test::expect_estimate(a + b, {{2.02, 0.2236067977, 5.00001e-10, 5.00001e-11}}, "a + b");
  alps_test::expect_estimate(a - b, {{-0.4, 0.2236067977, 5.00001e-11, 5.00001e-11}}, "a - b");
  alps_test::expect_estimate(a * b, {{0.9801, 0.2022003956, 5.00001e-11, 5.00001e-11}}, "a * b");
  alps_test::expect_estimate(a / b, {{0.6694214876, 0.1381055909, 5.00001e-11, 5.00001e-11}}, "a / b");

  alps_test::expect_estimate(a + 2., {{2.81, 0.1, 5.00001e-10, 5.00001e-11}}, "a + 2.");
  alps_test::expect_estimate(a - 2., {{-1.19, 0.1, 5.00001e-10, 5.00001e-11}}, "a - 2.");
  alps_test::expect_estimate(a * 2., {{1.62, 0.2, 5.00001e-10, 5.00001e-11}}, "a * 2.");
  alps_test::expect_estimate(a / 2., {{0.405, 0.05, 5.00001e-11, 5.00001e-12}}, "a / 2.");

  alps_test::expect_estimate(2. + a, {{2.81, 0.1, 5.00001e-10, 5.00001e-11}}, "2. + a");
  alps_test::expect_estimate(2. - a, {{1.19, 0.1, 5.00001e-10, 5.00001e-11}}, "2. - a");
  alps_test::expect_estimate(2. * a, {{1.62, 0.2, 5.00001e-10, 5.00001e-11}}, "2. * a");
  alps_test::expect_estimate(2. / a, {{2.469135802, 0.3048315806, 5.00001e-10, 5.00001e-11}}, "2. / a");

  mcdata<double> res;

  res = pow(a,2.71);

  alps_test::expect_estimate(res, {{0.5649296918, 0.1890073413, 5.00001e-11, 5.00001e-11}}, "res");

  alps_test::expect_estimate(sq(a), {{0.6561, 0.162, 5.00001e-11, 5.00001e-11}}, "sq(a)");
  alps_test::expect_estimate(cb(a), {{0.531441, 0.19683, 5.00001e-11, 5.00001e-11}}, "cb(a)");
  alps_test::expect_estimate(sqrt(a), {{0.9, 0.05555555556, 5.00001e-11, 5.00001e-12}}, "sqrt(a)");
  alps_test::expect_estimate(cbrt(a), {{0.9321697518, 0.03836089514, 5.00001e-11, 5.00001e-12}}, "cbrt(a)");
  alps_test::expect_estimate(exp(a), {{2.247907987, 0.2247907987, 5.00001e-10, 5.00001e-11}}, "exp(a)");
  alps_test::expect_estimate(log(a), {{-0.2107210313, 0.1234567901, 5.00001e-11, 5.00001e-11}}, "log(a)");

  alps_test::expect_estimate(sin(a), {{0.7242871744, 0.0689498433, 5.00001e-11, 5.00001e-12}}, "sin(a)");
  alps_test::expect_estimate(cos(a), {{0.689498433, 0.07242871744, 5.00001e-11, 5.00001e-12}}, "cos(a)");
  alps_test::expect_estimate(tan(a), {{1.050455142, 0.2103456006, 5.00001e-10, 5.00001e-11}}, "tan(a)");

  alps_test::expect_estimate(sinh(a), {{0.9015249602, 0.1346383026, 5.00001e-11, 5.00001e-11}}, "sinh(a)");
  alps_test::expect_estimate(cosh(a), {{1.346383026, 0.09015249602, 5.00001e-10, 5.00001e-12}}, "cosh(a)");
  alps_test::expect_estimate(tanh(a), {{0.6695902596, 0.05516488842, 5.00001e-11, 5.00001e-12}}, "tanh(a)");

}

TEST(AleaMcdata2, VectorNumericalScenarios)
{
  using namespace alps::alea;

  mcdata<std::vector<double> > vecA(std::vector<double>(10,0.81),std::vector<double>(10,0.1));
  alps_test::expect_uniform_estimate(vecA, {0.81, 0.1, 5.00001e-11, 5.00001e-11}, "vecA");

  mcdata<std::vector<double> > vecB(std::vector<double>(10,1.21),std::vector<double>(10,0.2));
  alps_test::expect_uniform_estimate(vecB, {1.21, 0.2, 5.00001e-10, 5.00001e-11}, "vecB");

  mcdata<std::vector<double> > vecC(std::vector<double>(10,-4.5),std::vector<double>(10,0.3));
  alps_test::expect_uniform_estimate(vecC, {-4.5, 0.3, 5.00001e-10, 5.00001e-11}, "vecC");

  std::vector<double> vec2(10,2.);

  alps_test::expect_uniform_estimate(+vecA, {0.81, 0.1, 5.00001e-11, 5.00001e-11}, "+vecA");
  alps_test::expect_uniform_estimate((-vecA), {-0.81, 0.1, 5.00001e-11, 5.00001e-11}, "(-vecA)");
  alps_test::expect_uniform_estimate(abs(vecC), {4.5, 0.3, 5.00001e-10, 5.00001e-11}, "abs(vecC)");

  mcdata<std::vector<double> > vecE = vecA;

  vecA =  vecE;
  vecA += vecB;
  alps_test::expect_uniform_estimate(vecA, {2.02, 0.2236067977, 5.00001e-10, 5.00001e-11}, "vecA");

  vecA = vecE;
  vecA -= vecB;
  alps_test::expect_uniform_estimate(vecA, {-0.4, 0.2236067977, 5.00001e-11, 5.00001e-11}, "vecA");

  vecA = vecE;
  vecA *= vecB;
  alps_test::expect_uniform_estimate(vecA, {0.9801, 0.2022003956, 5.00001e-11, 5.00001e-11}, "vecA");

  vecA = vecE;
  vecA /= vecB;
  alps_test::expect_uniform_estimate(vecA, {0.6694214876, 0.1381055909, 5.00001e-11, 5.00001e-11}, "vecA");

  vecA = vecE;
  vecA += vec2;
  alps_test::expect_uniform_estimate(vecA, {2.81, 0.1, 5.00001e-10, 5.00001e-11}, "vecA");

  vecA = vecE;
  vecA -= vec2;
  alps_test::expect_uniform_estimate(vecA, {-1.19, 0.1, 5.00001e-10, 5.00001e-11}, "vecA");

  vecA = vecE;
  vecA *= vec2;
  alps_test::expect_uniform_estimate(vecA, {1.62, 0.2, 5.00001e-10, 5.00001e-11}, "vecA");

  vecA = vecE;
  vecA /= vec2;
  alps_test::expect_uniform_estimate(vecA, {0.405, 0.05, 5.00001e-11, 5.00001e-12}, "vecA");

  vecA = vecE;

  alps_test::expect_uniform_estimate(vecA + vecB, {2.02, 0.2236067977, 5.00001e-10, 5.00001e-11}, "vecA + vecB");
  alps_test::expect_uniform_estimate(vecA - vecB, {-0.4, 0.2236067977, 5.00001e-11, 5.00001e-11}, "vecA - vecB");
  alps_test::expect_uniform_estimate(vecA * vecB, {0.9801, 0.2022003956, 5.00001e-11, 5.00001e-11}, "vecA * vecB");
  alps_test::expect_uniform_estimate(vecA / vecB, {0.6694214876, 0.1381055909, 5.00001e-11, 5.00001e-11}, "vecA / vecB");

  alps_test::expect_uniform_estimate(vecA + vec2, {2.81, 0.1, 5.00001e-10, 5.00001e-11}, "vecA + vec2");
  alps_test::expect_uniform_estimate(vecA - vec2, {-1.19, 0.1, 5.00001e-10, 5.00001e-11}, "vecA - vec2");
  alps_test::expect_uniform_estimate(vecA * vec2, {1.62, 0.2, 5.00001e-10, 5.00001e-11}, "vecA * vec2");
  alps_test::expect_uniform_estimate(vecA / vec2, {0.405, 0.05, 5.00001e-11, 5.00001e-12}, "vecA / vec2");

  alps_test::expect_uniform_estimate(vec2 + vecA, {2.81, 0.1, 5.00001e-10, 5.00001e-11}, "vec2 + vecA");
  alps_test::expect_uniform_estimate(vec2 - vecA, {1.19, 0.1, 5.00001e-10, 5.00001e-11}, "vec2 - vecA");
  alps_test::expect_uniform_estimate(vec2 * vecA, {1.62, 0.2, 5.00001e-10, 5.00001e-11}, "vec2 * vecA");
  alps_test::expect_uniform_estimate(vec2 / vecA, {2.469135802, 0.3048315806, 5.00001e-10, 5.00001e-11}, "vec2 / vecA");

  alps_test::expect_uniform_estimate(vecA + 2., {2.81, 0.1, 5.00001e-10, 5.00001e-11}, "vecA + 2.");
  alps_test::expect_uniform_estimate(vecA - 2., {-1.19, 0.1, 5.00001e-10, 5.00001e-11}, "vecA - 2.");
  alps_test::expect_uniform_estimate(vecA * 2., {1.62, 0.2, 5.00001e-10, 5.00001e-11}, "vecA * 2.");
  alps_test::expect_uniform_estimate(vecA / 2., {0.405, 0.05, 5.00001e-11, 5.00001e-12}, "vecA / 2.");

  alps_test::expect_uniform_estimate(2. + vecA, {2.81, 0.1, 5.00001e-10, 5.00001e-11}, "2. + vecA");
  alps_test::expect_uniform_estimate(2. - vecA, {1.19, 0.1, 5.00001e-10, 5.00001e-11}, "2. - vecA");
  alps_test::expect_uniform_estimate(2. * vecA, {1.62, 0.2, 5.00001e-10, 5.00001e-11}, "2. * vecA");
  alps_test::expect_uniform_estimate(2. / vecA, {2.469135802, 0.3048315806, 5.00001e-10, 5.00001e-11}, "2. / vecA");

  mcdata<std::vector<double> > resV;

  resV = pow(vecA,2.71);

  alps_test::expect_uniform_estimate(resV, {0.5649296918, 0.1890073413, 5.00001e-11, 5.00001e-11}, "resV");

  alps_test::expect_uniform_estimate(sq(vecA), {0.6561, 0.162, 5.00001e-11, 5.00001e-11}, "sq(vecA)");
  alps_test::expect_uniform_estimate(cb(vecA), {0.531441, 0.19683, 5.00001e-11, 5.00001e-11}, "cb(vecA)");
  alps_test::expect_uniform_estimate(sqrt(vecA), {0.9, 0.05555555556, 5.00001e-11, 5.00001e-12}, "sqrt(vecA)");
  alps_test::expect_uniform_estimate(cbrt(vecA), {0.9321697518, 0.03836089514, 5.00001e-11, 5.00001e-12}, "cbrt(vecA)");
  alps_test::expect_uniform_estimate(exp(vecA), {2.247907987, 0.2247907987, 5.00001e-10, 5.00001e-11}, "exp(vecA)");
  alps_test::expect_uniform_estimate(log(vecA), {-0.2107210313, 0.1234567901, 5.00001e-11, 5.00001e-11}, "log(vecA)");

  alps_test::expect_uniform_estimate(sin(vecA), {0.7242871744, 0.0689498433, 5.00001e-11, 5.00001e-12}, "sin(vecA)");
  alps_test::expect_uniform_estimate(cos(vecA), {0.689498433, 0.07242871744, 5.00001e-11, 5.00001e-12}, "cos(vecA)");
  alps_test::expect_uniform_estimate(tan(vecA), {1.050455142, 0.2103456006, 5.00001e-10, 5.00001e-11}, "tan(vecA)");

  alps_test::expect_uniform_estimate(sinh(vecA), {0.9015249602, 0.1346383026, 5.00001e-11, 5.00001e-11}, "sinh(vecA)");
  alps_test::expect_uniform_estimate(cosh(vecA), {1.346383026, 0.09015249602, 5.00001e-10, 5.00001e-12}, "cosh(vecA)");
  alps_test::expect_uniform_estimate(tanh(vecA), {0.6695902596, 0.05516488842, 5.00001e-11, 5.00001e-12}, "tanh(vecA)");

}
