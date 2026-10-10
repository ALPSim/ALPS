/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <gtest/gtest.h>
#include <alps/model.h>
#include <alps/lattice.h>
#include <ostream>

namespace {
struct SignScenario {
  const char* name;
  const char* lattice;
  int exchange;
  const char* diagonal_exchange;
  const char* transverse_field;
  bool has_sign_problem;
};
void PrintTo(const SignScenario& scenario, std::ostream* out) { *out << scenario.name; }
class ModelSignProblem : public ::testing::TestWithParam<SignScenario> {};

TEST_P(ModelSignProblem, DetectsHistoricalScenario) {
  const auto& scenario = GetParam();
  alps::Parameters parameters;
  parameters["MODEL_LIBRARY"] = "models.xml";
  parameters["LATTICE_LIBRARY"] = "lattices.xml";
  parameters["MODEL"] = "spin";
  parameters["LATTICE"] = scenario.lattice;
  parameters["L"] = 4;
  parameters["J"] = scenario.exchange;
  if (scenario.diagonal_exchange) parameters["J'"] = scenario.diagonal_exchange;
  if (scenario.transverse_field) parameters["Gamma"] = scenario.transverse_field;

  alps::ModelLibrary models(parameters);
  alps::graph_helper<> lattice(parameters);
  auto hamiltonian = models.get_hamiltonian(lattice, parameters);
  parameters.copy_undefined(hamiltonian.default_parameters());
  hamiltonian.set_parameters(parameters);
  EXPECT_EQ(alps::has_sign_problem(hamiltonian, lattice, parameters), scenario.has_sign_problem);
}

// Expected classifications are the eight original example8 reference results.
INSTANTIATE_TEST_SUITE_P(HistoricalModels, ModelSignProblem, ::testing::Values(
  SignScenario{"SquareAntiferromagnet", "square lattice", 1, nullptr, nullptr, false},
  SignScenario{"FrustratedSquareDefaultCoupling", "frustrated square lattice", 1, nullptr, nullptr, false},
  SignScenario{"FrustratedSquarePositiveDiagonal", "frustrated square lattice", 1, "1", nullptr, true},
  SignScenario{"FrustratedSquareNegativeDiagonal", "frustrated square lattice", 1, "-1", nullptr, false},
  SignScenario{"TriangularAntiferromagnet", "triangular lattice", 1, nullptr, nullptr, true},
  SignScenario{"ChainAntiferromagnet", "chain lattice", 1, nullptr, nullptr, false},
  SignScenario{"ChainAntiferromagnetInTransverseField", "chain lattice", 1, nullptr, "1", true},
  SignScenario{"ChainFerromagnetInTransverseField", "chain lattice", -1, nullptr, "1", false}),
  [](const ::testing::TestParamInfo<SignScenario>& info) { return info.param.name; });
} // namespace
