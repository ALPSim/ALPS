/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2004 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Axel Grzesik <axel@th.physik.uni-bonn.de>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <gtest/gtest.h>
#include <alps/model.h>
#include <alps/parser/xslt_path.h>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace model_test {
alps::ModelLibrary library() {
  std::ifstream input(alps::search_xml_library_path("models.xml"));
  if (!input) throw std::runtime_error("Cannot open the model library");
  return alps::ModelLibrary(input);
}

alps::SiteBasisDescriptor<short> descriptor(
    const std::string& name, const alps::Parameters& parameters = {}) {
  auto result = library().get_site_basis(name);
  result.set_parameters(parameters);
  return result;
}

using States = std::vector<std::vector<double>>;

template<class State>
void expect_states(const alps::site_basis<short, State>& basis,
                   const std::vector<std::string>& names, const States& expected) {
  ASSERT_EQ(basis.basis().size(), names.size());
  ASSERT_EQ(basis.size(), expected.size());
  ASSERT_FALSE(basis.empty());
  EXPECT_TRUE(basis.check_sort());
  for (std::size_t column = 0; column < names.size(); ++column)
    EXPECT_EQ(basis.basis()[column].name(), names[column]);
  for (std::size_t row = 0; row < expected.size(); ++row) {
    SCOPED_TRACE(::testing::Message() << "state " << row);
    ASSERT_EQ(expected[row].size(), names.size());
    EXPECT_EQ(basis.index(basis[row]), row);
    for (std::size_t column = 0; column < names.size(); ++column)
      EXPECT_EQ(alps::get_quantumnumber(basis[row], column),
                alps::half_integer<short>(expected[row][column])) << names[column];
  }
}

// Enumerate the historical spinful-boson contract in integer quantum numbers:
// 0 <= N <= cutoff, 0 <= J <= N*spin, and -J <= Jz <= J, in that order.
States spinful_states(int cutoff, int spin) {
  States states;
  for (int n = 0; n <= cutoff; ++n)
    for (int j = 0; j <= n * spin; ++j)
      for (int m = -j; m <= j; ++m)
        states.push_back({double(n), double(j), double(m)});
  return states;
}
} // namespace model_test

TEST(ModelBasis, FermionEnumeratesBothOccupations) {
  const alps::site_basis<short> basis(model_test::descriptor("fermion"));
  model_test::expect_states(basis, {"Nup", "Ndown"}, {{0, 0}, {0, 1}, {1, 0}, {1, 1}});
}

TEST(ModelBasis, HardcoreBosonHasEmptyAndOccupiedStates) {
  const alps::site_basis<short> basis(model_test::descriptor("hardcore boson"));
  model_test::expect_states(basis, {"N"}, {{0}, {1}});
}

namespace {
void expect_spin_basis(double spin) {
  alps::Parameters parameters;
  parameters["local_spin"] = spin;
  const alps::site_basis<short> basis(model_test::descriptor("spin", parameters));
  model_test::States expected;
  // Exactly 2*S+1 projections, in ascending order, with fixed total spin S.
  for (int i = 0; i <= int(2 * spin); ++i) expected.push_back({spin, -spin + i});
  model_test::expect_states(basis, {"S", "Sz"}, expected);
}
} // namespace

TEST(ModelBasis, SpinHalfProjections) { expect_spin_basis(0.5); }
TEST(ModelBasis, SpinOneProjections) { expect_spin_basis(1); }
TEST(ModelBasis, SpinThreeHalvesProjections) { expect_spin_basis(1.5); }
TEST(ModelBasis, SpinTwoProjections) { expect_spin_basis(2); }

TEST(ModelBasis, DefaultSpinfulBosonStates) {
  const alps::site_basis<short> basis(model_test::descriptor("spinful boson"));
  model_test::expect_states(basis, {"N", "J", "Jz"},
                            {{0, 0, 0}, {1, 0, 0}, {1, 1, -1}, {1, 1, 0}, {1, 1, 1}});
}

TEST(ModelBasis, SpinfulBosonHonorsSpinAndOccupationCutoff) {
  alps::Parameters parameters;
  parameters["boson_spin"] = 2;
  parameters["NMax"] = 2;
  const alps::site_basis<short> basis(model_test::descriptor("spinful boson", parameters));
  ASSERT_EQ(basis.size(), 35u);
  model_test::expect_states(basis, {"N", "J", "Jz"}, model_test::spinful_states(2, 2));
}

TEST(ModelBasis, TJHasOnlyVacuumAndSingleParticleStates) {
  const alps::site_basis<short> basis(model_test::descriptor("t-J"));
  model_test::expect_states(basis, {"N", "S", "Sz"}, {{0, 0, 0}, {1, 0.5, -0.5}, {1, 0.5, 0.5}});
}

TEST(ModelBasis, AlternativeTJExcludesDoubleOccupancy) {
  const alps::site_basis<short> basis(model_test::descriptor("alternative t-J"));
  model_test::expect_states(basis, {"Nup", "Ndown"}, {{0, 0}, {0, 1}, {1, 0}});
}

TEST(ModelBasis, SpinfulBosonOccupationTwoStates) {
  alps::Parameters parameters;
  parameters["NMax"] = 2;
  const alps::site_basis<short> basis(model_test::descriptor("spinful boson", parameters));
  ASSERT_EQ(basis.size(), 14u);
  model_test::expect_states(basis, {"N", "J", "Jz"}, model_test::spinful_states(2, 1));
}

TEST(ModelBasis, AddingQuantumNumberDescriptorsExtendsOccupationRange) {
  alps::Parameters parameters;
  parameters["NMax"] = 2;
  auto descriptor = model_test::descriptor("spinful boson", parameters);
  const auto other = model_test::library().get_site_basis("spinful boson");
  ASSERT_EQ(descriptor.size(), other.size());
  for (std::size_t i = 0; i < descriptor.size(); ++i) descriptor[i] += other[i];
  const alps::site_basis<short> basis(descriptor);
  // Preserve the old example7 scenario: descriptor addition extends N to four.
  // This is enumeration of the combined ranges, not a tensor-product basis.
  ASSERT_EQ(basis.size(), 55u);
  model_test::expect_states(basis, {"N", "J", "Jz"}, model_test::spinful_states(4, 1));
}

TEST(ModelBasis, SingleQuantumNumberHardcoreBoson) {
  const alps::site_basis<short, alps::single_qn_site_state<short>> basis(
      model_test::descriptor("hardcore boson"));
  ASSERT_NO_FATAL_FAILURE(model_test::expect_states(basis, {"N"}, {{0}, {1}}));
  EXPECT_FALSE(alps::is_fermionic(basis, 0));
  EXPECT_FALSE(alps::is_fermionic(basis, 1));
}

TEST(ModelBasis, SingleQuantumNumberSpinlessFermion) {
  const alps::site_basis<short, alps::single_qn_site_state<short>> basis(
      model_test::descriptor("spinless fermion"));
  ASSERT_NO_FATAL_FAILURE(model_test::expect_states(basis, {"N"}, {{0}, {1}}));
  EXPECT_FALSE(alps::is_fermionic(basis, 0));
  EXPECT_TRUE(alps::is_fermionic(basis, 1));
}
