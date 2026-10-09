#include <alps/testing/stream_fixture.hpp>
/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2004 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/model.h>
#include <iostream>
#include <string>

namespace {
void check_basis_states(const char* fixture) {
  alps::testing::StreamFixture transcript(std::string(ALPS_TEST_SOURCE_DIR) + "/" + fixture + ".input");
  { // Flush serialization objects before checking the captured stream.
    alps::Parameters parms;
    std::cin >> parms;
    alps::ModelLibrary models(parms);
    alps::graph_helper<> lattices(parms);
    alps::HamiltonianDescriptor<short> ham(models.get_hamiltonian(parms["MODEL"]));
    parms.copy_undefined(ham.default_parameters());
    ham.set_parameters(parms);
    alps::basis_states_descriptor<short> basis(ham.basis(),lattices.graph());
    alps::basis_states<short> states(basis);
    std::cout << "Built states:\n" << states << std::endl;
  }
  transcript.expect_output(std::string(ALPS_TEST_SOURCE_DIR) + "/" + fixture + ".output");
}
} // namespace

TEST(ModelSerialization, Example10) { check_basis_states("example10"); }
TEST(ModelSerialization, Example12) { check_basis_states("example12"); }
