#include <alps/testing/stream_fixture.hpp>
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

/* $Id$ */

#include <alps/model.h>
#include <alps/lattice.h>
#include <alps/parameter.h>
#include <iostream>
#include <string>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

TEST(ModelSerialization, Example17) {
  alps::testing::StreamFixture transcript(ALPS_TEST_SOURCE_DIR "/example17.input");
  { // Flush serialization objects before checking the captured stream.
    alps::ParameterList parms;
    std::cin >> parms;
    for (int i=0;i<parms.size();++i) {
      alps::ModelLibrary models(parms[i]);
      alps::graph_helper<> lattice(parms[i]);
      alps::HamiltonianDescriptor<short> ham(models.get_hamiltonian(lattice,parms[i],true));
      std::cout << ham;
    }
  }
  transcript.expect_output(ALPS_TEST_SOURCE_DIR "/example17.output");
}
