#include <alps/testing/stream_fixture.hpp>
/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/lattice.h>
#include <iostream>
#include <fstream>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

TEST(LatticeSerialization, Example9) {
  alps::testing::StreamFixture transcript(ALPS_TEST_SOURCE_DIR "/example9.input");
  { // Flush serialization objects before checking the captured stream.
    // read parameters
    alps::Parameters parameters;
    std::cin >> parameters;
    // create a graph factory with default graph type
    alps::graph_helper<> lattice(parameters);

    std::cout << "Volume of the lattice is: " << volume(lattice.lattice()) << std::endl
              << "Number of sites in a unit cell is: "
              << num_sites(lattice.unit_cell().graph()) << std::endl
              << "Number of sites is: " << num_sites(lattice.graph()) << std::endl
              << "Number of bonds is: " << num_bonds(lattice.graph()) << std::endl;

    // write the graph created from the input in XML
    std::cout << lattice.graph();
  }
  transcript.expect_output(ALPS_TEST_SOURCE_DIR "/example9.output");
}
