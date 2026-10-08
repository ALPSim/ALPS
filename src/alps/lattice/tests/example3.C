#include <alps/testing/stream_fixture.hpp>
/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2003 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/parser/xslt_path.h>
#include <alps/lattice.h>
#include <iostream>
#include <fstream>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

TEST(LatticeSerialization, Example3) {
  alps::testing::StreamFixture transcript;
  { // Flush serialization objects before checking the captured stream.
    // create the library from an XML file
    std::ifstream in(alps::search_xml_library_path("lattices.xml"));
    alps::LatticeLibrary lib(in);

    // write one of the lattices in XML
    std::cout << lib.lattice_descriptor("square lattice 3x3");
  }
  transcript.expect_output(ALPS_TEST_SOURCE_DIR "/example3.output");
}
