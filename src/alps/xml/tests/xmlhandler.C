#include <alps/testing/stream_fixture.hpp>
/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2003 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/parser/xmlhandler.h>
#include <alps/parser/xmlparser.h>
#include <cstdlib>
#include <iostream>
#include <string>

TEST(XmlSerialization, Xmlhandler) {
  alps::testing::StreamFixture transcript(ALPS_TEST_SOURCE_DIR "/xmlhandler.input");
  { // Flush serialization objects before checking the captured stream.
    double v0;
    alps::SimpleXMLHandler<double> handler0("VALUE0", v0);

    double v1;
    alps::SimpleXMLHandler<double> handler1("VALUE1", v1, "value");

    alps::CompositeXMLHandler handler("TEST");
    handler.add_handler(handler0);
    handler.add_handler(handler1);

    alps::XMLParser parser(handler);

    parser.parse(std::cin);

    std::cout << v0 << std::endl
              << v1 << std::endl;
  }
  transcript.expect_output(ALPS_TEST_SOURCE_DIR "/xmlhandler.output");
}
