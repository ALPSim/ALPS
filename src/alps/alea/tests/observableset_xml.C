/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2006 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/testing/stream_fixture.hpp>
#include <iostream>
#include <alps/alea/observableset_p.h>

TEST(AleaXml, observableset_xml)
{
  alps::testing::StreamFixture transcript(std::string(ALPS_TEST_SOURCE_DIR) + "/observableset_xml.input");
  { // Finish XML stream destruction before comparing the transcript.

  alps::ObservableSet obs;
  alps::ObservableSetXMLHandler handler(obs);
  alps::XMLParser parser(handler);
  parser.parse(std::cin);
  
  alps::oxstream oxs;
  obs.write_xml(oxs);

  }
  transcript.expect_output(std::string(ALPS_TEST_SOURCE_DIR) + "/observableset_xml.output");
}
