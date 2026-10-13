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

#include <alps/parameter/parameterlist_p.h>
#include <iostream>
#include <cstdlib>

#include <alps/testing/stream_fixture.hpp>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>

TEST(Serialization, parameterlist_xml) {
  alps::testing::StreamFixture transcript(ALPS_TEST_SOURCE_DIR "/parameterlist_xml.input");
  alps::testing::TemporaryDirectory directory;
  {
    alps::ParameterList parameterlist;
    alps::ParameterListXMLHandler handler(parameterlist);

    alps::XMLParser parser(handler);
    parser.parse(std::cin);

    std::cout << "[Output to std::cout]\n";
    std::cout << parameterlist;

    std::cout << "\n[Output to alps::oxstream]\n";
    alps::oxstream oxs;
    oxs << parameterlist;

  } // Flush oxstream before comparing its serialized bytes.
  transcript.expect_output(ALPS_TEST_SOURCE_DIR "/parameterlist_xml.output");
}
