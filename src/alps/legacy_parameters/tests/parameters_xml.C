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

#include <alps/parameter/parameters_p.h>
#include <alps/osiris/xdrdump.h>
#include <iostream>
#include <cstdlib>

#include <alps/testing/stream_fixture.hpp>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>

TEST(Serialization, parameters_xml) {
  alps::testing::StreamFixture transcript(ALPS_TEST_SOURCE_DIR "/parameters_xml.input");
  alps::testing::TemporaryDirectory directory;
  {
    alps::Parameters parameters;
    alps::ParametersXMLHandler handler(parameters);

    alps::XMLParser parser(handler);
    parser.parse(std::cin);

    std::cout << parameters;

    {
      alps::OXDRFileDump od(
          boost::filesystem::path((directory.path() / "parameters.dump").string()));
      od << parameters;
    }

    parameters.clear();

    {
      alps::IXDRFileDump id(
          boost::filesystem::path((directory.path() / "parameters.dump").string()));
      id >> parameters;
    }

    std::cout << parameters;

    alps::oxstream oxs;
    oxs << parameters;

  } // Flush oxstream before comparing its serialized bytes.
  transcript.expect_output(ALPS_TEST_SOURCE_DIR "/parameters_xml.output");
}
