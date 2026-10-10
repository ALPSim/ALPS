/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/clone_info_p.h>
#include <alps/parser/xmlparser.h>
#include <boost/filesystem/operations.hpp>
#include <iostream>

#include <alps/testing/stream_fixture.hpp>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>

TEST(Serialization, clone_phase) {
  alps::testing::StreamFixture transcript(ALPS_TEST_SOURCE_DIR "/clone_phase.ip");
  alps::testing::TemporaryDirectory directory;
  {
    alps::clone_phase phase;
    alps::clone_phase_xml_handler handler(phase);

    alps::XMLParser parser(handler);
    parser.parse(std::cin);

    alps::oxstream ox(std::cout);

    ox << phase;

    boost::filesystem::path xdrpath((directory.path() / "clone_phase.xdr").string());
    {
      alps::OXDRFileDump dp(xdrpath);
      dp << phase;
    }
    phase = alps::clone_phase();
    {
      alps::IXDRFileDump dp(xdrpath);
      dp >> phase;
    }
    ox << phase;
    boost::filesystem::remove(xdrpath);

    boost::filesystem::path h5path((directory.path() / "clone_phase.h5").string());
#pragma omp critical(hdf5io)
    {
      alps::hdf5::archive ar(h5path.string(), "a");
      ar["/phase"] << phase;
    }
    phase = alps::clone_phase();
#pragma omp critical(hdf5io)
    {
      alps::hdf5::archive ar(h5path.string());
      ar["/phase"] >> phase;
    }
    ox << phase;
    boost::filesystem::remove(h5path);

  } // Flush oxstream before comparing its serialized bytes.
  transcript.expect_output(ALPS_TEST_SOURCE_DIR "/clone_phase.op");
}
