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

#include <alps/parapack/clone_info.h>
#include <alps/parapack/clone_info_p.h>
#include <boost/filesystem/operations.hpp>
#include <iostream>

#include <alps/testing/stream_fixture.hpp>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>

TEST(Serialization, clone_info) {
  alps::testing::StreamFixture transcript(ALPS_TEST_SOURCE_DIR "/clone_info.ip");
  alps::testing::TemporaryDirectory directory;
  {
    alps::clone_info info;
    alps::clone_info_xml_handler handler(info);

    alps::XMLParser parser(handler);
    parser.parse(std::cin);

    alps::oxstream ox(std::cout);

    ox << info;

    boost::filesystem::path xdrpath((directory.path() / "clone_info.xdr").string());
    {
      alps::OXDRFileDump dp(xdrpath);
      dp << info;
    }
    info = alps::clone_info();
    {
      alps::IXDRFileDump dp(xdrpath);
      dp >> info;
    }
    ox << info;
    boost::filesystem::remove(xdrpath);

    boost::filesystem::path h5path((directory.path() / "clone_info.h5").string());
#pragma omp critical(hdf5io)
    {
      alps::hdf5::archive ar(h5path.string(), "a");
      ar["/info"] << info;
    }
    info = alps::clone_info();
#pragma omp critical(hdf5io)
    {
      alps::hdf5::archive ar(h5path.string());
      ar["/info"] >> info;
    }
    ox << info;
    boost::filesystem::remove(h5path);

  } // Flush oxstream before comparing its serialized bytes.
  transcript.expect_output(ALPS_TEST_SOURCE_DIR "/clone_info.op");
}
