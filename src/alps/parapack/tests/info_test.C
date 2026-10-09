/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/clone_info.h>
#include <gtest/gtest.h>

TEST(CloneInfo, TracksPhasesSeedAndProgress) {
  alps::Parameters params;
  params["SEED"] = 29832;
  alps::clone_info info(0, params, "info_test");
  EXPECT_EQ(info.clone_id(), 0);
  EXPECT_TRUE(info.has_seed());
  alps::clone_info same_clone(0, params, "same"), other_clone(1, params, "other");
  EXPECT_EQ(info.worker_seed(), same_clone.worker_seed());
  EXPECT_NE(info.worker_seed(), other_clone.worker_seed());
  EXPECT_DOUBLE_EQ(info.progress(), 0);
  info.start("test 1");
  EXPECT_EQ(info.phase(), "test 1");
  info.stop();
  info.start("test 2");
  EXPECT_EQ(info.phase(), "test 2");
  info.stop();
  ASSERT_EQ(info.phases().size(), 2);
  EXPECT_GE(info.elapsed().total_microseconds(), 0);
  info.set_progress(.593483);
  EXPECT_DOUBLE_EQ(info.progress(), .593483);
  info.set_progress(1);
  EXPECT_DOUBLE_EQ(info.progress(), 1);
  EXPECT_EQ(info.dumpfile(), "info_test.clone1");
  EXPECT_EQ(info.dumpfile_h5(), "info_test.clone1.h5");
  EXPECT_EQ(info.dumpfile_xdr(), "info_test.clone1.xdr");
}
