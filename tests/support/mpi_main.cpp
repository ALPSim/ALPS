// SPDX-License-Identifier: MIT
#include <gtest/gtest.h>
#include <mpi.h>
#include <string>
#include <iostream>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::GTEST_FLAG(output) = "xml:gtest-rank-" + std::to_string(rank) + ".xml";
  int local = RUN_ALL_TESTS();
  if (::testing::UnitTest::GetInstance()->test_to_run_count() == 0) {
    std::cerr << "No GoogleTest cases selected on rank " << rank << '\n';
    local = 1;
  }
  int result = 0;
  MPI_Allreduce(&local, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Finalize();
  return result;
}
