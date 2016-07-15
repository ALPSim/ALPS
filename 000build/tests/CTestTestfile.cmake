# CMake generated Testfile for 
# Source directory: /Users/hterlets/Desktop/work/DMFT/tests
# Build directory: /Users/hterlets/Desktop/work/DMFT/000build/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(cluster_helper_test "/Users/hterlets/Desktop/work/DMFT/000build/tests/cluster_helper_test")
add_test(cluster_transform_test "/Users/hterlets/Desktop/work/DMFT/000build/tests/cluster_transform_test")
add_test(double_cell_transform_test "/Users/hterlets/Desktop/work/DMFT/000build/tests/double_cell_transform_test")
add_test(normal_state_transform_test "/Users/hterlets/Desktop/work/DMFT/000build/tests/normal_state_transform_test")
subdirs(../gtest)
