"""Exercise discovery and failure propagation with a real, tiny C++ suite."""
import json
import os
from pathlib import Path
import subprocess

import pytest

SOURCE = Path(__file__).resolve().parents[2]
pytestmark = pytest.mark.skipif(
    not os.environ.get("ALPS_DIR"), reason="requires configured C++ dependencies")


def run(*args):
    return subprocess.run(args, capture_output=True, text=True, timeout=180)


@pytest.mark.parametrize("generator", ["Ninja", "Ninja Multi-Config"])
def test_discovery_isolation_and_failure_propagation(tmp_path, generator):
    (tmp_path / "CMakeLists.txt").write_text(f'''
cmake_minimum_required(VERSION 3.27)
project(testing_contract CXX)
set(CMAKE_CXX_STANDARD 17)
set(ALPS_BUILD_TESTING ON)
enable_testing()
include("{SOURCE.as_posix()}/cmake/ALPSTesting.cmake")
add_executable(probe probe.cpp)
alps_add_gtest(probe COMPONENT harness LABELS compatibility TIMEOUT 10)
''')
    (tmp_path / "probe.cpp").write_text('''
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include <filesystem>
#include <fstream>
TEST(Harness, IsolatedFiles) {
    std::filesystem::path removed;
    {
        alps::testing::TemporaryDirectory first, second;
        EXPECT_NE(first.path(), second.path());
        removed = first.path();
        std::ofstream(first.path() / "shared-name") << "one";
        EXPECT_FALSE(std::filesystem::exists(second.path() / "shared-name"));
    }
    EXPECT_FALSE(std::filesystem::exists(removed));
}
TEST(Harness, DeliberateFailure) { EXPECT_EQ(1, 2); }
''')
    build = tmp_path / "build"
    configured = run("cmake", "-S", str(tmp_path), "-B", str(build), "-G", generator,
                     *json.loads(os.environ.get("ALPS_TEST_CMAKE_ARGS", "[]")))
    assert configured.returncode == 0, configured.stdout + configured.stderr
    compiled = run("cmake", "--build", str(build), "--config", "Release", "-j", "2")
    assert compiled.returncode == 0, compiled.stdout + compiled.stderr
    listing = run("ctest", "--test-dir", str(build), "-C", "Release", "--show-only=json-v1")
    assert listing.returncode == 0, listing.stderr
    cases = json.loads(listing.stdout)["tests"]
    assert {case["name"] for case in cases} == {
        "harness.probe.Harness.IsolatedFiles", "harness.probe.Harness.DeliberateFailure"}
    for case in cases:
        props = {prop["name"]: prop["value"] for prop in case["properties"]}
        assert set(props["LABELS"]) == {"harness", "compatibility"}
        assert props["TIMEOUT"] == 10
        assert "OMP_NUM_THREADS=set:1" in props["ENVIRONMENT_MODIFICATION"]
        assert "OPENBLAS_NUM_THREADS=set:1" in props["ENVIRONMENT_MODIFICATION"]
    passed = run("ctest", "--test-dir", str(build), "-C", "Release",
                 "-R", "IsolatedFiles", "--no-tests=error", "--output-on-failure")
    assert passed.returncode == 0, passed.stdout + passed.stderr
    failed = run("ctest", "--test-dir", str(build), "-C", "Release",
                 "-R", "DeliberateFailure", "--no-tests=error", "--output-on-failure")
    assert failed.returncode != 0
    assert "DeliberateFailure" in failed.stdout
    # Losing every case in one executable must fail even when the repository
    # still has other registered tests. Whole-suite --no-tests cannot catch it.
    (tmp_path / "probe.cpp").write_text("#include <gtest/gtest.h>\n")
    compiled = run("cmake", "--build", str(build), "--config", "Release", "-j", "2")
    assert compiled.returncode == 0, compiled.stdout + compiled.stderr
    empty = run("ctest", "--test-dir", str(build), "-C", "Release", "--show-only=json-v1")
    assert empty.returncode != 0
    assert "zero GoogleTest cases" in empty.stdout + empty.stderr


@pytest.mark.parametrize("testing", ["OFF", "ON"])
def test_unused_googletest_is_not_found_or_fetched(tmp_path, testing):
    (tmp_path / "CMakeLists.txt").write_text(f'''
cmake_minimum_required(VERSION 3.27)
project(no_testing NONE)
set(ALPS_BUILD_TESTING {testing})
set(FETCHCONTENT_FULLY_DISCONNECTED ON)
include("{SOURCE.as_posix()}/cmake/ALPSTesting.cmake")
if(TARGET GTest::gtest OR DEFINED googletest_SOURCE_DIR)
  message(FATAL_ERROR "Unused test dependency was found or fetched")
endif()
''')
    configured = run("cmake", "-S", str(tmp_path), "-B", str(tmp_path / "build"))
    assert configured.returncode == 0, configured.stdout + configured.stderr
