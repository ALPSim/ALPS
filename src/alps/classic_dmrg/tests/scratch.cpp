// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Scratch storage of the classic DMRG solver: every solve writes its block
// files into a private directory below the scratch parent and removes that
// directory when it returns or throws, leaving unrelated files alone, also
// for a relative parent. These are the guarantees the dmrg application's
// scratch-file tests checked for its former copy of the engine.
#include "support.hpp"
#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <random>
#include <string>
#include <vector>

using namespace classic_test;
namespace fs = std::filesystem;

namespace {

struct Directory {
    fs::path path;
    fs::path previous = fs::current_path();
    Directory() {
        std::random_device random;
        path = fs::temp_directory_path() / ("alps-classic-dmrg-test-" + std::to_string(random()));
        fs::create_directories(path);
    }
    ~Directory() {
        std::error_code ignored;
        fs::current_path(previous, ignored);
        fs::remove_all(path, ignored);
    }
};

std::vector<std::string> entries(fs::path const& directory) {
    std::vector<std::string> names;
    for (auto const& entry : fs::directory_iterator(directory)) names.push_back(entry.path().filename().string());
    return names;
}

// The private solver directories below `parent` that already hold files.
int populated_private_directories(fs::path const& parent) {
    int count = 0;
    for (auto const& entry : fs::directory_iterator(parent))
        if (entry.is_directory() && entry.path().filename().string().rfind("alps-classic-dmrg-", 0) == 0
            && !fs::is_empty(entry.path()))
            ++count;
    return count;
}

mdl::Hamiltonian chain() { return mdl::heisenberg(lat::chain(8, 1.0, lat::Boundary::open), 1.0); }

} // namespace

TEST(ClassicDmrgScratch, PrivateDirectoryIsUsedAndRemoved) {
    Directory parent;
    std::ofstream(parent.path / "block_ALPS_unrelated") << "keep me";
    auto options = exact_options(16, {{"Sz", 0.}});
    options.scratch_parent = parent.path;
    int seen = 0;
    auto const result = cd::solve(chain(), options, {}, [&](cd::Iteration const&) {
        seen = std::max(seen, populated_private_directories(parent.path));
    });
    EXPECT_EQ(seen, 1) << "the solve did not write into one private directory";
    EXPECT_NEAR(result.energies.at(0), ed_energies(chain(), {{"Sz", 0.}}).at(0), 1e-9);
    EXPECT_EQ(entries(parent.path), (std::vector<std::string>{"block_ALPS_unrelated"}));
    std::ifstream kept(parent.path / "block_ALPS_unrelated");
    EXPECT_EQ(std::string(std::istreambuf_iterator<char>(kept), {}), "keep me");
}

TEST(ClassicDmrgScratch, AbortedSolveRemovesItsDirectory) {
    Directory parent;
    std::ofstream(parent.path / "block_ALPS_unrelated") << "keep me";
    auto options = exact_options(16, {{"Sz", 0.}});
    options.scratch_parent = parent.path;
    struct Stop {};
    int calls = 0;
    EXPECT_THROW(cd::solve(chain(), options, {}, [&](cd::Iteration const&) {
        if (++calls == 3) throw Stop{};
    }), Stop);
    EXPECT_EQ(entries(parent.path), (std::vector<std::string>{"block_ALPS_unrelated"}));
}

TEST(ClassicDmrgScratch, RelativeParentIsResolvedOnce) {
    Directory parent;
    fs::create_directories(parent.path / "scratch");
    fs::current_path(parent.path);
    auto options = exact_options(16, {{"Sz", 0.}});
    options.scratch_parent = "scratch";
    // Changing the working directory during the solve must not redirect the
    // remaining files or the cleanup.
    bool moved = false;
    cd::solve(chain(), options, {}, [&](cd::Iteration const&) {
        if (!moved) fs::current_path(parent.path.parent_path());
        moved = true;
    });
    EXPECT_TRUE(fs::is_empty(parent.path / "scratch"));
}
