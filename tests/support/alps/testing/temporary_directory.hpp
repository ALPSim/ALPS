// SPDX-License-Identifier: MIT
#pragma once

#include <atomic>
#include <chrono>
#include <filesystem>
#include <random>
#include <stdexcept>
#include <string>

namespace alps::testing {
// Each instance owns a new directory, including when CTest launches cases
// from the same executable concurrently. Never changes process-wide cwd.
class TemporaryDirectory {
public:
    TemporaryDirectory() {
        static std::atomic<unsigned long> sequence{0};
        std::random_device random;
        const auto root = std::filesystem::temp_directory_path();
        for (int attempt = 0; attempt < 100; ++attempt) {
            auto candidate = root / ("alps-test-" + std::to_string(random()) + "-" +
                                     std::to_string(sequence.fetch_add(1)));
            std::error_code error;
            if (std::filesystem::create_directory(candidate, error)) {
                path_ = std::move(candidate);
                return;
            }
            if (error && error != std::errc::file_exists)
                throw std::filesystem::filesystem_error("create test directory", candidate, error);
        }
        throw std::runtime_error("Could not allocate a unique ALPS test directory");
    }
    TemporaryDirectory(const TemporaryDirectory&) = delete;
    TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;
    ~TemporaryDirectory() {
        std::error_code ignored;
        std::filesystem::remove_all(path_, ignored);
    }
    const std::filesystem::path& path() const { return path_; }
private:
    std::filesystem::path path_;
};
} // namespace alps::testing
