// SPDX-License-Identifier: MIT
#pragma once
#include <gtest/gtest.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace alps::testing {
// Preserve established textual serialization fixtures while tests migrate to
// GoogleTest. Numerical algorithms should assert values directly instead.
class StreamFixture {
public:
    explicit StreamFixture(const std::string& input = "")
        : old_flags_(std::cout.flags()), old_precision_(std::cout.precision()) {
        if (!input.empty()) {
            input_.open(input);
            if (!input_) throw std::runtime_error("Cannot read test input: " + input);
            old_input_ = std::cin.rdbuf(input_.rdbuf());
        }
        old_output_ = std::cout.rdbuf(output_.rdbuf());
    }
    ~StreamFixture() { restore(); }
    StreamFixture(const StreamFixture&) = delete;
    StreamFixture& operator=(const StreamFixture&) = delete;
    void expect_output(const std::string& path) {
        restore();
        std::ifstream input(path);
        ASSERT_TRUE(input.good()) << "Missing reference fixture: " << path;
        const std::string expected{std::istreambuf_iterator<char>(input), {}};
        EXPECT_EQ(output_.str(), expected) << "Serialization fixture: " << path;
    }
private:
    void restore() {
        if (old_output_) {
            std::cout.rdbuf(old_output_);
            std::cout.flags(old_flags_);
            std::cout.precision(old_precision_);
            old_output_ = nullptr;
        }
        if (old_input_) {
            std::cin.rdbuf(old_input_);
            std::cin.clear();
            old_input_ = nullptr;
        }
    }
    std::ifstream input_;
    std::ostringstream output_;
    std::streambuf* old_input_ = nullptr;
    std::streambuf* old_output_ = nullptr;
    std::ios::fmtflags old_flags_;
    std::streamsize old_precision_;
};
} // namespace alps::testing
