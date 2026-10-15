// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/invocation.hpp. Adaptation: the
// expression-stream boundary (alpsx's disorder context) is not ported; the
// scratch parent defaults to std::filesystem::temp_directory_path().
//
// This private header is included only by classic.cpp. The dmtk engine keeps
// process-wide state (its quantum-number registry, scratch file list, work
// buffers, random seeds and output streams); one invocation claims all of it
// under a mutex, installs a fresh-process state and restores the previous
// state afterwards. No dmtk object may escape the invocation's lifetime.
#pragma once

#include "dmtk/dmtk.hpp"

#include <filesystem>
#include <mutex>
#include <random>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <system_error>
#include <utility>

namespace dmtk {
namespace invocation_detail {

inline thread_local bool active = false;

class Claim {
public:
    Claim() {
        if (active) throw std::logic_error("classic DMRG cannot be invoked recursively");
        active = true;
    }
    ~Claim() noexcept { active = false; }
    Claim(Claim const&) = delete;
    Claim& operator=(Claim const&) = delete;
};

inline std::mutex& mutex() {
    static std::mutex value;
    return value;
}

class ScratchDirectory {
public:
    explicit ScratchDirectory(std::filesystem::path parent) {
        if (parent.empty()) parent = std::filesystem::temp_directory_path();
        parent = std::filesystem::absolute(parent);
        if (!std::filesystem::is_directory(parent))
            throw std::invalid_argument("classic DMRG scratch parent is not a directory");
        std::random_device random;
        for (int attempt = 0; attempt < 256; ++attempt) {
            std::ostringstream name;
            name << "alps-classic-dmrg-" << std::hex << random() << '-' << random();
            auto candidate = parent / name.str();
            std::error_code error;
            if (std::filesystem::create_directory(candidate, error)) {
                path_.swap(candidate); // No allocation after the directory is created.
                return;
            }
            if (error && error != std::errc::file_exists)
                throw std::filesystem::filesystem_error("cannot create classic DMRG scratch directory",
                                                        candidate, error);
        }
        throw std::runtime_error("cannot create unique classic DMRG scratch directory");
    }
    ~ScratchDirectory() noexcept {
        try {
            std::error_code ignored;
            std::filesystem::remove_all(path_, ignored);
        } catch (...) {
            // Cleanup must never replace a solver exception. Only this
            // invocation's uniquely created directory is removed.
        }
    }
    std::filesystem::path const& path() const noexcept { return path_; }
    ScratchDirectory(ScratchDirectory const&) = delete;
    ScratchDirectory& operator=(ScratchDirectory const&) = delete;

private:
    std::filesystem::path path_;
};

class NullBuffer : public std::streambuf {
    int_type overflow(int_type value) override { return traits_type::not_eof(value); }
    std::streamsize xsputn(char const*, std::streamsize count) override { return count; }
};

// Routes the engine's diagnostics of the calling thread to the caller's log
// stream (or discards them). A separate ostream owns all formatting state, so
// the caller's stream flags are left untouched.
class OutputScope {
public:
    explicit OutputScope(std::ostream* log)
        : stream_(log && log->rdbuf() ? log->rdbuf() : &discard_),
          previous_output_(output_stream), previous_error_(error_stream) {
        output_stream = &stream_;
        error_stream = &stream_;
    }
    ~OutputScope() noexcept {
        output_stream = previous_output_;
        error_stream = previous_error_;
    }
    OutputScope(OutputScope const&) = delete;
    OutputScope& operator=(OutputScope const&) = delete;

private:
    NullBuffer discard_;
    std::ostream stream_;
    std::ostream* previous_output_;
    std::ostream* previous_error_;
};

} // namespace invocation_detail

class InvocationGuard {
public:
    explicit InvocationGuard(std::filesystem::path scratch_parent = {}, std::ostream* log = nullptr)
        : lock_(invocation_detail::mutex()), scratch_(std::move(scratch_parent)), output_(log),
          files_(scratch_.path().string().c_str()),
          real_random_(Vector<double>::random_state()),
          complex_random_(Vector<std::complex<double>>::random_state()) {
        // All potentially throwing setup is complete. These swaps install a
        // fresh-process state and retain the previous state for restoration.
        QN::swap_registry(registry_);
        tmp_files.swap(files_);
        globals_double.swap(real_buffers_);
        globals_complex.swap(complex_buffers_);
        Vector<double>::random_state() = 42;
        Vector<std::complex<double>>::random_state() = 42;
    }
    ~InvocationGuard() noexcept {
        // The caller must have destroyed every engine object before this point.
        Vector<double>::random_state() = real_random_;
        Vector<std::complex<double>>::random_state() = complex_random_;
        globals_complex.swap(complex_buffers_);
        globals_double.swap(real_buffers_);
        tmp_files.swap(files_);
        QN::swap_registry(registry_);
    }
    std::filesystem::path const& scratch_directory() const noexcept { return scratch_.path(); }
    InvocationGuard(InvocationGuard const&) = delete;
    InvocationGuard& operator=(InvocationGuard const&) = delete;

private:
    // Reverse destruction restores streams, removes scratch, then unlocks.
    // Claim is first so recursive calls fail before acquiring the mutex.
    invocation_detail::Claim claim_;
    std::unique_lock<std::mutex> lock_;
    invocation_detail::ScratchDirectory scratch_;
    invocation_detail::OutputScope output_;
    FileList files_;
    QN::Registry registry_;
    DMTKglobals<double> real_buffers_;
    DMTKglobals<std::complex<double>> complex_buffers_;
    long real_random_;
    long complex_random_;
};

} // namespace dmtk
