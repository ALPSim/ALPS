// SPDX-License-Identifier: MIT
#include <gtest/gtest.h>
#include <alps/ngs/mcoptions.hpp>
#include <alps/parseargs.hpp>
#include <alps/cli/options.hpp>

#include <initializer_list>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

template<class Options>
Options parse(std::initializer_list<char const*> args) {
    std::vector<std::string> storage{"cli-contract"};
    for (auto const* arg : args) storage.emplace_back(arg);
    std::vector<char*> argv;
    for (auto& arg : storage) argv.push_back(&arg[0]);
    argv.push_back(nullptr);
    return Options(static_cast<int>(storage.size()), argv.data());
}

template<class Exception>
void rejects(std::initializer_list<char const*> args, char const* message) {
    bool rejected = false;
    try { parse<alps::mcoptions>(args); }
    catch (Exception const&) { rejected = true; }
    ASSERT_TRUE((rejected)) << message;
}

struct capture_stdout {
    std::ostringstream text;
    std::streambuf* previous = std::cout.rdbuf(text.rdbuf());
    ~capture_stdout() { std::cout.rdbuf(previous); }
};

TEST(CommandLine, MonteCarloOptions) {
    auto defaults = parse<alps::mcoptions>({"sample.in.h5"});
    ASSERT_TRUE((defaults.valid && !defaults.resume && defaults.time_limit == 0)) << "mcoptions default validity, resume and time limit";
    ASSERT_TRUE((defaults.type == alps::mcoptions::SINGLE && defaults.checkpoint_file.empty())) << "mcoptions default execution and checkpoint";
    ASSERT_TRUE((defaults.input_file == "sample.in.h5" && defaults.output_file == "sample.out.h5")) << "mcoptions default input/output names";
    ASSERT_TRUE((parse<alps::mcoptions>({"sample.out.h5"}).output_file == "sample.out.h5")) << "Existing output checkpoint keeps its name";
    ASSERT_TRUE((parse<alps::mcoptions>({"sample.h5"}).output_file == "sample.out.h5")) << "General input extension is replaced";
    ASSERT_TRUE((parse<alps::mcoptions>({"sample"}).output_file == "sample.out.h5")) << "Extensionless input gets an output suffix";
    // These are existing substring rules, not a request to normalize filenames.
    ASSERT_TRUE((parse<alps::mcoptions>({"sample.out.h5.backup"}).output_file == "sample.out.h5.backup")) << "Existing output-name substring behavior";

    auto positional = parse<alps::mcoptions>({"input.h5", "output.h5", "checkpoint.h5", "-c", "-T", "120"});
    ASSERT_TRUE((positional.input_file == "input.h5" && positional.output_file == "output.h5"
            && positional.checkpoint_file == "checkpoint.h5" && positional.resume
            && positional.time_limit == 120)) << "mcoptions positional filenames and short options";
    auto named = parse<alps::mcoptions>({"--input-file", "input.h5", "--output-file", "output.h5",
            "--checkpoint-file", "checkpoint.h5", "--continue", "--time-limit", "17", "--single"});
    ASSERT_TRUE((named.input_file == "input.h5" && named.output_file == "output.h5"
            && named.checkpoint_file == "checkpoint.h5" && named.resume && named.time_limit == 17
            && named.type == alps::mcoptions::SINGLE)) << "mcoptions named options";

    {
        capture_stdout output;
        auto help = parse<alps::mcoptions>({"--help"});
        ASSERT_TRUE((!help.valid && help.input_file.empty())) << "Help does not require an input file";
        ASSERT_TRUE((output.text.str().find("--input-file") != std::string::npos
                && output.text.str().find("--time-limit") != std::string::npos)) << "Help describes the current grammar";
    }
    rejects<std::invalid_argument>({}, "Missing input is rejected");
    rejects<alps::cli::error>({"sample.h5", "--unknown"},
            "Unknown mcoptions option is rejected");
    rejects<alps::cli::error>({"sample.h5", "-T", "not-a-number"},
            "Nonnumeric time limit is rejected");
    rejects<alps::cli::error>({"a", "b", "c", "d"},
            "Extra positional filename is rejected");

#ifndef ALPS_NGS_SINGLE_THREAD
    ASSERT_TRUE((parse<alps::mcoptions>({"sample.h5", "--threaded"}).type == alps::mcoptions::THREADED)) << "Threaded execution option";
#else
    rejects<alps::cli::error>({"sample.h5", "--threaded"},
            "Single-thread build rejects threaded execution");
#endif
#ifdef ALPS_HAVE_MPI
    ASSERT_TRUE((parse<alps::mcoptions>({"sample.h5", "--mpi"}).type == alps::mcoptions::MPI)) << "MPI execution option";
#ifndef ALPS_NGS_SINGLE_THREAD
    ASSERT_TRUE((parse<alps::mcoptions>({"sample.h5", "--mpi", "--threaded"}).type == alps::mcoptions::HYBRID)) << "Combined MPI and threaded execution";
#endif
#else
    rejects<alps::cli::error>({"sample.h5", "--mpi"},
            "Serial build rejects MPI execution");
#endif
}

// The grammar every ALPS program shares, shown on mcoptions.
TEST(CommandLine, OptionGrammar) {
    auto prefix = parse<alps::mcoptions>({"sample.h5", "--time", "5", "--cont"});
    ASSERT_TRUE((prefix.time_limit == 5 && prefix.resume)) << "A long name may be shortened to a unique prefix";
    rejects<alps::cli::error>({"sample.h5", "--c"}, "A prefix of both --continue and --checkpoint-file is rejected");
    rejects<alps::cli::error>({"sample.h5", "--Continue"}, "Option names are case-sensitive");

    auto equals = parse<alps::mcoptions>({"--input-file=sample.h5", "--time-limit=7"});
    ASSERT_TRUE((equals.input_file == "sample.h5" && equals.time_limit == 7)) << "--name=value";
    ASSERT_TRUE((parse<alps::mcoptions>({"sample.h5", "-T5"}).time_limit == 5)) << "-T5";
    auto grouped = parse<alps::mcoptions>({"sample.h5", "-cT5"});
    ASSERT_TRUE((grouped.resume && grouped.time_limit == 5)) << "Short flags group, and the last may take a value";

    auto ended = parse<alps::mcoptions>({"--", "-c", "--unknown"});
    ASSERT_TRUE((ended.input_file == "-c" && ended.output_file == "--unknown" && !ended.resume)) << "-- ends the options";

    rejects<alps::cli::error>({"sample.h5", "-c", "-c"}, "A repeated flag is rejected");
    rejects<alps::cli::error>({"sample.h5", "-T", "5", "--time-limit", "6"}, "A repeated value is rejected");
    rejects<alps::cli::error>({"sample.h5", "--input-file", "other.h5"}, "Input given both by position and by name is rejected");
    rejects<alps::cli::error>({"sample.h5", "-T"}, "A missing value is rejected");
    rejects<alps::cli::error>({"sample.h5", "-T", "-c"}, "A short option is not taken as a value");
    rejects<alps::cli::error>({"sample.h5", "-T", "-5"}, "A negative time limit is rejected");
    rejects<alps::cli::error>({"sample.h5", "--continue=1"}, "A flag takes no value");
}

TEST(CommandLine, LegacyParseargs) {
    auto defaults = parse<alps::parseargs>({"sample.in.h5"});
    ASSERT_TRUE((!defaults.resume && defaults.timelimit == 0 && defaults.tmin == 1 && defaults.tmax == 600)) << "parseargs scheduling defaults";
    // The two parsers intentionally retain different option and filename rules.
    ASSERT_TRUE((defaults.input_file == "sample.in.h5" && defaults.output_file == "sample.in.out.h5")) << "parseargs replaces only the final extension";
    ASSERT_TRUE((parse<alps::parseargs>({"sample"}).output_file == "sample.out.h5")) << "parseargs extensionless filename";
    auto positional = parse<alps::parseargs>({"input.xml", "output.h5", "-c", "-T", "120", "-i", "2", "-a", "30"});
    ASSERT_TRUE((positional.input_file == "input.xml" && positional.output_file == "output.h5"
            && positional.resume && positional.timelimit == 120 && positional.tmin == 2 && positional.tmax == 30)) << "parseargs positional filenames and short options";
    auto named = parse<alps::parseargs>({"--inputfile", "input.xml", "--outputfile", "output.h5",
            "--continue", "--timelimit", "17", "--Tmin", "3", "--Tmax", "40"});
    ASSERT_TRUE((named.input_file == "input.xml" && named.output_file == "output.h5"
            && named.resume && named.timelimit == 17 && named.tmin == 3 && named.tmax == 40)) << "parseargs named options retain their spelling";
    auto empty = parse<alps::parseargs>({});
    ASSERT_TRUE((empty.input_file.empty() && empty.output_file == ".out.h5")) << "parseargs continues to accept an empty command line";
}
}
