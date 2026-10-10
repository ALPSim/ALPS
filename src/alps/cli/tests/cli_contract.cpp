// SPDX-License-Identifier: MIT
#include <alps/ngs/mcoptions.hpp>
#include <alps/parseargs.hpp>

#include <boost/program_options/errors.hpp>

#include <initializer_list>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition, char const* message) {
    if (!condition) throw std::runtime_error(message);
}

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
    require(rejected, message);
}

struct capture_stdout {
    std::ostringstream text;
    std::streambuf* previous = std::cout.rdbuf(text.rdbuf());
    ~capture_stdout() { std::cout.rdbuf(previous); }
};

void mcoptions_contract() {
    auto defaults = parse<alps::mcoptions>({"sample.in.h5"});
    require(defaults.valid && !defaults.resume && defaults.time_limit == 0,
            "mcoptions default validity, resume and time limit");
    require(defaults.type == alps::mcoptions::SINGLE && defaults.checkpoint_file.empty(),
            "mcoptions default execution and checkpoint");
    require(defaults.input_file == "sample.in.h5" && defaults.output_file == "sample.out.h5",
            "mcoptions default input/output names");
    require(parse<alps::mcoptions>({"sample.out.h5"}).output_file == "sample.out.h5",
            "Existing output checkpoint keeps its name");
    require(parse<alps::mcoptions>({"sample.h5"}).output_file == "sample.out.h5",
            "General input extension is replaced");
    require(parse<alps::mcoptions>({"sample"}).output_file == "sample.out.h5",
            "Extensionless input gets an output suffix");
    // These are existing substring rules, not a request to normalize filenames.
    require(parse<alps::mcoptions>({"sample.out.h5.backup"}).output_file == "sample.out.h5.backup",
            "Existing output-name substring behavior");

    auto positional = parse<alps::mcoptions>({"input.h5", "output.h5", "checkpoint.h5", "-c", "-T", "120"});
    require(positional.input_file == "input.h5" && positional.output_file == "output.h5"
            && positional.checkpoint_file == "checkpoint.h5" && positional.resume
            && positional.time_limit == 120, "mcoptions positional filenames and short options");
    auto named = parse<alps::mcoptions>({"--input-file", "input.h5", "--output-file", "output.h5",
            "--checkpoint-file", "checkpoint.h5", "--continue", "--time-limit", "17", "--single"});
    require(named.input_file == "input.h5" && named.output_file == "output.h5"
            && named.checkpoint_file == "checkpoint.h5" && named.resume && named.time_limit == 17
            && named.type == alps::mcoptions::SINGLE, "mcoptions named options");

    {
        capture_stdout output;
        auto help = parse<alps::mcoptions>({"--help"});
        require(!help.valid && help.input_file.empty(), "Help does not require an input file");
        require(output.text.str().find("--input-file") != std::string::npos
                && output.text.str().find("--time-limit") != std::string::npos,
                "Help describes the current grammar");
    }
    rejects<std::invalid_argument>({}, "Missing input is rejected");
    rejects<boost::program_options::unknown_option>({"sample.h5", "--unknown"},
            "Unknown mcoptions option is rejected");
    rejects<boost::program_options::invalid_option_value>({"sample.h5", "-T", "not-a-number"},
            "Nonnumeric time limit is rejected");
    rejects<boost::program_options::error>({"a", "b", "c", "d"},
            "Extra positional filename is rejected");

#ifndef ALPS_NGS_SINGLE_THREAD
    require(parse<alps::mcoptions>({"sample.h5", "--threaded"}).type == alps::mcoptions::THREADED,
            "Threaded execution option");
#else
    rejects<boost::program_options::unknown_option>({"sample.h5", "--threaded"},
            "Single-thread build rejects threaded execution");
#endif
#ifdef ALPS_HAVE_MPI
    require(parse<alps::mcoptions>({"sample.h5", "--mpi"}).type == alps::mcoptions::MPI,
            "MPI execution option");
#ifndef ALPS_NGS_SINGLE_THREAD
    require(parse<alps::mcoptions>({"sample.h5", "--mpi", "--threaded"}).type == alps::mcoptions::HYBRID,
            "Combined MPI and threaded execution");
#endif
#else
    rejects<boost::program_options::unknown_option>({"sample.h5", "--mpi"},
            "Serial build rejects MPI execution");
#endif
}

void parseargs_contract() {
    auto defaults = parse<alps::parseargs>({"sample.in.h5"});
    require(!defaults.resume && defaults.timelimit == 0 && defaults.tmin == 1 && defaults.tmax == 600,
            "parseargs scheduling defaults");
    // The two parsers intentionally retain different option and filename rules.
    require(defaults.input_file == "sample.in.h5" && defaults.output_file == "sample.in.out.h5",
            "parseargs replaces only the final extension");
    require(parse<alps::parseargs>({"sample"}).output_file == "sample.out.h5",
            "parseargs extensionless filename");
    auto positional = parse<alps::parseargs>({"input.xml", "output.h5", "-c", "-T", "120", "-i", "2", "-a", "30"});
    require(positional.input_file == "input.xml" && positional.output_file == "output.h5"
            && positional.resume && positional.timelimit == 120 && positional.tmin == 2 && positional.tmax == 30,
            "parseargs positional filenames and short options");
    auto named = parse<alps::parseargs>({"--inputfile", "input.xml", "--outputfile", "output.h5",
            "--continue", "--timelimit", "17", "--Tmin", "3", "--Tmax", "40"});
    require(named.input_file == "input.xml" && named.output_file == "output.h5"
            && named.resume && named.timelimit == 17 && named.tmin == 3 && named.tmax == 40,
            "parseargs named options retain their spelling");
    auto empty = parse<alps::parseargs>({});
    require(empty.input_file.empty() && empty.output_file == ".out.h5",
            "parseargs continues to accept an empty command line");
}
}

int main() {
    try {
        mcoptions_contract();
        parseargs_contract();
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
