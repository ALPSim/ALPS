// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Wrote legacy_observables.h5 and legacy_observables.expected, once, with
// ALEA as it was when the statistics tests were added. It is kept as the
// record of how the file was made and is not run by any test: the point of the
// file is that its bytes no longer change.
//
//   statistics_write_legacy_observables <directory>
#include "legacy_observables.hpp"

#include <alps/hdf5.hpp>

#include <cstdio>
#include <string>

namespace {

void write(std::FILE * out, const char * group, const statistics_test::numbers & values) {
    // Hexadecimal floating point: every bit of each double.
    for (const auto & value : values)
        std::fprintf(out, "%s/%s %a\n", group, value.first.c_str(), value.second);
}

}

int main(int argc, char ** argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: %s <directory>\n", argv[0]);
        return 2;
    }
    const std::string directory = argv[1];
    alps::ObservableSet run = statistics_test::legacy_run(1001, 20000);
    alps::ObservableSet merged = statistics_test::legacy_run(1001, 20000, false);
    merged << statistics_test::legacy_run(2002, 12000, false);
    // Merged runs with a signed observable: ALEA writes them but cannot load
    // them back (see frozen_files.cpp), so no numbers are recorded for them.
    alps::ObservableSet merged_signed = statistics_test::legacy_run(1001, 20000);
    merged_signed << statistics_test::legacy_run(2002, 12000);
    // "w" opens an existing file for update, so start from nothing.
    std::remove((directory + "/legacy_observables.h5").c_str());
    {
        alps::hdf5::archive archive(directory + "/legacy_observables.h5", "w");
        archive["/run"] << run;
        archive["/merged"] << merged;
        archive["/merged_signed"] << merged_signed;
    }
    // Describe what a reader gets back, not what was in memory.
    alps::ObservableSet run_loaded, merged_loaded;
    {
        alps::hdf5::archive archive(directory + "/legacy_observables.h5", "r");
        archive["/run"] >> run_loaded;
        archive["/merged"] >> merged_loaded;
    }
    std::FILE * out = std::fopen((directory + "/legacy_observables.expected").c_str(), "w");
    if (!out)
        return 1;
    write(out, "run", statistics_test::describe_run(run_loaded));
    write(out, "merged", statistics_test::describe_merged(merged_loaded));
    return std::fclose(out);
}
