// SPDX-License-Identifier: MIT
#include <alps/ngs/mcoptions.hpp>
#include <alps/parseargs.hpp>
#include <stdexcept>

int main() {
    char executable[] = "consumer";
    char input[] = "sample.in.h5";
    char output[] = "chosen.out.h5";
    char resume[] = "--continue";
    char* argv[] = {executable, input, output, resume};
    alps::mcoptions current(4, argv);
    alps::parseargs older(4, argv);
    if (!current.valid || !current.resume || !older.resume
        || current.input_file != input || older.input_file != input
        || current.output_file != output || older.output_file != output)
        throw std::runtime_error("Installed command-line component contract failed");
}
