/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2008 by Lukas Gamper <mistral@student.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/config.h>

#include <string>
#include <iostream>
#include <stdexcept>
#include <fstream>

#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

#include "archive_sqlite.hpp"
#include "archive_index.hpp"
#include "archive_node.hpp"
#include "archive_xml.hpp"

#include "archive_plot.hpp"

int main(int ac, char* av[]) {
        SQLite db = SQLite();
        try {
                po::options_description generic("Generic options");
                generic.add_options()
                        ("version", "print version string")
                        ("verbose,v", "Prints what's done")
                        ("help,h", "produce help message");

                po::options_description config("Configuration");
                config.add_options()
                        ("db-file,d", po::value<std::string>()->default_value("archive.db"), "database file")
                        ("xml-file,x", po::value<std::string>(), "simulation XML file")
                        ("plot-file,p", po::value<std::string>(), "plot XML file")
                        ("output-path,o", po::value<std::string>()->default_value("."), "output path")
                        ("pattern-file,i", po::value<std::string>(), "index file, only used by installation and if compiled with -DUSEPATTERN")
                        ("full-list,l", "plots the whole list");

                po::options_description command("Commands");
                command.add_options()
                        ("command,c", po::value<std::string>(), "possible commands: 'help(h)', 'install(i)', 'plot'(p), 'append'(a), 'rebuilt'(r), 'list'(l)");

                po::options_description hidden("Hidden");
                hidden.add_options()
                        ("positional-arg", po::value<std::string>(), "positional argument");

                po::options_description cmdline_options;
                cmdline_options.add(generic).add(config).add(command).add(hidden);

                po::options_description config_file_options;
                config_file_options.add(config);

                po::options_description visible("Allowed options");
                visible.add(generic).add(config).add(command);

                po::positional_options_description posOpt;
                posOpt.add("positional-arg", 1);

                po::variables_map vm;
                store(po::command_line_parser(ac, av).options(cmdline_options).positional(posOpt).run(), vm);

                std::ifstream ifs("conf/config.cfg");
                po::store(po::parse_config_file(ifs, config_file_options), vm);
                po::notify(vm);

                if (vm.count("version")) {
                        std::cout << "ALPS Plot Generator, version 1.0" << std::endl;
                        exit(0);

                } else if (vm.count("help") || !vm.count("command") || vm["command"].as<std::string>() == "help" || vm["command"].as<std::string>() == "h") {
                        std::cout << visible << std::endl;
                        exit(0);

                } else if (!vm.count("db-file"))
                        throw std::runtime_error("no db-file specified!");

                else {
                        std::string const cmd = vm["command"].as<std::string>();

                        // A positional file stands in for -p (plot) or -x (append, rebuild).
                        // Resolve it once, so a file given both ways or to a command that
                        // reads none is an error instead of being silently ignored.
                        std::string positional;
                        if (vm.count("positional-arg")) {
                                positional = vm["positional-arg"].as<std::string>();
                                if (vm.count("plot-file") || vm.count("xml-file"))
                                        throw std::runtime_error("give the input file either as a positional argument or with -p/-x, not both");
                                if (cmd != "plot" && cmd != "p" && cmd != "append" && cmd != "a" && cmd != "rebuild" && cmd != "r")
                                        throw std::runtime_error(std::string("command '") + cmd + "' takes no input file");
                        }
                        auto inputFile = [&](char const* option, char const* what) {
                                if (positional.empty() && !vm.count(option))
                                        throw std::runtime_error(std::string("no ") + what + " specified");
                                return absolute(fs::path(positional.empty() ? vm[option].as<std::string>() : positional));
                        };

                        db.setVerbose(vm.count("verbose"));
                        db.open(absolute(fs::path(vm["db-file"].as<std::string>())));

                        if (cmd == "install" || cmd == "i") {
                                #ifdef USEPATTERN
                                        if (!vm.count("pattern-file"))
                                                throw std::runtime_error("no index-file specified!");
                                        Index(db, vm.count("verbose")).install(complte(fs::path(vm["pattern-file"].as<std::string>())));
                                #else
                                        Index(db, vm.count("verbose")).install();
                                #endif

                        } else if (cmd == "list" || cmd == "l") {
                                Index(db, vm.count("verbose")).list(vm.count("full-list"));

                        } else if (cmd == "plot" || cmd == "p") {
                                fs::path plotPath = inputFile("plot-file", "plot-file");
                                if (!vm.count("output-path"))
                                        throw std::runtime_error("no output-path specified");
                                fs::path outputPath = absolute(fs::path(vm["output-path"].as<std::string>()));
                                Plot(outputPath, db).exec(XML(true)(plotPath, true), plotPath.filename().string());

                        } else {
                                fs::path xmlPath = inputFile("xml-file", "xml-file");
                                if (cmd == "rebuild" || cmd == "r") {
                                        if(db.clear() == false)
                                                throw std::runtime_error("Could not clear tables!");
                                } else if (cmd != "append" && cmd != "a")
                                        throw std::runtime_error(std::string("Unknown command '") + cmd + std::string("'"));
                                Index(db, vm.count("verbose")).exec(xmlPath);
                        }
                }
        } catch(std::exception& e) {
                std::cerr << e.what() << std::endl;
                db.close(false);
                return 1;
        }
        db.close(true);
        return 0;
}
