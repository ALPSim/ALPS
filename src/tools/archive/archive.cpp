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

#include <alps/cli/options.hpp>
#include <filesystem>
#include <map>

namespace fs = std::filesystem;

#include "archive_sqlite.hpp"
#include "archive_index.hpp"
#include "archive_node.hpp"
#include "archive_xml.hpp"

#include "archive_plot.hpp"

namespace {

// Settings from a file of "name = value" lines, where # starts a comment.
// Only the Configuration options can be set this way.
std::map<std::string, std::string> readConfigFile(char const* file) {
        std::map<std::string, std::string> settings;
        auto trim = [](std::string const& text) {
                std::string::size_type first = text.find_first_not_of(" \t\r");
                return first == std::string::npos ? std::string() : text.substr(first, text.find_last_not_of(" \t\r") - first + 1);
        };
        std::ifstream in(file);
        for (std::string line; std::getline(in, line);) {
                line = trim(line.substr(0, line.find('#')));
                std::string::size_type equals = line.find('=');
                if (line.empty())
                        continue;
                if (equals == std::string::npos)
                        throw std::runtime_error(std::string(file) + " contains an invalid line '" + line + "'");
                std::string name = trim(line.substr(0, equals));
                if (name != "db-file" && name != "xml-file" && name != "plot-file" && name != "output-path"
                                && name != "pattern-file" && name != "full-list")
                        throw std::runtime_error(std::string(file) + ": unrecognised option '" + name + "'");
                if (!settings.emplace(name, trim(line.substr(equals + 1))).second)
                        throw std::runtime_error(std::string(file) + ": option '" + name + "' cannot be specified more than once");
        }
        return settings;
}

}

int main(int ac, char* av[]) {
        SQLite db = SQLite();
        try {
                std::string dbFile, outputPath;
                alps::cli::options options("Generic options");
                options.flag("version", "print version string")
                        .flag("verbose,v", "Prints what's done")
                        .flag("help,h", "produce help message")
                        .section("Configuration")
                        .value("db-file,d", dbFile, "archive.db", "database file")
                        .value<std::string>("xml-file,x", "simulation XML file")
                        .value<std::string>("plot-file,p", "plot XML file")
                        .value("output-path,o", outputPath, ".", "output path")
                        .value<std::string>("pattern-file,i", "index file, only used by installation and if compiled with -DUSEPATTERN")
                        .flag("full-list,l", "plots the whole list")
                        .section("Commands")
                        .value<std::string>("command,c", "possible commands: 'help(h)', 'install(i)', 'plot'(p), 'append'(a), 'rebuilt'(r), 'list'(l)")
                        .value<std::string>("positional-arg", "positional argument").hidden()
                        .positional("positional-arg");
                options.parse(ac, av);

                // The command line takes precedence over conf/config.cfg.
                std::map<std::string, std::string> const config = readConfigFile("conf/config.cfg");
                auto given = [&](char const* name) { return options.count(name) || config.count(name); };
                auto value = [&](char const* name) {
                        return options.count(name) ? options.get<std::string>(name) : config.at(name);
                };
                if (!options.count("db-file") && config.count("db-file"))
                        dbFile = config.at("db-file");
                if (!options.count("output-path") && config.count("output-path"))
                        outputPath = config.at("output-path");

                if (options.count("version")) {
                        std::cout << "ALPS Plot Generator, version 1.0" << std::endl;
                        exit(0);

                } else if (options.count("help") || !options.count("command") || value("command") == "help" || value("command") == "h") {
                        std::cout << "Allowed options:\n\n" << options << std::endl;
                        exit(0);

                } else {
                        std::string const cmd = value("command");

                        // A positional file stands in for -p (plot) or -x (append, rebuild).
                        // Resolve it once, so a file given both ways or to a command that
                        // reads none is an error instead of being silently ignored.
                        std::string positional;
                        if (options.count("positional-arg")) {
                                positional = value("positional-arg");
                                if (given("plot-file") || given("xml-file"))
                                        throw std::runtime_error("give the input file either as a positional argument or with -p/-x, not both");
                                if (cmd != "plot" && cmd != "p" && cmd != "append" && cmd != "a" && cmd != "rebuild" && cmd != "r")
                                        throw std::runtime_error(std::string("command '") + cmd + "' takes no input file");
                        }
                        auto inputFile = [&](char const* option, char const* what) {
                                if (positional.empty() && !given(option))
                                        throw std::runtime_error(std::string("no ") + what + " specified");
                                return absolute(fs::path(positional.empty() ? value(option) : positional));
                        };

                        db.setVerbose(options.count("verbose"));
                        db.open(absolute(fs::path(dbFile)));

                        if (cmd == "install" || cmd == "i") {
                                #ifdef USEPATTERN
                                        if (!given("pattern-file"))
                                                throw std::runtime_error("no index-file specified!");
                                        Index(db, options.count("verbose")).install(complte(fs::path(value("pattern-file"))));
                                #else
                                        Index(db, options.count("verbose")).install();
                                #endif

                        } else if (cmd == "list" || cmd == "l") {
                                Index(db, options.count("verbose")).list(given("full-list"));

                        } else if (cmd == "plot" || cmd == "p") {
                                fs::path plotPath = inputFile("plot-file", "plot-file");
                                Plot(absolute(fs::path(outputPath)), db).exec(XML(true)(plotPath, true), plotPath.filename().string());

                        } else {
                                fs::path xmlPath = inputFile("xml-file", "xml-file");
                                if (cmd == "rebuild" || cmd == "r") {
                                        if(db.clear() == false)
                                                throw std::runtime_error("Could not clear tables!");
                                } else if (cmd != "append" && cmd != "a")
                                        throw std::runtime_error(std::string("Unknown command '") + cmd + std::string("'"));
                                Index(db, options.count("verbose")).exec(xmlPath);
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
