/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2008 by Lukas Gamper <mistral@student.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
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
                cmdline_options.add(generic).add(config).add(command);

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

                if (vm.count("help") || !vm.count("command") || vm["command"].as<std::string>() == "help" || vm["command"].as<std::string>() == "h") {
                        std::cout << visible << std::endl;
                        exit(0);

                } else if (vm.count("version")) {
                        std::cout << "ALPS Plot Generator, version 1.0" << std::endl;
                        exit(0);

                } else if (!vm.count("db-file"))
                        throw std::runtime_error("no db-file specified!");

                else {
                        db.setVerbose(vm.count("verbose"));
                        db.open(absolute(fs::path(vm["db-file"].as<std::string>())));

                        if (vm["command"].as<std::string>() == "install" || vm["command"].as<std::string>() == "i") {
                                #ifdef USEPATTERN
                                        if (!vm.count("pattern-file"))
                                                throw std::runtime_error("no index-file specified!");
                                        Index(db, vm.count("verbose")).install(complte(fs::path(vm["pattern-file"].as<std::string>())));
                                #else
                                        Index(db, vm.count("verbose")).install();
                                #endif

                        } else if (vm["command"].as<std::string>() == "list" || vm["command"].as<std::string>() == "l") {
                                Index(db, vm.count("verbose")).list(vm.count("full-list"));

                        } else if (vm["command"].as<std::string>() == "plot" || vm["command"].as<std::string>() == "p") {
                                if (!vm.count("plot-file") && !vm.count("positional-arg"))
                                        throw std::runtime_error("no plot-file specified");
                                else if (!vm.count("output-path"))
                                        throw std::runtime_error("no output-path specified");
                                fs::path plotPath = absolute(fs::path(vm["plot-file"].as<std::string>()));
                                if (vm.count("positional-arg"))
                                  plotPath = absolute(fs::path(vm["positional-arg"].as<std::string>()));
                                fs::path outputPath = absolute(fs::path(vm["output-path"].as<std::string>()));
                                Plot(outputPath, db).exec(XML(true)(plotPath, true), plotPath.filename().string());

                        } else {
                                if (!vm.count("xml-file") && !vm.count("positional-arg"))
                                        throw std::runtime_error("No xml-file specified!");
                                else if (vm["command"].as<std::string>() == "rebuild" || vm["command"].as<std::string>() == "r") {
                                        if(db.clear() == false)
                                                throw std::runtime_error("Could not clear tables!");
                                } else if (vm["command"].as<std::string>() != "append" && vm["command"].as<std::string>() != "a")
                                        throw std::runtime_error(std::string("Unknown command '") + vm["command"].as<std::string>() + std::string("'"));
                                Index(db, vm.count("verbose")).exec(absolute(fs::path(vm["xml-file"].as<std::string>())));
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
