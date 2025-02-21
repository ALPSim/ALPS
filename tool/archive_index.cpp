/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2006-2013 by Lukas Gamper <mistral@student.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

/* $Id$ */

#include <alps/config.h>

#include <fstream>
#include <string>
#include <list>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <boost/lexical_cast.hpp>

#include "archive_index.hpp"

#include "archive_node.hpp"
#include "archive_xml.hpp"

#include <alps/parapack/parapack.h>

#ifdef USEPATTERN
        #include <stdexcept>
        #include <stack>

        void Index::install(fs::path inPatternFile) {
                mDB("COMMIT;", true);
                mDB("CREATE TABLE pattern (name TEXT PRIMARY KEY, type TEXT, value TEXT);", true);
                cretateTables();
                mDB("BEGIN;", true);
                Node xmlRoot = XML(mVerbose)(inPatternFile);
                std::list<Node> parameters = xmlRoot.getElement("pattern").getElement("parameters").nodeTest("parameter");
                std::list<Node> measurements = xmlRoot.getElement("pattern").getElement("measurements").nodeTest("measurement");
                for (std::list<Node>::iterator it = parameters.begin(); it != parameters.end(); it++) {
                        std::string name = it->getAttribute("name");
                        if (name.empty())
                                name = it->string();
                        mDB(std::string("INSERT INTO pattern (name, type, value) VALUES ('") + mDB.quote(name) + "', 'parameter', '" + mDB.quote(it->string()) + "');", true);
                }
                for (std::list<Node>::iterator it = measurements.begin(); it != measurements.end(); it++) {
                        std::string name = it->getAttribute("name");
                        if (name.empty())
                                name = it->string();
                        mDB(std::string("INSERT INTO pattern (name, type, value) VALUES ('") + mDB.quote(name) + "', 'measurement', '" + mDB.quote(it->string()) + "');", true);
                }
                std::cout << "Patterns read" << std::endl;
        }

        bool Index::patternFilter(std::string type, std::string name) {
                return (mDB(std::string("SELECT * FROM pattern WHERE type='") + type + "' AND value='" + mDB.quote(name) +"';", true).empty() == false);
        }

#else
        #include <stdexcept>
        #include <stack>
        #include <ctime>

        void Index::install() {
                mDB("COMMIT;", true);
                cretateTables();
                mDB("BEGIN;", true);
        }

        bool Index::patternFilter(std::string, std::string) {
                return true;
        }
#endif

void Index::cretateTables() {
        mDB("CREATE TABLE uri (path TEXT UNIQUE, fID INTEGER PRIMARY KEY AUTOINCREMENT);", true);
        mDB("CREATE TABLE parameter (fID INTEGER, name TEXT, value TEXT, PRIMARY KEY (fID, name));", true);
        mDB("CREATE TABLE measurement (fID INTEGER, indexvalue TEXT, name TEXT, count TEXT, mean TEXT, error TEXT, variance TEXT, autocorr TEXT, PRIMARY KEY (fID, indexvalue, name));", true);
        mDB("CREATE VIEW cnt AS SELECT name FROM measurement GROUP BY name, fID;", true);
        mDB("PRAGMA synchronous=OFF;", true);
        std::cout << "Tables created" << std::endl;
}

void Index::list(bool inFullList) {
        std::list<std::map<std::string, std::string> > rs = mDB("SELECT name, count(name) AS count FROM parameter GROUP BY name ORDER BY name", true);
        if (rs.empty())
                std::cout << "No parameter avalable!";
        else {
                std::cout << std::setw(7) << "Total" << " Parameters {Values(Occurrence)*}" << std::endl;
                for (std::list<std::map<std::string, std::string> >::iterator it = rs.begin(); it != rs.end(); it++) {
                        std::cout << std::setw(7) << (*it)["count"] << " " << mDB.unQuote((*it)["name"]) << " {";
                        std::list<std::map<std::string, std::string> > subRS = mDB("SELECT value, count(value) AS count FROM parameter WHERE name='"
                                + (*it)["name"] + "' GROUP BY value ORDER BY value", true);
                        unsigned int count = 0;
                        for (std::list<std::map<std::string, std::string> >::iterator subIT = subRS.begin(); subIT != subRS.end(); subIT++, count ++) {
                                if (!inFullList && count == 6 && subRS.size() > 10)
                                        std::cout << ", ...";
                                if (inFullList || (count < 5 || count > subRS.size() - 6)) {
                                        if (subIT != subRS.begin())
                                                std::cout << ", ";
                                        std::cout << (*subIT)["value"] << "(" << (*subIT)["count"] << ")";
                                }
                        }
                        std::cout << "}" << std::endl;
                }
        }
        rs = mDB("SELECT name, count(name) AS count FROM measurement GROUP BY name ORDER BY name", true);
        if (rs.empty())
                std::cout << std::endl << "No measurements avalable!" << std::endl;
        else {
                std::cout << std::endl << "Measurements" << std::endl;
                for (std::list<std::map<std::string, std::string> >::iterator it = rs.begin(); it != rs.end(); it++)
                        std::cout << std::setw(6) << (*it)["count"] << " " << mDB.unQuote((*it)["name"]) << std::endl;
        }
}

void Index::exec(fs::path xmlPath) {
        if (!fs::exists(xmlPath))
                throw std::runtime_error(std::string("Not found: ") + xmlPath.string());
        fs::path basedir = xmlPath.parent_path();
        std::string file_in_str;
        std::string file_out_str;
        std::vector<fs::path> files;
        int t = alps::parapack::load_filename(xmlPath, file_in_str, file_out_str);
        if (t == 1) {
          std::vector<alps::task> tasks;
          fs::path file_in = absolute(fs::path(file_in_str), basedir);
          fs::path file_out = absolute(fs::path(file_out_str), basedir);
          std::string simname;
          alps::parapack::load_tasks(file_in, file_out, basedir, simname, tasks);
          BOOST_FOREACH(alps::task& t, tasks) {
            fs::path f = absolute(fs::path(t.file_out_str()), basedir);
            if (fs::exists(f)) {
              files.push_back(f);
            } else {
              std::cerr << "Warning: " << f.string() << " not found\n";
            }
          }
        } else {
          files.push_back(xmlPath);
        }

        long fileCnt = 0;
        #ifdef DEBUG
                std::clock_t timer;
        #endif
        std::clock_t fulltime = std::clock();
        BOOST_FOREACH(fs::path const& path, files) {
                if (mDB(std::string("SELECT * FROM uri WHERE path='") + path.string() + "';").empty()) {
                        #ifdef DEBUG
                                timer = std::clock();
                        #endif
                        Node xmlRoot = XML(mVerbose)(path, false);
                        #ifdef DEBUG
                                timer = std::clock() - timer;
                        #endif
                        std::list<Node> parameters = xmlRoot.getElement("simulation").getElement("parameters").nodeTest("parameter");
                        std::list<Node> averages = xmlRoot.getElement("simulation").nodeTest("averages");
                        if (parameters.size() == 0) {
                                std::cerr << "The file '" << path.string() << "' is not a valid XML document!" << std::endl;
                        } else {
                          int aidx = 0; // average index
                          for (std::list<Node>::iterator av = averages.begin(); av != averages.end(); ++av, ++aidx) {
                            std::string pathname = path.string();
                            if (aidx > 0) pathname += "-" + boost::lexical_cast<std::string>(aidx);
                            mDB(std::string("INSERT INTO uri (path) VALUES ('") + pathname + "');", true);
                            std::string fID = mDB(std::string("SELECT fID FROM uri WHERE path='") + pathname + "';", true).front()["fID"];
                            std::list<Node> scalarAverage = av->nodeTest("scalar_average");
                            std::list<Node> vectorAverage = av->nodeTest("vector_average");
                            std::list<Node> histogram = av->nodeTest("histogram");
                            mDB(std::string("INSERT INTO parameter (fID, name, value) VALUES ('") + fID + "', '" + mDB.quote("average index") +
                                "', '" + mDB.quote(boost::lexical_cast<std::string>(aidx)) + "');", true);
                            for (std::list<Node>::iterator it = parameters.begin(); it != parameters.end(); ++it)
                              if (patternFilter("parameter", it->getAttribute("name")))
                                mDB(std::string("INSERT INTO parameter (fID, name, value) VALUES ('") + fID + "', '" + mDB.quote(it->getAttribute("name")) +
                                    "', '" + mDB.quote(it->string()) + "');", true);
                            for (std::list<Node>::iterator it = scalarAverage.begin(); it != scalarAverage.end(); it++)
                              if (patternFilter("measurement", it->getAttribute("name")))
                                mDB(std::string("INSERT INTO measurement (fID, indexvalue, name, count, mean, error, variance, autocorr) ")
                                    + "VALUES ('" + fID + "', '', '"
                                    + mDB.quote(it->getAttribute("name")) + "', '"
                                    + mDB.quote(it->getElement("count").string()) + "', '"
                                    + mDB.quote(it->getElement("mean").string()) + "', '"
                                    + mDB.quote(it->getElement("error").string()) + "', '"
                                    + mDB.quote(it->getElement("variance").string()) + "', '"
                                    + mDB.quote(it->getElement("autocorr").string()) + "')", true);
                            for (std::list<Node>::iterator vIt = vectorAverage.begin(); vIt != vectorAverage.end(); vIt++)
                              if (patternFilter("measurement", vIt->getAttribute("name"))) {
                                std::list<Node> scalarAvg = vIt->nodeTest("scalar_average");
                                for (std::list<Node>::iterator it = scalarAvg.begin(); it != scalarAvg.end(); it++)
                                  mDB(std::string("INSERT INTO measurement (fID, indexvalue, name, count, mean, error, variance, autocorr) ")
                                      + "VALUES ('" + fID + "', '"
                                      + it->getAttribute("indexvalue") + "', '"
                                      + mDB.quote(vIt->getAttribute("name")) + "', '"
                                      + mDB.quote(it->getElement("count").string()) + "', '"
                                      + mDB.quote(it->getElement("mean").string()) + "', '"
                                      + mDB.quote(it->getElement("error").string()) + "', '"
                                      + mDB.quote(it->getElement("variance").string()) + "', '"
                                      + mDB.quote(it->getElement("autocorr").string()) + "')", true);
                              }
                            for (std::list<Node>::iterator vIt = histogram.begin(); vIt != histogram.end(); vIt++)
                              if (patternFilter("measurement", vIt->getAttribute("name"))) {
                                std::list<Node> entry = vIt->nodeTest("entry");
                                for (std::list<Node>::iterator it = entry.begin(); it != entry.end(); it++) {
                                  uint64_t count = boost::lexical_cast<uint64_t>(it->getElement("count").string());
                                  double value = boost::lexical_cast<double>(it->getElement("value").string());
                                  double mean = value / count;
                                  double error = std::sqrt(value) / count;
                                  mDB(std::string("INSERT INTO measurement (fID, indexvalue, name, count, mean, error, variance, autocorr) ")
                                      + "VALUES ('" + fID + "', '"
                                      + it->getAttribute("indexvalue") + "', '"
                                      + mDB.quote(vIt->getAttribute("name")) + "', '"
                                      + mDB.quote(boost::lexical_cast<std::string>(count)) + "', '"
                                      + mDB.quote(boost::lexical_cast<std::string>(mean)) + "', '"
                                      + mDB.quote(boost::lexical_cast<std::string>(error)) + "', '', '')", true);
                                }
                              }
                          }
                          fileCnt++;
                          if (mVerbose && fileCnt % 10 == 0)
                            std::cout << fileCnt << " Files scanned" << std::endl;
                          #ifdef DEBUG
                          if (mVerbose)
                            std::cout << "Index created on " << path.string() << " (" << std::setiosflags(std::ios_base::fixed) << std::setprecision (2) << std::difftime(std::clock(), timer)/1000000 << "s)" << std::endl;
                          #endif
                        }
                } else {
                  std::cout << "File " << path.string() << " already indexed" << std::endl;
                }
        }
        // if (mVerbose)
                std::cout << fileCnt << " files scanned in " << std::setiosflags(std::ios_base::fixed) << std::setprecision (2) << std::difftime(std::clock(), fulltime) / 1000000 << "s" << std::endl;
}
