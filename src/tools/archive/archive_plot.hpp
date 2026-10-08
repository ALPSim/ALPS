/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005 by Lukas Gamper <mistral@student.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef _PLOT_HPP_
#define _PLOT_HPP_

#include <alps/config.h>

#include <string>
#include <map>
#include <list>
#include <vector>
#include <boost/filesystem/path.hpp>

#include "archive_sqlite.hpp"
#include "archive_node.hpp"

namespace fs = boost::filesystem;

class Plot {
        SQLite &mDB;
        fs::path mOutPath;
        bool mVerbose;

        std::string strToLower(std::string inStr);
        void writeFile(fs::path inOutFile, std::string inBuffer);

        public:
                Plot(fs::path inOutPath, SQLite &inDB, bool verbose = true): mDB(inDB), mOutPath(inOutPath), mVerbose(verbose) {}
                void setOutPath(fs::path inOutPath) { mOutPath = inOutPath; }
                void setDB(SQLite &inDB) { mDB = inDB; }
                void exec(Node inNode, std::string inInFile);
};

#endif //_PLOT_HPP_
