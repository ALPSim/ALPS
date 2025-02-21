/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005 by Lukas Gamper <mistral@student.ethz.ch>
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
