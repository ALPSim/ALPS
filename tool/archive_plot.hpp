/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005 by Lukas Gamper <mistral@student.ethz.ch>
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
