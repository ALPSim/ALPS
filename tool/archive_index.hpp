/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2009 by Lukas Gamper <mistral@student.ethz.ch>,
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

#ifndef _INDEX_H_
#define _INDEX_H_

#include "archive_sqlite.hpp"
#include <boost/filesystem/operations.hpp>

namespace fs = boost::filesystem;

class Index {
        SQLite &mDB;
        bool mVerbose;
        void cretateTables();
        bool patternFilter(std::string, std::string);

        public:
                Index(SQLite &inDB, bool verbose = true): mDB(inDB), mVerbose(verbose) {}

                #ifdef USEPATTERN
                        void install(fs::path inPatternFile);
                #else
                        void install();
                #endif

                /**
                 * Sets a reference to the database
                 *
                 * @param &inDB Reference of the database
                 */
                void setDB(SQLite &inDB) { mDB = inDB; }

                /**
                 * Lists all the parameters and measurements indexed in the database
                 *
                 * @param inFullList Should the whole list be displaied
                 */
                void list(bool inFullList);

                /**
                 * Adds all files in the directory and subdirectories to the index it they do not already exit.
                 *
                 * @param xmlPath Path to scan
                 */
                void exec(fs::path xmlPath);
};

#endif //_INDEX_H_
