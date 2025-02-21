/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2009 by Lukas Gamper <mistral@student.ethz.ch>,
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
