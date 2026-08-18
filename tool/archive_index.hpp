/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2009 by Lukas Gamper <mistral@student.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
