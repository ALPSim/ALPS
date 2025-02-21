/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2006 by Lukas Gamper <mistral@student.ethz.ch>,
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

#ifndef _SQLITE_H_
#define _SQLITE_H_

#include <alps/config.h>

#include <string>
#include <list>
#include <map>
#include <stdexcept>
#include <iostream>
#include <ctime>
#include <boost/filesystem/path.hpp>
#include <sqlite3.h>

namespace fs = boost::filesystem;

class SQLite {
    sqlite3 *mDB;
        char *mErrMsg;
        int mCnt;
        bool mVerbose;
        double timer;
        std::list<std::map<std::string, std::string> > mRS;

        public:

                /**
                 * Defaultconstructor
                 * 
                 * @param verbose If verbose ist true, the database says what's doing
                 */
                SQLite(bool verbose = true): mDB(NULL), mCnt(0), mVerbose(verbose), timer(0) {}

                /**
                 * Constructor to initialize the database. If the Database given in filename not exist, a new database will 
                 * be created.
                 * 
                 * @param filename Path to the file that contains the database.
                 * @param verbose If verbose ist true, the database says what's doing
                 */
            SQLite(fs::path filename, bool verbose = true): mCnt(0), mVerbose(verbose), timer(0) { open(filename); }

                /**
                 * Sets the verbosemode of the Database.
                 * 
                 * @param verbose If verbose ist true, the database says what's doing
                 */
                void setVerbose(bool verbose = true) { mVerbose = verbose; }

                /**
                 * Return the last Errormessage of the Database. If no Error accures a empty String will be returned.
                 * 
                 * @return Last errormsg of the database
                 */
                std::string getErrorMsg();

                /**
                 * Open the database and starts the transaction. If the File not exist, a new Database will be created.
                 * 
                 * @param filename Path to the database file.
                 */
            void open(fs::path filename);

                /**
                 * Close the database and commit or rollback the transaction
                 * 
                 * @param commit Commit on true, rollback on false
                 */
            void close(bool commit);

                /**
                 * Deletes the Content of the Databases
                 */
                bool clear();

                /**
                 * Internal Function to communicate with the C-interface of the sqlite. This function should not
                 * used from outside!!
                 */
                void addTupel(std::map<std::string, std::string>);

                /**
                 * Replace all ' by &apos; to avoid problems in queries.
                 * 
                 * @param str String to Quote
                 * @return Quoted string
                 */
                static std::string quote(std::string str);

                /**
                 * Replace all &apos; by ' make results von database readable.
                 * 
                 * @param str String to be unquoted
                 * @return Unquoted string
                 */
                static std::string unQuote(std::string str);

                /**
                 * Send a Query to the database. The Result is saved as list<name => value>
                 * 
                 * @param query Query to send to the database
                 * @return Results of the query
                 */
                std::list<std::map<std::string, std::string> > operator()(std::string query);

                /**
                 * Send a Query to the database. The Result is saved as list<name => value>. If the Query fails,
                 * and throwError ist true an Error is thrown.
                 * 
                 * @param query Query to send to the database
                 * @param throwError days it an error should be thrown
                 * @return Results of the query
                 */
                std::list<std::map<std::string, std::string> > operator()(std::string query, bool throwError);
};
#endif //_SQLITE_H_
