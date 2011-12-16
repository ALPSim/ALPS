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

#include <iomanip>
#include "archive_sqlite.hpp"

int SQLiteCallback(void *db, int argc, char **argv, char **azColName){
        std::map<std::string, std::string> tupel;
        for(int loop = 0; loop < argc; loop++)
                tupel[azColName[loop]] = (argv[loop] ? argv[loop] : "NULL");
        ((SQLite *)db)->addTupel(tupel);
        return 0;
}

void SQLite::addTupel(std::map<std::string, std::string> tupel) {
        mRS.push_back(tupel);
}

std::string SQLite::getErrorMsg() {
        return std::string(mErrMsg);
}

void SQLite::open(fs::path filename) {
        #ifdef DEBUG
                if (mVerbose)
                        std::cout << "Open database '" << filename.string() << "'" << std::endl;
        #endif
        if(sqlite3_open(filename.string().c_str(), &mDB))
                throw std::runtime_error(std::string("Can't open database: ") + sqlite3_errmsg(mDB));
        operator()("BEGIN;");
 }

void SQLite::close(bool commit) {
        #ifdef DEBUG
                if(mCnt && mVerbose)
                        std::cout << "Number of executed queries: " << mCnt << " (" << std::setiosflags(std::ios_base::fixed) << std::setprecision(2) << timer / 1000 << "s)" << std::endl;
        #endif
         if(mDB) {
                if (commit)
                        operator()("COMMIT;");
                else
                        operator()("ROLLBACK;");
                #ifdef DEBUG
                        if (mVerbose)
                                std::cout << "Close database" << std::endl;
                #endif
                 sqlite3_close(mDB);
         } else if (mVerbose) {
                #ifdef DEBUG
                         std::cerr << "Databse not opened!" << std::endl;
                 #endif
         }
}

bool SQLite::clear() {
        return sqlite3_exec(mDB, "DELETE FROM uri;", NULL, 0, &mErrMsg) == SQLITE_OK &&
                sqlite3_exec(mDB, "DELETE FROM parameter;", NULL, 0, &mErrMsg) == SQLITE_OK &&
                sqlite3_exec(mDB, "DELETE FROM measurement;", NULL, 0, &mErrMsg) == SQLITE_OK;
}

std::string SQLite::quote(std::string str) {
        std::string::size_type pos = str.find("'");
        if (pos != std::string::npos)
                str.replace(pos, 1, "&apos;" );
        return str;
}

std::string SQLite::unQuote(std::string str) {
        std::string::size_type pos = str.find("&apos;");
        if (pos != std::string::npos)
                str.replace(pos, 6, "'" );
        return str;
}

std::list<std::map<std::string, std::string> > SQLite::operator()(std::string query) {
        return operator()(query, false);
}

std::list<std::map<std::string, std::string> > SQLite::operator()(std::string query, bool throwError) {
        mRS = std::list<std::map<std::string, std::string> >();
        mCnt++;
        clock_t tmp = std::clock();
        if (sqlite3_exec(mDB, query.c_str(), SQLiteCallback, this, &mErrMsg) != SQLITE_OK && throwError)
                throw std::runtime_error(std::string("Error in Query '") + query + "': " + sqlite3_errmsg(mDB));
        timer += std::difftime(std::clock(), tmp);
        if (mVerbose) std::clog << "Exectured query " << query << "\nGot " << mRS.size() << " results\n";
        return mRS;
}
