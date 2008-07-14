/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2008 by Lukas Gamper <mistral@student.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>,
*                            Niall Moran <nmoran@thphys.nuim.ie>
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

/* $Id: archive_plot.cpp 2628 2007-11-20 02:20:42Z wistaria $ */

#include "archive_plot.hpp"

#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <iostream>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#define DEBUG

std::string Plot::strToLower(std::string inStr) {
        std::transform (inStr.begin(), inStr.end(), inStr.begin(), tolower);
        return inStr;
}

void Plot::writeFile(fs::path inOutFile, std::string inBuffer) {
        if (mVerbose)
                std::cout << "Writing Datafile '" << inOutFile.string() << "'";
        std::ofstream fout(inOutFile.native_directory_string().c_str());
    if (!fout.good())
            throw std::runtime_error(std::string("Could not open file: ") + inOutFile.string());
    fout << inBuffer;
        fout.close();
        if (mVerbose)
                std::cout << ": OK" << std::endl;
}

std::string operatorToString(std::string const& in) {
  if (in == "lessthan")
    return "<";
  else if (in == "lessorequalthan")
    return "<=";
  else if (in == "greaterthan")
    return ">";
  else if (in == "greaterorequalthan")
    return ">=";
  else if (in == "notequal")
    return "!=";
  else if (in == "equal")
    return "=";
  else
    throw std::runtime_error(std::string("Unknown Operator: ") + in);
  return "";
}

void append_where(std::string& where, std::string const& w) {
  if (where.size()) where += " AND ";
  where += w;
}

void append_from(std::string& from, std::string const& f) {
  if (from.size()) from += ", ";
  from += f;
}



void Plot::exec(Node inNode, std::string inInFile) {

  typedef std::list<std::map<std::string, std::string> > result_type;

  // initialize variables
  Node bodyNode = inNode.nodeTest("plot").front();

  int alias_number = 0;
  std::map<std::string, std::string> parameter_aliases;
  std::map<std::string, std::string> measurement_aliases;

  std::map<std::string, std::list<std::map<std::string, std::string> > > constraints;
  std::vector<std::string> constraintNames;
  std::list<Node> constrainNodes;

  std::string buffer = "";
  std::string bufferBody = "";
  int plot_num = 0 ;

  // find Node with axis tags and Check if all Attributes are there

  // parse <constraint> tags
  std::string constraintSQLFrom;
  std::string constraintSQLWhere;
  Node context = inNode.nodeTest("plot").front();
  while (context.nodeTest("for-each").size() > 0) context = context.nodeTest("for-each").front();
  BOOST_FOREACH(Node const& constraint, context.nodeTest("constraint")) {
    std::string name = constraint.getAttribute("name");
    std::string type = strToLower(constraint.getAttribute("type"));
    std::string value = constraint.getAttribute("value");
    std::string alias;
    if (type == "parameter") {
      if (parameter_aliases.find(name) == parameter_aliases.end()) {
        // assign a new alias name for parameter
        alias = std::string("x") + boost::lexical_cast<std::string>(alias_number++);
        append_from(constraintSQLFrom, "parameter AS " + alias);
        if (alias != "x0") append_where(constraintSQLWhere, "x0.fID=" + alias + ".fID");
        parameter_aliases[name] = alias;
      } else {
        alias = parameter_aliases[name];
      }
      // generate SQL string for constraint
      append_where(constraintSQLWhere, alias +".name='" + SQLite::quote(name) + "' AND " + alias
                   + ".value" + operatorToString(strToLower(constraint.getAttribute("operator")))
                   + value);
    } else {
      if (measurement_aliases.find(name) == measurement_aliases.end()) {
        // assign a new alias name for measurement
        alias = std::string("x") + boost::lexical_cast<std::string>(alias_number++);
        append_from(constraintSQLFrom, "measurement AS " + alias);
        if (alias != "x0") append_where(constraintSQLWhere, "x0.fID=" + alias + ".fID");
        measurement_aliases[name] = alias;
      } else {
        alias = measurement_aliases[name];
      }
      // generate SQL string for constraint
      append_where(constraintSQLWhere, alias +".name='" + SQLite::quote(name) + "' AND " + alias
                   + ".value" + operatorToString(constraint.getAttribute("operator")) + value);
    }
  }

  std::cout << "constraintSQLFrom: " << constraintSQLFrom << std::endl
            << "constraintSQLWhere: " << constraintSQLWhere << std::endl;

  // parse <for-each> tags
  std::string foreachSQLFrom;
  std::string foreachSQLWhere;
  std::vector<std::string> foreachName;
  std::vector<std::vector<std::string> > foreachData;
//   std::vector<std::vector<std::string> > forData;
//   std::vector<std::string> forVector;
//   std::vector<std::vector<std::string> > forMeta;
//   std::vector<unsigned int> forCount;
  context = inNode.nodeTest("plot").front();
  for (unsigned int pos = 0; context.nodeTest("for-each").size() > 0; ++pos) {
    context = context.nodeTest("for-each").front();
    std::string name = context.getAttribute("name");
    std::string from = constraintSQLFrom;
    std::string where = constraintSQLWhere;
    std::string alias;
    // forData.push_back(std::vector<std::string>());
    // forCount.push_back(0);
    // forMeta.push_back(std::vector<std::string>());
    // forVector.push_back(context.getAttribute("name"));
    if (parameter_aliases.find(name) == parameter_aliases.end()) {
      // assign a new alias name for parameter
      alias = std::string("x") + boost::lexical_cast<std::string>(alias_number++);
      append_from(foreachSQLFrom, "parameter AS " + alias);
      append_from(from, "parameter AS " + alias);
      append_where(where, alias + ".name='" + name + "'");
      if (alias != "x0") {
        append_where(foreachSQLWhere, "x0.fID=" + alias + ".fID");
        append_where(where, "x0.fID=" + alias + ".fID");
      }
      parameter_aliases[name] = alias;
    } else {
      alias = parameter_aliases[name];
    }
    std::string req("SELECT DISTINCT " + alias + ".value FROM " + from + " WHERE " + where + ";");
    result_type rs = mDB(req, true);
    foreachName.push_back(name);
    foreachData.push_back(std::vector<std::string>());
    for (result_type::iterator it = rs.begin(); it != rs.end(); it++)
      foreachData[pos].push_back((*it)["value"]);
  }

  std::cout << "foreachSQLFrom: " << foreachSQLFrom << std::endl
            << "foreachSQLWhere: " << foreachSQLWhere << std::endl;

  // define Order
  std::string orderSQLFrom;
  std::string orderSQLWhere;
  std::string orderSQLOrder;
  context = bodyNode.nodeTest("xaxis").front();
  std::string xname = context.getAttribute("type");
  std::string xtype = strToLower(context.getAttribute("type"));
  if (xtype == "parameter") {
    std::string alias;
    if (parameter_aliases.find(xname) == parameter_aliases.end()) {
      // assign a new alias name for parameter
      alias = std::string("x") + boost::lexical_cast<std::string>(alias_number++);
      append_from(orderSQLFrom, "parameter AS " + alias);
      if (alias != "x0") append_where(orderSQLWhere, "x0.fID=" + alias + ".fID");
      parameter_aliases[xname] = alias;
    } else {
      alias = parameter_aliases[xname];
    }
    orderSQLOrder = alias + ".value";
  } else {
    std::string alias;
    if (measurement_aliases.find(xname) == measurement_aliases.end()) {
      // assign a new alias name for measurement
      alias = std::string("x") + boost::lexical_cast<std::string>(alias_number++);
      append_from(orderSQLFrom, "measurement AS " + alias);
      if (alias != "x0") append_where(orderSQLWhere, "x0.fID=" + alias + ".fID");
      measurement_aliases[xname] = alias;
    } else {
      alias = measurement_aliases[xname];
    }
    if (xtype == "index") {
      orderSQLOrder = alias + ".indexvalue";
    } else {
      orderSQLOrder = alias + "." + xtype;
    }
  }

  std::cout << "orderSQLFrom: " << orderSQLFrom << std::endl
            << "orderSQLWhere: " << orderSQLWhere << std::endl
            << "orderSQLOrder: " << orderSQLOrder << std::endl;

//         preSQLFrom += std::string(preSQLFrom == "" ? "" : ", ") + (type == "parameter" ? "parameter" : "measurement") + " AS o";
//         preSQLWhere += std::string(preSQLWhere == "" ? "" : " AND ") + "x0.fID=o.fID AND o.name='" + SQLite::quote(bodyNode.nodeTest("xaxis").front().getAttribute("name")) + "'";
//         preSQLOrder = "o." + (type == "parameter" ? "value" : (type == "index" ? "indexvalue" : type));

//         std::cout << "preSQLFrom: " << preSQLFrom << std::endl
//                   << "preSQLWhere: " << preSQLWhere << std::endl
//                   << "preSQLOrder " << preSQLOrder << std::endl;

//         // find file ID's
//         if (forData.size() > 0) { //while there is at lease one for loop
//                 unsigned int table = 0;
//                 while(forCount[0] < forData[0].size()) {

//                         std::string where = preSQLWhere;
//                         std::string from = preSQLFrom;
//                         for (unsigned int pos = 0; pos < forData.size(); pos++) {
//                                 where += std::string(preSQLWhere == "" ? "" : " AND ") + "x0.fID=" + forNames[forVector[pos]] + ".fID AND "
//                                         + forNames[forVector[pos]] + ".name='" + SQLite::quote(forVector[pos]) + "' AND "
//                                         + forNames[forVector[pos]] + ".value='" + forData[pos][forCount[pos]] + "'";
//                                 forMeta[pos].push_back(forVector[pos] + "=" + forData[pos][forCount[pos]]);
//                         }
//                         forCount[forCount.size() - 1]++; //increment the last element in the list by one
//                         for (unsigned int pos = forCount.size() - 1; pos > 0 && forCount[pos] == forData[pos].size(); pos--) {
//                                 forCount[pos] = 0;
//                                 forCount[pos - 1]++;
//                         }
//                         std::string req = std::string("SELECT DISTINCT x0.fID AS id FROM ") + from + " "
//                           + "WHERE " + where + " ORDER BY " + preSQLOrder;
//                         result_type rs = mDB(req, true);
// #if 0
//                         if (mVerbose) std::clog << "Executed Query: " << req << std::endl
//                                                 << rs.size() << " files selected" << std::endl;
// #endif
//                         fIDs.push_back(std::list<std::string>());
//                         for (result_type::iterator it = rs.begin(); it != rs.end(); it++)
//                                 fIDs[table].push_back((*it)["id"]);
//                         table++;
//                         plot_num = table;
//                 }
//         } else {
//                 fIDs.push_back(std::list<std::string>());
//                 std::string req = std::string("SELECT x0.fID AS id FROM ") + preSQLFrom + " WHERE "
//                   + preSQLWhere + " ORDER BY " + preSQLOrder;
//                 result_type rs = mDB(req, true);
// #if 0
//                 if (mVerbose) std::clog << "Executed Query: " << req << std::endl
//                                         << rs.size() << " entries selected" << std::endl;
// #endif
//                 for (result_type::iterator it = rs.begin(); it != rs.end(); it++)
//                         fIDs[0].push_back((*it)["id"]);
//                 plot_num = 1 ;
//         }

//         // Parameters cannot habe indeces or errors
//         if ((bodyNode.nodeTest("xaxis").front().getAttribute("type") == "parameter" && (bodyNode.nodeTest("xaxis").front().getAttribute("error") != ""
//                         || bodyNode.nodeTest("xaxis").front().getAttribute("error") != ""))
//                         || (bodyNode.nodeTest("yaxis").front().getAttribute("type") == "parameter" && (bodyNode.nodeTest("yaxis").front().getAttribute("error") != ""
//                         || bodyNode.nodeTest("yaxis").front().getAttribute("error") != "")))
//                 throw std::runtime_error("Only measurements can habe errorbars and indeces");

//         // initialize variables
//         forCount.clear();
//         for        (std::vector<std::vector<std::string> >::iterator it = forMeta.begin(); it != forMeta.end(); it++)
//                 forCount.push_back(0);

//         // Create Header
//         if (inNode.nodeTest("plot").front().getAttribute("output") == "text")
//                 buffer += std::string("# ") + inNode.nodeTest("plot").front().getAttribute("name") + "\n";
//         else if (inNode.nodeTest("plot").front().getAttribute("output") == "html")
//                 buffer += std::string("<h1>") + inNode.nodeTest("plot").front().getAttribute("name") + "</h1>\n";
//         else if (inNode.nodeTest("plot").front().getAttribute("output") == "xmgr") {
//                 buffer += std::string("@g0 on\n@with g0\n@    legend on\n@    title \"" + inNode.nodeTest("plot").front().getAttribute("name") + "\""
//                         + "\n@    xaxis  label \"") + bodyNode.nodeTest("xaxis").front().getAttribute("name") + "\""
//                         + "\n@    yaxis  label \"" + bodyNode.nodeTest("yaxis").front().getAttribute("name") + "\"";
//         }

//         // Loop over For-Loops
//         unsigned int table = 0;
//         bool firstLoop = true;
//         while( ((forCount[0] < forMeta[0].size()) || (forMeta.size() == 0 && firstLoop) ) && table < plot_num ) {

//                 if (inNode.nodeTest("plot").front().getAttribute("output") == "text") {
//                         buffer += std::string("\n#");
//                         for (unsigned int pos = 0; pos < forMeta.size(); pos++)
//                                 buffer += std::string(" ") + forMeta[pos][forCount[pos]];
//                         buffer += std::string("\n# ") + bodyNode.nodeTest("xaxis").front().getAttribute("name") + (bodyNode.nodeTest("xaxis").front().getAttribute("error") != ""
//                                         ? "\terror(" + bodyNode.nodeTest("xaxis").front().getAttribute("name") + ")\t" : "\t")
//                                 + bodyNode.nodeTest("yaxis").front().getAttribute("name") + (bodyNode.nodeTest("yaxis").front().getAttribute("error") != ""
//                                         ? "\terror(" + bodyNode.nodeTest("yaxis").front().getAttribute("name") + ")" : "") + std::string("\n");
//                 } else if (inNode.nodeTest("plot").front().getAttribute("output") == "html") {
//                         buffer += std::string("\n<h2>");
//                         for (unsigned int pos = 0; pos < forMeta.size(); pos++)
//                                 buffer += std::string(" ") + forMeta[pos][forCount[pos]];
//                         buffer += std::string("</h2>\n<table><tr><td>") + bodyNode.nodeTest("xaxis").front().getAttribute("name")
//                                 + (bodyNode.nodeTest("xaxis").front().getAttribute("error") != ""
//                                         ? "</td><td>error(" + bodyNode.nodeTest("xaxis").front().getAttribute("name") + ")</td><td>" : "</td><td>")
//                                 + bodyNode.nodeTest("yaxis").front().getAttribute("name") + (bodyNode.nodeTest("yaxis").front().getAttribute("error") != ""
//                                         ? "</td><td>error(" + bodyNode.nodeTest("yaxis").front().getAttribute("name") + ")</td></tr>" : "</td></tr>");
//                 } else if (inNode.nodeTest("plot").front().getAttribute("output") == "xmgr") {
//                         buffer += std::string("\n@    s") + boost::lexical_cast<std::string>(table) + " type xy"
//                                 + (bodyNode.nodeTest("xaxis").front().getAttribute("error") != "" ? "dx" : "")
//                                 + (bodyNode.nodeTest("yaxis").front().getAttribute("error") != "" ? "dy" : "")
//                                 + "\n@    s" + boost::lexical_cast<std::string>(table) + " legend \"";
//                                                    for (unsigned int pos = 0; pos < forMeta.size(); pos++)
//                         buffer += std::string(" ") + forMeta[pos][forCount[pos]];
//                                               buffer += std::string("\"");
//                                          bufferBody += std::string("\n@target G0.S") + boost::lexical_cast<std::string>(table) + "\n@type xy"
//                                         + (bodyNode.nodeTest("xaxis").front().getAttribute("error") != "" ? "dx" : "")
//                                 + (bodyNode.nodeTest("yaxis").front().getAttribute("error") != "" ? "dy" : "") + "\n";
//                 } //end if
//                 /*forCount[forCount.size() - 1]++;
//                 for (unsigned int pos = forCount.size() - 1; pos > 0 && forCount[pos] >= forMeta[pos].size() ; pos--) {
//                         forCount[pos] = 0;
//                         forCount[pos - 1]++;
//                 }*/

//                 for ( unsigned int pos = 0 ; pos < forCount.size() ; pos++ ){
//                         forCount[pos]++;
//                 }

//                                          //forCount[forCount.size() - 1]++;

//                 for (std::list<std::string>::iterator it = fIDs[table].begin(); it != fIDs[table].end(); it++) {
//                         std::string req = std::string("SELECT u.") + (bodyNode.nodeTest("xaxis").front().getAttribute("type") == "parameter"
//                                                 ? "value" : (bodyNode.nodeTest("xaxis").front().getAttribute("type") == "index"
//                                                 ? "indexvalue" : bodyNode.nodeTest("xaxis").front().getAttribute("type"))) + " AS x, "
//                                         + (bodyNode.nodeTest("xaxis").front().getAttribute("error") != "" ? "u.error AS ex, " : "")
//                                         + "v." + (bodyNode.nodeTest("yaxis").front().getAttribute("type") == "parameter"
//                                                 ? "value" : (bodyNode.nodeTest("yaxis").front().getAttribute("type") == "index"
//                                                 ? "indexvalue" : bodyNode.nodeTest("yaxis").front().getAttribute("type"))) + " AS y "
//                                         + (bodyNode.nodeTest("yaxis").front().getAttribute("error") != "" ? ", v.error AS ey " : "")
//                                 + "FROM "
//                                         + (bodyNode.nodeTest("xaxis").front().getAttribute("type") == "parameter" ? "parameter" : "measurement") + " AS u, "
//                                         + (bodyNode.nodeTest("yaxis").front().getAttribute("type") == "parameter" ? "parameter" : "measurement") + " AS v "
//                                 + "WHERE u.fID='" + (*it) + "' AND v.FID='" + (*it) + "'"
//                                         + " AND u.name='" + SQLite::quote(bodyNode.nodeTest("xaxis").front().getAttribute("name")) + "'"
//                                         + " AND v.name='" + SQLite::quote(bodyNode.nodeTest("yaxis").front().getAttribute("name")) + "'"
//                                         + (bodyNode.nodeTest("xaxis").front().getAttribute("index") != "" ? " AND u.indexvalue='"
//                                                 + bodyNode.nodeTest("xaxis").front().getAttribute("index") + "'" : "")
//                                         + (bodyNode.nodeTest("yaxis").front().getAttribute("index") != "" ? " AND v.indexvalue='"
//                                            + bodyNode.nodeTest("yaxis").front().getAttribute("index") + "'" : "");
//                         result_type rs = mDB(req);
// #if 0
//                         if (mVerbose) std::clog << "Executed Query: " << req << std::endl
//                                                 << rs.size() << " entries selected\n";
// #endif
//                         if (rs.size() > 0) {
//                                 if (rs.size() != 1)
//                                         throw std::runtime_error("The Result ob the Query is not unique!");
//                                 if (inNode.nodeTest("plot").front().getAttribute("output") == "text")
//                                         buffer += std::string(rs.front()["x"]) + (bodyNode.nodeTest("xaxis").front().getAttribute("error") != "" ? "\t" + rs.front()["ex"] + "\t" : "\t")
//                                                 + rs.front()["y"] + (bodyNode.nodeTest("yaxis").front().getAttribute("error") != "" ? "\t" + rs.front()["ey"] : "") + "\n";
//                                 else if (inNode.nodeTest("plot").front().getAttribute("output") == "html")
//                                         buffer += std::string("<tr><td allign=\"right\">") + rs.front()["x"] + (bodyNode.nodeTest("xaxis").front().getAttribute("error") != "" ? "</td><td>" + rs.front()["ex"] + "</td><td allign=\"right\">" : "</td><td allign=\"right\">")
//                                                 + rs.front()["y"] + (bodyNode.nodeTest("yaxis").front().getAttribute("error") != "" ? "</td><td allign=\"right\">" + rs.front()["ey"] : "") + "</td></tr>\n";
//                                 else if (inNode.nodeTest("plot").front().getAttribute("output") == "xmgr") {
//                                         bufferBody += std::string(rs.front()["x"])
//                                                 + (bodyNode.nodeTest("xaxis").front().getAttribute("error") != "" ? " " + rs.front()["ex"] + " " : " ")
//                                                 + rs.front()["y"]
//                                                 + (bodyNode.nodeTest("yaxis").front().getAttribute("error") != "" ? " " + rs.front()["ey"] : "") + " \n";
//                                 }
//                         }
//                 }
//                 if (inNode.nodeTest("plot").front().getAttribute("output") == "html")
//                         buffer += "</table>\n";
//                 else if (inNode.nodeTest("plot").front().getAttribute("output") == "xmgr")
//                         bufferBody += "&";
//                 table++;
//                 firstLoop = false;
//         }
//         if (inNode.nodeTest("plot").front().getAttribute("output") == "xmgr")
//                 buffer += bufferBody;

//         // Datei schreiben
//         writeFile(mOutPath / (inInFile + "." + inNode.nodeTest("plot").front().getAttribute("output")), buffer);
}
