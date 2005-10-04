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

#include "archive_plot.hpp"

#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <iostream>

#include <boost/lexical_cast.hpp>

std::string Plot::strToLower(std::string inStr) {
        transform (inStr.begin(), inStr.end(), inStr.begin(), tolower);
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

void Plot::exec(Node inNode, std::string inInFile) {

        // initialize variables
        std::vector<std::vector<std::string> > forData = std::vector<std::vector<std::string> >();
        std::vector<std::string> forVector = std::vector<std::string>();
        std::vector<std::vector<std::string> > forMeta = std::vector<std::vector<std::string> >();
        std::vector<unsigned int> forCount = std::vector<unsigned int>();
        std::map<std::string, std::string> forNames = std::map<std::string, std::string>();
        std::list<std::map<std::string, std::string> > rs = std::list<std::map<std::string, std::string> >();
        std::map<std::string, std::list<std::map<std::string, std::string> > > constraints = std::map<std::string, std::list<std::map<std::string, std::string> > >();
        std::vector<std::string> constraintNames = std::vector<std::string>();
        std::list<Node> constrainNodes;
        std::vector<std::list<std::string> > fIDs = std::vector<std::list<std::string> >();
        Node bodyNode = inNode.nodeTest("plot").front();
        int distinguish = 0;
        std::string preSQLWhere = "";
        std::string preSQLFrom = "";
        std::string preSQLOrder;
        std::string buffer = "";
        std::string bufferBody = "";

        // find Node with axis tags and Check if all Attributes are there
        while (bodyNode.nodeTest("for-each").size() > 0)
                bodyNode = bodyNode.nodeTest("for-each").front();

        // parse Constraints
        constrainNodes = bodyNode.nodeTest("constraint");
        for (std::list<Node>::iterator it = constrainNodes.begin(); it != constrainNodes.end(); it++) {
                std::map<std::string, std::string> constraint = std::map<std::string, std::string>();
                if (constraints.find(it->getAttribute("name")) == constraints.end()) {
                        constraints[it->getAttribute("name")] = std::list<std::map<std::string, std::string> >();
                        constraintNames.push_back(it->getAttribute("name"));
                }
                constraint["name"] = it->getAttribute("name");
                constraint["type"] = it->getAttribute("type");
                constraint["value"] = it->getAttribute("value");
                if (strToLower(it->getAttribute("operator")) == "lessthan") 
                        constraint["operator"] = "<";
                else if (strToLower(it->getAttribute("operator")) == "lessorequalthan") 
                        constraint["operator"] = "<=";
                else if (strToLower(it->getAttribute("operator")) == "greaterthan") 
                        constraint["operator"] = ">";
                else if (strToLower(it->getAttribute("operator")) == "greaterorequalthan") 
                        constraint["operator"] = ">=";
                else if (strToLower(it->getAttribute("operator")) == "notequal") 
                        constraint["operator"] = "!=";
                else if (strToLower(it->getAttribute("operator")) == "equal") 
                        constraint["operator"] = "=";
                else 
                        throw std::runtime_error(std::string("Unbekannter Operator: ") + it->getAttribute("operator"));
                constraints[it->getAttribute("name")].push_back(constraint);
        }

        // constrains to sql
        for        (std::vector<std::string>::iterator it = constraintNames.begin(); it != constraintNames.end(); it++) {
                bool create = false;
                std::string type = "";
                for        (std::list<std::map<std::string, std::string> >::iterator subIT = constraints[*it].begin(); subIT != constraints[*it].end(); subIT++) {
                        if (type == "")
                                type = (*subIT)["type"];
                        if ((type == "parameter" || (*subIT)["type"] == "parameter") && type != (*subIT)["type"])
                                throw std::runtime_error(std::string("The type of the variable '") + (*subIT)["name"] + "' is not unique!");
                        if ((*subIT)["type"] != "index") {
                                preSQLWhere += std::string(preSQLWhere == "" ? "x" : " AND x") + boost::lexical_cast<std::string>(distinguish) + "." 
                                        + ((*subIT)["type"] == "parameter" ? "value" : (*subIT)["type"])
                                        + (*subIT)["operator"] + "'" + (*subIT)["value"] + "'";
                                create = true;
                        }
                }
                if (create) {
                        forNames[*it] = std::string("x") + boost::lexical_cast<std::string>(distinguish);
                        preSQLFrom += std::string(preSQLFrom == "" ? "" : ", ") 
                                + (constraints[*it].front()["type"] == "parameter" ? "parameter" : "mesurement") + " AS x" + boost::lexical_cast<std::string>(distinguish);
                        preSQLWhere += std::string(preSQLWhere == "" ? "" : " AND x") + boost::lexical_cast<std::string>(distinguish) + ".name='" + *it + "'";
                        if (distinguish > 0)
                                preSQLWhere += std::string(preSQLWhere == "" ? "" : " AND ") + "x0.fID=x" + boost::lexical_cast<std::string>(distinguish) + ".fID";
                        distinguish++;
                }
        }

        // define Order
        std::string type = strToLower(bodyNode.nodeTest("xaxis").front().getAttribute("type"));
        preSQLFrom += std::string(preSQLFrom == "" ? "" : ", ") + (type == "parameter" ? "parameter" : "mesurement") + " AS o";
        preSQLWhere += std::string(preSQLWhere == "" ? "" : " AND ") + "x0.fID=o.fID AND o.name='" + bodyNode.nodeTest("xaxis").front().getAttribute("name") + "'";
        preSQLOrder = "o." + (type == "parameter" ? "value" : (type == "index" ? "indexvalue" : type));

        // parse For-Tags
        Node context = inNode.nodeTest("plot").front();
        for (unsigned int pos = 0; context.nodeTest("for-each").size() > 0; pos++) {
                context = context.nodeTest("for-each").front();
                forData.push_back(std::vector<std::string>());
                forCount.push_back(0);
                forMeta.push_back(std::vector<std::string>());
                forVector.push_back(context.getAttribute("name"));
                if (forNames.find(context.getAttribute("name")) == forNames.end()) {
                        forNames[context.getAttribute("name")] = std::string("x") + boost::lexical_cast<std::string>(distinguish);
                        preSQLFrom += std::string(preSQLFrom == "" ? "" : ", ") + "parameter AS x" + boost::lexical_cast<std::string>(distinguish);
                        distinguish++;
                }
                std::string where = "";
                if (constraints.find(context.getAttribute("name")) != constraints.end()) {
                        for        (std::list<std::map<std::string, std::string> >::iterator it = constraints[context.getAttribute("name")].begin(); 
                                        it != constraints[context.getAttribute("name")].end(); it++)
                                where += std::string(" AND value") + (*it)["operator"] + "'" + (*it)["value"] + "'";
                        constraints.erase(context.getAttribute("name"));
                }
                rs = mDB(std::string("SELECT DISTINCT value FROM parameter WHERE name='") + context.getAttribute("name")  + "'" + where + ";", true);
                for (std::list<std::map<std::string, std::string> >::iterator it = rs.begin(); it != rs.end(); it++)
                        forData[pos].push_back((*it)["value"]);
        }
        
        // find file ID's
        if (forData.size() > 0) {
                unsigned int table = 0;
                while(forCount[0] < forData[0].size()) {
                        
                        std::string where = preSQLWhere;
                        std::string from = preSQLFrom;
                        for (unsigned int pos = 0; pos < forData.size(); pos++) {
                                where += std::string(preSQLWhere == "" ? "" : " AND ") + "x0.fID=" + forNames[forVector[pos]] + ".fID AND " 
                                        + forNames[forVector[pos]] + ".name='" + forVector[pos] + "' AND "
                                        + forNames[forVector[pos]] + ".value='" + forData[pos][forCount[pos]] + "'";
                                forMeta[pos].push_back(forVector[pos] + "=" + forData[pos][forCount[pos]]);
                        }
                        forCount[forCount.size() - 1]++;
                        for (unsigned int pos = forCount.size() - 1; pos > 0 && forCount[pos] == forData[pos].size(); pos--) {
                                forCount[pos] = 0;
                                forCount[pos - 1]++;
                        }
                        #ifdef DEBUG
                                std::cout << "Execute Query: " << std::string("SELECT DISTINCT x0.fID AS id\n\tFROM ") + from
                                        + "\n\tWHERE " + where + "\n\tORDER BY " + preSQLOrder << std::endl;
                        #endif
                                rs = mDB(std::string("SELECT DISTINCT x0.fID AS id FROM ") + from + " "
                                        + "WHERE " + where + " ORDER BY " + preSQLOrder, true);
                        #ifdef DEBUG
                                std::cout << rs.size() << " Files selected" << std::endl;
                        #endif
                        fIDs.push_back(std::list<std::string>());
                        for (std::list<std::map<std::string, std::string> >::iterator it = rs.begin(); it != rs.end(); it++) 
                                fIDs[table].push_back((*it)["id"]);
                        table++;
                }
        } else {
                fIDs.push_back(std::list<std::string>());
                rs = mDB(std::string("SELECT x0.fID AS id FROM ") + preSQLFrom + " WHERE " 
                        + preSQLWhere + " ORDER BY " + preSQLOrder, true);
                for (std::list<std::map<std::string, std::string> >::iterator it = rs.begin(); it != rs.end(); it++)
                        fIDs[0].push_back((*it)["id"]);
        }

        // Parameters cannot habe indeces or errors
        if ((bodyNode.nodeTest("xaxis").front().getAttribute("type") == "parameter" && (bodyNode.nodeTest("xaxis").front().getAttribute("error") != "" 
                        || bodyNode.nodeTest("xaxis").front().getAttribute("error") != ""))
                        || (bodyNode.nodeTest("yaxis").front().getAttribute("type") == "parameter" && (bodyNode.nodeTest("yaxis").front().getAttribute("error") != "" 
                        || bodyNode.nodeTest("yaxis").front().getAttribute("error") != "")))
                throw std::runtime_error("Only measurements can habe errorbars and indeces");

        // initialize variables
        rs = std::list<std::map<std::string, std::string> >();
        forCount = std::vector<unsigned int>();
        for        (std::vector<std::vector<std::string> >::iterator it = forMeta.begin(); it != forMeta.end(); it++)
                forCount.push_back(0);

        // Create Header                
        if (inNode.nodeTest("plot").front().getAttribute("output") == "text")                                        
                buffer += std::string("# ") + inNode.nodeTest("plot").front().getAttribute("name") + "\n";
        else if (inNode.nodeTest("plot").front().getAttribute("output") == "html")
                buffer += std::string("<h1>") + inNode.nodeTest("plot").front().getAttribute("name") + "</h1>\n";
        else if (inNode.nodeTest("plot").front().getAttribute("output") == "xmgr") {
                buffer += std::string("@g0 on\n@with g0\n@    legend on\n@    title \"" + inNode.nodeTest("plot").front().getAttribute("name") + "\""
                        + "\n@    xaxis  label \"") + bodyNode.nodeTest("xaxis").front().getAttribute("name") + "\""
                        + "\n@    yaxis  label \"" + bodyNode.nodeTest("yaxis").front().getAttribute("name") + "\"";
        }
        
        // Loop over For-Loops
        unsigned int table = 0;
        bool firstLoop = true;
        while((forCount[0] < forMeta[0].size()) || (forMeta.size() == 0 && firstLoop)) {

                if (inNode.nodeTest("plot").front().getAttribute("output") == "text") {                        
                        buffer += std::string("\n#");
                        for (unsigned int pos = 0; pos < forMeta.size(); pos++)
                                buffer += std::string(" ") + forMeta[pos][forCount[pos]];
                        buffer += std::string("\n#") + bodyNode.nodeTest("xaxis").front().getAttribute("name") + (bodyNode.nodeTest("xaxis").front().getAttribute("error") != "" 
                                        ? "\terror(" + bodyNode.nodeTest("xaxis").front().getAttribute("name") + ")\t" : "\t")
                                + bodyNode.nodeTest("yaxis").front().getAttribute("name") + (bodyNode.nodeTest("yaxis").front().getAttribute("error") != "" 
                                        ? "\terror(" + bodyNode.nodeTest("yaxis").front().getAttribute("name") + ")" : "") + std::string("\n");
                } else if (inNode.nodeTest("plot").front().getAttribute("output") == "html") {
                        buffer += std::string("\n<h2>");
                        for (unsigned int pos = 0; pos < forMeta.size(); pos++)
                                buffer += std::string(" ") + forMeta[pos][forCount[pos]];
                        buffer += std::string("</h2>\n<table><tr><td>") + bodyNode.nodeTest("xaxis").front().getAttribute("name")
                                + (bodyNode.nodeTest("xaxis").front().getAttribute("error") != "" 
                                        ? "</td><td>error(" + bodyNode.nodeTest("xaxis").front().getAttribute("name") + ")</td><td>" : "</td><td>")
                                + bodyNode.nodeTest("yaxis").front().getAttribute("name") + (bodyNode.nodeTest("yaxis").front().getAttribute("error") != "" 
                                        ? "</td><td>error(" + bodyNode.nodeTest("yaxis").front().getAttribute("name") + ")</td></tr>" : "</td></tr>");
                } else if (inNode.nodeTest("plot").front().getAttribute("output") == "xmgr") {
                        buffer += std::string("\n@    s") + boost::lexical_cast<std::string>(table) + " type xy"
                                + (bodyNode.nodeTest("xaxis").front().getAttribute("error") != "" ? "dx" : "")
                                + (bodyNode.nodeTest("yaxis").front().getAttribute("error") != "" ? "dy" : "")
                                + "\n@    s" + boost::lexical_cast<std::string>(table) + " legend \"";
      for (unsigned int pos = 0; pos < forMeta.size(); pos++)
                          buffer += std::string(" ") + forMeta[pos][forCount[pos]];                                
            buffer += std::string("\"");                                
                        bufferBody += std::string("\n@target G0.S") + boost::lexical_cast<std::string>(table) + "\n@type xy"
                                + (bodyNode.nodeTest("xaxis").front().getAttribute("error") != "" ? "dx" : "")
                                + (bodyNode.nodeTest("yaxis").front().getAttribute("error") != "" ? "dy" : "") + "\n";
                }
                forCount[forCount.size() - 1]++;
                for (unsigned int pos = forCount.size() - 1; pos > 0 && forCount[pos] == forMeta[pos].size(); pos--) {
                        forCount[pos] = 0;
                        forCount[pos - 1]++;
                }
                for (std::list<std::string>::iterator it = fIDs[table].begin(); it != fIDs[table].end(); it++) {
                        rs = mDB(std::string("SELECT u.") + (bodyNode.nodeTest("xaxis").front().getAttribute("type") == "parameter"
                                                ? "value" : (bodyNode.nodeTest("xaxis").front().getAttribute("type") == "index" 
                                                ? "indexvalue" : bodyNode.nodeTest("xaxis").front().getAttribute("type"))) + " AS x, "
                                        + (bodyNode.nodeTest("xaxis").front().getAttribute("error") != "" ? "u.error AS ex, " : "")
                                        + "v." + (bodyNode.nodeTest("yaxis").front().getAttribute("type") == "parameter"  
                                                ? "value" : (bodyNode.nodeTest("yaxis").front().getAttribute("type") == "index" 
                                                ? "indexvalue" : bodyNode.nodeTest("yaxis").front().getAttribute("type"))) + " AS y "
                                        + (bodyNode.nodeTest("yaxis").front().getAttribute("error") != "" ? ", v.error AS ey " : "")
                                + "FROM "
                                        + (bodyNode.nodeTest("xaxis").front().getAttribute("type") == "parameter" ? "parameter" : "mesurement") + " AS u, "
                                        + (bodyNode.nodeTest("yaxis").front().getAttribute("type") == "parameter" ? "parameter" : "mesurement") + " AS v "
                                + "WHERE u.fID='" + (*it) + "' AND v.FID='" + (*it) + "'"
                                        + " AND u.name='" + bodyNode.nodeTest("xaxis").front().getAttribute("name") + "'"
                                        + " AND v.name='" + bodyNode.nodeTest("yaxis").front().getAttribute("name") + "'"
                                        + (bodyNode.nodeTest("xaxis").front().getAttribute("index") != "" ? " AND u.indexvalue='" 
                                                + bodyNode.nodeTest("xaxis").front().getAttribute("index") + "'" : "")
                                        + (bodyNode.nodeTest("yaxis").front().getAttribute("index") != "" ? " AND v.indexvalue='" 
                                                + bodyNode.nodeTest("yaxis").front().getAttribute("index") + "'" : "")
                                );
                        if (rs.size() > 0) {
                                if (rs.size() != 1)
                                        throw std::runtime_error("The Result ob the Query is not unique!"); 
                                if (inNode.nodeTest("plot").front().getAttribute("output") == "text")
                                        buffer += std::string(rs.front()["x"]) + (bodyNode.nodeTest("xaxis").front().getAttribute("error") != "" ? "\t" + rs.front()["ex"] + "\t" : "\t")
                                                + rs.front()["y"] + (bodyNode.nodeTest("yaxis").front().getAttribute("error") != "" ? "\t" + rs.front()["ey"] : "") + "\n";
                                else if (inNode.nodeTest("plot").front().getAttribute("output") == "html")
                                        buffer += std::string("<tr><td allign=\"right\">") + rs.front()["x"] + (bodyNode.nodeTest("xaxis").front().getAttribute("error") != "" ? "</td><td>" + rs.front()["ex"] + "</td><td allign=\"right\">" : "</td><td allign=\"right\">")
                                                + rs.front()["y"] + (bodyNode.nodeTest("yaxis").front().getAttribute("error") != "" ? "</td><td allign=\"right\">" + rs.front()["ey"] : "") + "</td></tr>\n";
                                else if (inNode.nodeTest("plot").front().getAttribute("output") == "xmgr") {
                                        bufferBody += std::string(rs.front()["x"])
                                                + (bodyNode.nodeTest("xaxis").front().getAttribute("error") != "" ? " " + rs.front()["ex"] + " " : " ")
                                                + rs.front()["y"]                                                
                                                + (bodyNode.nodeTest("yaxis").front().getAttribute("error") != "" ? " " + rs.front()["ey"] : "") + " \n";
                                }
                        }
                }
                if (inNode.nodeTest("plot").front().getAttribute("output") == "html")
                        buffer += "</table>\n";
                else if (inNode.nodeTest("plot").front().getAttribute("output") == "xmgr")
                        bufferBody += "&";
                table++;
                firstLoop = false;
        }
        if (inNode.nodeTest("plot").front().getAttribute("output") == "xmgr")
                buffer += bufferBody;
        
        // Datei schreiben
        writeFile(mOutPath / (inInFile + "." + inNode.nodeTest("plot").front().getAttribute("output")), buffer);
}
