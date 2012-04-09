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

/* $Id$ */

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
        std::ofstream fout(inOutFile.string().c_str());
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
  if (w.size()) {
    if (where.size()) where += " AND ";
    where += w;
  }
}

void append_from(std::string& from, std::string const& f) {
  if (f.size()) {
    if (from.size()) from += ", ";
    from += f;
  }
}

void check_duplicate(std::list<std::map<std::string, std::string> >& rs) {
  std::string old = "";
  for (std::list<std::map<std::string, std::string> >::iterator itr = rs.begin();
       itr != rs.end(); ++itr) {
    if (old != "" && (*itr)["x"] == old)
      boost::throw_exception(std::runtime_error("duplicated entries"));
    old = (*itr)["x"];
  }
}

void create_header(std::string& buffer, std::string const& output_type,
                   std::string const& plot_name, std::string const& xlabel,
                   std::string const& ylabel) {
  if (output_type == "text")
    buffer += "# " + plot_name + "\n";
  else if (output_type == "html")
    buffer += "<h1>" + plot_name + "</h1>\n";
  else if (output_type == "gnuplot")
    buffer += "set title \"" + plot_name + "\"\n"
      + "set xlabel \"" + xlabel + "\"\n"
      + "set ylabel \"" + ylabel + "\"\n"
      + "plot ";
  else if (output_type == "xmgr") {
    buffer += "@g0 on\n@with g0\n@    legend on\n@    title \"" + plot_name + "\""
      "\n@    xaxis  label \"" + xlabel + "\""
      "\n@    yaxis  label \"" + ylabel + "\"";
  } else if (output_type == "xml") {
    buffer += "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
      "<!DOCTYPE plot SYSTEM \"conf/plotdata.dtd\">\n"
      "<plotdata";
    if (plot_name.size()) buffer += " name=\"" + plot_name + "\"";
    buffer += ">\n"
      "  <xaxis label=\"" + xlabel + "\"/>\n"
      "  <yaxis label=\"" + ylabel + "\"/>\n";
  }
}

void create_trailer(std::string& buffer, std::string const& output_type) {
  if (output_type == "xml") {
    buffer += "</plotdata>\n";
  }
}

void write_result(std::string& buffer, std::string& bufferBody,
                  std::list<std::map<std::string, std::string> >& rs,
                  std::string const& output_type, int table,
                  std::string const& xlabel, std::string const& ylabel, bool xerror, bool yerror,
                  std::vector<std::string> const& foreachName,
                  std::vector<std::vector<std::string> > const& foreachData,
                  std::vector<unsigned int> const& foreachOffset) {
  std::string label;
  int p = table;
  for (std::size_t v = 0; v < foreachData.size(); ++v, p %= foreachOffset[v-1]) {
    if (label.size()) label += " ";
    label += foreachName[v] + "=" + foreachData[v][p / foreachOffset[v]];
  }
  if (output_type == "text") {
    buffer += "\n# " + label;
    buffer += "\n# " + xlabel + (xerror ? "\terror(" + xlabel + ")" : "")
      + "\t" + ylabel + (yerror ? "\terror(" + ylabel + ")" : "") + "\n";
    for (std::list<std::map<std::string, std::string> >::iterator itr = rs.begin();
         itr != rs.end(); ++itr) {
      buffer += std::string((*itr)["x"]) + (xerror ? "\t" + (*itr)["ex"] + "\t" : "\t")
        + (*itr)["y"] + (yerror ? "\t" + (*itr)["ey"] : "") + "\n";
    }
  } else if (output_type == "html") {
    buffer += "\n<h2>" + label + "</h2>\n";
    buffer += "<table><tr><td>" + xlabel + "</td>"
      + (xerror ? "<td>error(" + xlabel + ")</td>" : "")
      + "<td>" + ylabel + "</td>"
      + (yerror ? "<td>error(" + ylabel + ")</td>" : "")
      + "</tr>\n";
    for (std::list<std::map<std::string, std::string> >::iterator itr = rs.begin();
         itr != rs.end(); ++itr) {
      buffer += "<tr><td allign=\"right\">" + (*itr)["x"] + "</td>"
        + (xerror ? "<td allign=\"right\">" + (*itr)["ex"] + "</td>" : "")
        + "<tr><td allign=\"right\">" + (*itr)["y"] + "</td>"
        + (yerror ? "<td allign=\"right\">" + (*itr)["ey"] + "</td>" : "")
        + "</tr>\n";
    }
    buffer += "</table>\n";
  } else if (output_type == "gnuplot") {
    if (table != 0) buffer += ", ";
    if (xerror)
      if (yerror)
        buffer += "'-'  using 1:2:3:4 title \"" + label + "\" with xyerrorlines";
      else
        buffer += "'-'  using 1:2:3 title \"" + label + "\" with xerrorlines";
    else
      if (yerror)
        buffer += "'-'  using 1:2:3 title \"" + label + "\" with yerrorlines";
      else
        buffer += "'-'  using 1:2 title \"" + label + "\" with lines";
    if (table == 0) bufferBody += "\n";
    for (std::list<std::map<std::string, std::string> >::iterator itr = rs.begin();
         itr != rs.end(); ++itr)
      bufferBody += (*itr)["x"] + " " + (xerror ? (*itr)["ex"] + " " : "")
        + (*itr)["y"] + " " + (yerror ? (*itr)["ey"] + " " : "") + "\n";
    bufferBody += "e\n";
  } else if (output_type == "xmgr") {
    buffer += "\n@    s" + boost::lexical_cast<std::string>(table)
      + " type xy" + (xerror ? "dx" : "") + (yerror ? "dy" : "")
      + "\n@    s" + boost::lexical_cast<std::string>(table) + " legend \"" + label + "\"";
    bufferBody += std::string("\n@target G0.S") + boost::lexical_cast<std::string>(table)
      + "\n@type xy" +  (xerror ? "dx" : "") + (yerror ? "dy" : "") + "\n";
    for (std::list<std::map<std::string, std::string> >::iterator itr = rs.begin();
         itr != rs.end(); ++itr) {
      bufferBody += (*itr)["x"] + " " + (xerror ? (*itr)["ex"] + " " : "")
        + (*itr)["y"] + " " + (yerror ? (*itr)["ey"] + " " : "") + "\n";
      bufferBody += "&";
    }
  } else if (output_type == "xml") {
    buffer += "  <set label=\"" + label + "\">\n";
    for (std::list<std::map<std::string, std::string> >::iterator itr = rs.begin();
         itr != rs.end(); ++itr) {
      buffer += "    <point>"
        "<x>" + (*itr)["x"] + "</x>"
        + (xerror ? "<dx>" + (*itr)["ex"] + "</dx>" : "")
        + "<y>" + (*itr)["y"] + "</y>"
        + (yerror ? "<dy>" + (*itr)["ey"] + "</dy>" : "")
        + "</point>\n";
    }
    buffer += "  </set>\n";
  }
}

void Plot::exec(Node inNode, std::string inInFile) {

  typedef std::list<std::map<std::string, std::string> > result_type;

  std::string plot_name = inNode.nodeTest("plot").front().getAttribute("name");
  std::string output_type = inNode.nodeTest("plot").front().getAttribute("output");

  int alias_number = 0;
  std::map<std::string, std::string> parameter_aliases;
  std::map<std::string, std::string> measurement_aliases;

  // std::map<std::string, std::list<std::map<std::string, std::string> > > constraints;
  // std::vector<std::string> constraintNames;
  // std::list<Node> constrainNodes;

  std::string buffer;
  std::string bufferBody;

  // parse <constraint> tags
  std::string constraintSQLFrom;
  std::string constraintSQLWhere;
  Node context = inNode.nodeTest("plot").front();
  while (context.nodeTest("for-each").size() > 0) context = context.nodeTest("for-each").front();
  // BOOST_FOREACH(Node const& constraint, context.nodeTest("constraint")) {
  std::list<Node> constraints = context.nodeTest("constraint");
  for (std::list<Node>::const_iterator itr = constraints.begin(); itr != constraints.end(); ++itr) {
    const Node& constraint = *itr;
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
                   + "'" + value + "'");
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
  // std::cout << "constraintSQLFrom: " << constraintSQLFrom << std::endl
  //           << "constraintSQLWhere: " << constraintSQLWhere << std::endl;

  // parse <for-each> tags
  std::string from = constraintSQLFrom;
  std::string where = constraintSQLWhere;
  std::string foreachSQLFrom;
  std::string foreachSQLWhere;
  std::vector<std::string> foreachName;
  std::vector<std::vector<std::string> > foreachData;
  context = inNode.nodeTest("plot").front();
  for (unsigned int pos = 0; context.nodeTest("for-each").size() > 0; ++pos) {
    context = context.nodeTest("for-each").front();
    std::string name = context.getAttribute("name");
    std::string alias;
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
    std::string req("SELECT DISTINCT " + alias + ".value as value FROM " + from + " WHERE " + where + ";");
    result_type rs = mDB(req, true);
    foreachName.push_back(name);
    foreachData.push_back(std::vector<std::string>());
    for (result_type::iterator it = rs.begin(); it != rs.end(); it++)
      foreachData[pos].push_back((*it)["value"]);
  }

  int num_plots = 1;
  std::vector<unsigned int> foreachOffset(foreachData.size());
  for (int v = foreachData.size() - 1; v >= 0; --v) {
    foreachOffset[v] = num_plots;
    num_plots *= foreachData[v].size();
  }

  // std::cout << "foreachSQLFrom: " << foreachSQLFrom << std::endl
  //           << "foreachSQLWhere: " << foreachSQLWhere << std::endl;

  // xaxis and yaxis
  std::string axisSQLFrom;
  std::string axisSQLWhere;
  std::string axisSQLOrder;
  context = inNode.nodeTest("plot").front();
  while (context.nodeTest("for-each").size() > 0) context = context.nodeTest("for-each").front();
  context = context.nodeTest("xaxis").front();
  std::string xname = context.getAttribute("name");
  std::string xtype = strToLower(context.getAttribute("type"));
  std::string xlabel = xname;
  if (xtype != "parameter" && xtype != "mean") xlabel = xtype + "(" + xlabel + ")";
  if (xtype == "index") xtype = "indexvalue";
  std::string xalias;
  bool xerror = (context.getAttribute("error") != "");
  if (xtype == "parameter") {
    if (parameter_aliases.find(xname) == parameter_aliases.end()) {
      // assign a new alias name for parameter
      xalias = std::string("x") + boost::lexical_cast<std::string>(alias_number++);
      append_from(axisSQLFrom, "parameter AS " + xalias);
      if (xalias != "x0") append_where(axisSQLWhere, "x0.fID=" + xalias + ".fID");
      parameter_aliases[xname] = xalias;
    } else {
      xalias = parameter_aliases[xname];
    }
    axisSQLOrder = xalias + ".value";
  } else {
    if (measurement_aliases.find(xname) == measurement_aliases.end()) {
      // assign a new alias name for measurement
      xalias = std::string("x") + boost::lexical_cast<std::string>(alias_number++);
      append_from(axisSQLFrom, "measurement AS " + xalias);
      if (xalias != "x0") append_where(axisSQLWhere, "x0.fID=" + xalias + ".fID");
      measurement_aliases[xname] = xalias;
    } else {
      xalias = measurement_aliases[xname];
    }
    if (context.getAttribute("index") != "") {
      append_where(axisSQLWhere, xalias + ".indexvalue=\""
                   + SQLite::quote(context.getAttribute("index")) + "\"");
      xlabel += " " + context.getAttribute("index");
    }
    axisSQLOrder = xalias + "." + xtype;
  }
  append_where(axisSQLWhere, xalias + ".name='" + SQLite::quote(xname) + "'");
  if (xtype == "parameter") xtype = "value";

  context = inNode.nodeTest("plot").front();
  while (context.nodeTest("for-each").size() > 0) context = context.nodeTest("for-each").front();
  context = context.nodeTest("yaxis").front();
  std::string yname = context.getAttribute("name");
  std::string ytype = strToLower(context.getAttribute("type"));
  std::string ylabel = yname;
  if (ytype != "parameter" && ytype != "mean") ylabel = ytype + "(" + ylabel + ")";
  if (ytype == "index") ytype = "indexvalue";
  std::string yalias;
  bool yerror = (context.getAttribute("error") != "");
  if (ytype == "parameter") {
    if (parameter_aliases.find(yname) == parameter_aliases.end()) {
      // assign a new alias name for parameter
      yalias = std::string("x") + boost::lexical_cast<std::string>(alias_number++);
      append_from(axisSQLFrom, "parameter AS " + yalias);
      if (yalias != "x0") append_where(axisSQLWhere, "x0.fID=" + yalias + ".fID");
      parameter_aliases[yname] = yalias;
    } else {
      yalias = parameter_aliases[yname];
    }
    axisSQLOrder = yalias + ".value";
  } else {
    if (measurement_aliases.find(yname) == measurement_aliases.end()) {
      // assign a new alias name for measurement
      yalias = std::string("x") + boost::lexical_cast<std::string>(alias_number++);
      append_from(axisSQLFrom, "measurement AS " + yalias);
      if (yalias != "x0") append_where(axisSQLWhere, "x0.fID=" + yalias + ".fID");
      measurement_aliases[yname] = yalias;
    } else {
      yalias = measurement_aliases[yname];
    }
    if (context.getAttribute("index") != "") {
      append_where(axisSQLWhere, yalias + ".indexvalue=\""
                   + SQLite::quote(context.getAttribute("index")) + "\"");
      ylabel += " " + context.getAttribute("index");
    }
}
  append_where(axisSQLWhere, yalias + ".name='" + SQLite::quote(yname) + "'");
  if (ytype == "parameter") ytype = "value";

  // Parameters cannot have indeces or errors
  if (((xtype == "value" || xtype == "indexvalue") && xerror) ||
      ((ytype == "value" || ytype == "indexvalue") && yerror))
    throw std::runtime_error("Only measurements can have errorbars and indeces");

  // std::cout << "axisSQLFrom: " << axisSQLFrom << std::endl
  //           << "axisSQLWhere: " << axisSQLWhere << std::endl
  //           << "axisSQLOrder: " << axisSQLOrder << std::endl;

  std::string SQLFrom = constraintSQLFrom;
  append_from(SQLFrom, foreachSQLFrom);
  append_from(SQLFrom, axisSQLFrom);

  std::string SQLWhere = constraintSQLWhere;
  append_where(SQLWhere, foreachSQLWhere);
  append_where(SQLWhere, axisSQLWhere);

  std::string SQLOrder = axisSQLOrder;

  std::string SQLSelect("SELECT " + xalias + "." + xtype + " AS x");
  if (xerror) SQLSelect += ", " + xalias + ".error AS ex";
  SQLSelect += ", " + yalias + "." + ytype + " AS y";
  if (yerror) SQLSelect += ", " + yalias + ".error AS ey";

  //
  // Loop over <for-each> and output results
  //

  create_header(buffer, output_type, plot_name, xlabel, ylabel);

  if (foreachName.size() == 0) {
    if (SQLFrom.size()) SQLFrom = " FROM " + SQLFrom;
    if (SQLWhere.size()) SQLWhere = " WHERE " + SQLWhere;
    if (SQLOrder.size()) SQLOrder = " ORDER BY " + SQLOrder;
    result_type rs = mDB(SQLSelect + SQLFrom + SQLWhere + SQLOrder +";", true);
    check_duplicate(rs);
    write_result(buffer, bufferBody, rs, output_type, 0, xlabel, ylabel, xerror, yerror,
                 foreachName, foreachData, foreachOffset);
  } else {
    if (SQLFrom.size()) SQLFrom = " FROM " + SQLFrom;
    if (SQLOrder.size()) SQLOrder = " ORDER BY " + SQLOrder;
    for (int p = 0; p < num_plots; ++p) {
      std::string where = SQLWhere;
      int q = p;
      for (std::size_t v = 0; v < foreachData.size(); ++v, q = q % foreachOffset[v-1]) {
        std::string name = foreachName[v];
        append_where(where, parameter_aliases[name] + ".name='" + SQLite::quote(name) + "'");
        append_where(where, parameter_aliases[name] + ".value='"
                     + foreachData[v][q / foreachOffset[v]] + "'");
      }
      if (where.size()) where = " WHERE " + where;
      result_type rs = mDB(SQLSelect + SQLFrom + where + SQLOrder +";", true);
      check_duplicate(rs);
      write_result(buffer, bufferBody, rs, output_type, p, xlabel, ylabel, xerror, yerror,
                   foreachName, foreachData, foreachOffset);
    }
  }

  buffer += bufferBody;
  create_trailer(buffer, output_type);

  // Datei schreiben
  writeFile(mOutPath / (inInFile + "." + output_type), buffer);
}
