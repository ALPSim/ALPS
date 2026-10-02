/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2004 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/parser/xmlparser.h>
#include <alps/parser/parser.h>

#include <boost/filesystem/fstream.hpp>
#include <fstream>
#include <istream>
#include <string>

namespace alps {

XMLParser::XMLParser(XMLHandlerBase& h) : handler_(h) {}

XMLParser::~XMLParser() {}

void XMLParser::parse(std::istream& in)
{
  while (in) {
    char c;
    in >> c;
    if (!in)
      break;
    in.putback(c);
    if (c=='<') {
      XMLTag tag = parse_tag(in, false);
      if (tag.type == XMLTag::OPENING || tag.type== XMLTag::SINGLE) {
        // start tag
        handler_.start_element(tag.name, tag.attributes, xml::element);
      }
      if (tag.type == XMLTag::CLOSING || tag.type== XMLTag::SINGLE) {
        // end tag
        if (tag.type == XMLTag::CLOSING)
          tag.name.erase(0,1);
        handler_.end_element(tag.name, xml::element);
      }
      if (tag.type == XMLTag::PROCESSING) {
        // processing instruction
        if (tag.name == "xml-stylesheet") {
          handler_.start_element(tag.name, tag.attributes, xml::stylesheet);
          handler_.end_element(tag.name, xml::stylesheet);
        } else {
          handler_.start_element(tag.name, tag.attributes,
                                 xml::processing_instruction);
          handler_.end_element(tag.name, xml::processing_instruction);
        }
      }
    }
    else {
      std::string t = parse_content(in);
      int pi = 0;
      for (std::size_t p = 0; p <= t.length(); ++p) {
        if ( p == t.length() || t[p] == '\n') {
          std::string s = t.substr(pi, p - pi);
          // remove preceding and following blanks
          s = s.erase(0, s.find_first_not_of(' '));
          s = s.erase(s.find_last_not_of(' ') + 1);
          if (s.size()) handler_.text(s);
          pi = p + 1;
        }
      }
    }
  }
}
void XMLParser::parse(const std::string& file) {
  std::ifstream is(file.c_str());
  parse(is);
}
void XMLParser::parse(const boost::filesystem::path& file) {
  boost::filesystem::ifstream is(file);
  parse(is);
}

} // namespace alps
