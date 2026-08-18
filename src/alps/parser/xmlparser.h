/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_XML_XMLPARSER_H
#define ALPS_XML_XMLPARSER_H

#include <alps/config.h>
#include <alps/parser/xmlhandler.h>

#include <boost/filesystem/path.hpp>
#include <iosfwd>
#include <string>

#if defined(ALPS_HAVE_XERCES_PARSER)
# include <xercesc/parsers/SAXParser.hpp>
# include <xercesc/sax/HandlerBase.hpp>
#elif defined(ALPS_HAVE_EXPAT_PARSER)
# include <expat.h>
#endif

namespace alps {

class ALPS_DECL XMLParser
{
public:
  XMLParser(XMLHandlerBase&);
  ~XMLParser();

  void parse(std::istream& is);
  void parse(const std::string& file);
  void parse(const boost::filesystem::path& file);

private:
  XMLParser();

#if defined(ALPS_HAVE_XERCES_PARSER)
  XERCES_CPP_NAMESPACE_QUALIFIER SAXParser* parser_;
  XERCES_CPP_NAMESPACE_QUALIFIER HandlerBase* handler_;
#elif defined(ALPS_HAVE_EXPAT_PARSER)
  XML_Parser parser_;
#else
  XMLHandlerBase& handler_;
#endif
};

} // end namespace alps

#endif // ALPS_XML_XMLPARSER_H
