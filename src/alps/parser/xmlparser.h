/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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
