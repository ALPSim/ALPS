/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2009 by Lukas Gamper <mistral@student.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>,
*                            Niall Moran <nmoran@thphys.nuim.ie>
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

#include <fstream>
#include <stdexcept>
#include <iostream>
#include <ctime>
#include <iomanip>
#include <boost/classic_spirit.hpp>
#include <boost/throw_exception.hpp>

namespace sp = boost::spirit;

#include "archive_xml.hpp"

namespace XMLActions {

Node mRoot;
Node *mContext;
std::string mCurAttrName;
Node getRoot() { return mRoot; }

void newDoc(char const*, char const*)
{
  mRoot = Node();
  mContext = &mRoot;
}

void startElem(char const* str, char const* end)
{
  std::string name = std::string(str, end);
  std::transform (name.begin(), name.end(), name.begin(), tolower);
  mContext->addElement(name);
  mContext = mContext->getLastElement();
}

void endElem(char const*, char const*) { mContext = mContext->getParent();  }

void newText(char const* str, char const* end)
{ mContext->addText(std::string(str, end)); }

void storeAttrName(char const* str, char const* end)
{ mCurAttrName = std::string(str, end); }

void newAttr(char const* str, char const* end)
{ mContext->addAttribute(mCurAttrName, std::string(str, end)); }

} // end namespace XMLActions

struct XMLEbnf : public sp::grammar<XMLEbnf>
{
  template <typename ScannerT>
  struct definition
  {
    sp::rule<ScannerT> Constructor, QName, DirCommentConstructor,
      DirAttributeValue;
    sp::rule<ScannerT> Header, DirElemContent, DirAttributeList,
      DirElemConstructor, DocType;

    definition(XMLEbnf const&)
    {
      Constructor
        = sp::eps_p[&XMLActions::newDoc]
          >> *(Header | DocType)
          >> *DirCommentConstructor
          >> DirElemConstructor;
      Header
        = sp::as_lower_d[sp::str_p("<?xml")]
          >> *( (sp::anychar_p - '?')
              | (sp::ch_p('?') >> (sp::anychar_p - '>'))
              )
          >> "?>";
      DocType
        = sp::as_lower_d[sp::str_p("<!doctype")]
          >> *(sp::anychar_p - '>')
          >> ">";
      DirElemConstructor
        = sp::ch_p('<')
          >> QName[&XMLActions::startElem]
          >> *DirAttributeList
          >> ( ("/>" >> sp::eps_p[&XMLActions::endElem])
             | ( '>'
                 >> *( DirElemContent[&XMLActions::newText]
                     | DirElemConstructor
                     | DirCommentConstructor
                     )
                 >> "</"
                 >> QName
                 >> !sp::space_p
                 >> '>'
                 >> sp::eps_p[&XMLActions::endElem]
               )
             );
      DirElemContent
        = +~sp::ch_p('<');
      DirAttributeList
        = QName[&XMLActions::storeAttrName]
          >> '='
          >> DirAttributeValue;
      DirAttributeValue
        = ( sp::ch_p('"')
            >> ( *( sp::str_p("\"\"")
                  | (sp::anychar_p - sp::ch_p('"'))
                  )
               )[&XMLActions::newAttr]
            >> '"'
          )
        | ( sp::ch_p('\'')
            >> ( *( sp::str_p("''")
                  | (sp::anychar_p - sp::ch_p('\''))
                  )
               )[&XMLActions::newAttr]
            >> '\''
          );
      DirCommentConstructor
        = sp::str_p("<!--")
          >> *( (sp::anychar_p - '-')
              | (sp::ch_p('-') >> (sp::anychar_p - '-'))
              )
          >> "-->";
      QName
        = sp::lexeme_d
          [
            +( sp::range_p('a','z')
             | sp::range_p('A','Z')
             | sp::range_p('0','9')
             | '-'
             | '_'
             | '.'
             )
            >> !( sp::ch_p(':')
                  >> +( sp::range_p('a','z')
                      | sp::range_p('A','Z')
                      | sp::range_p('0','9')
                      | '-'
                      | '_'
                      | '.'
                      )
                )
          ];
    }

    sp::rule<ScannerT> const& start() const { return Constructor; }
  };
};

struct PlotDTD : public sp::grammar<PlotDTD>
{
  template <typename ScannerT>
  struct definition
  {
    sp::rule<ScannerT> Constructor, Header, DocType, AttributeValue;
    sp::rule<ScannerT> Plot, ForEach, XAxis, YAxis, Constraint, Comment;
    sp::rule<ScannerT> AtName, AtOutput, AtType, AtLabel, AtError, AtIndex,
      AtOperator, AtValue;

    definition(PlotDTD const&)
    {
      Constructor
        = *(Header | DocType) >> Plot;
      Header
        = sp::as_lower_d[sp::str_p("<?xml")]
          >> *( (sp::anychar_p - '?')
              | (sp::ch_p('?') >> (sp::anychar_p - '>'))
              )
          >> "?>";
      DocType
        = sp::as_lower_d[sp::str_p("<!doctype")]
          >> *(sp::anychar_p - '>')
          >> ">";
      Plot
        = *Comment
          >> ( sp::as_lower_d[sp::str_p("<plot")]
               >> *(AtName | AtOutput)
               >> '>'
               >> (ForEach | (XAxis >> YAxis >> *Constraint ))
               >> *Comment
               >> sp::as_lower_d[sp::str_p("</plot")]
               >> '>'
             );
      ForEach
        = *Comment
          >> sp::as_lower_d[sp::str_p("<for-each")] >> AtName >> '>'
          >> (( *Comment >> *Constraint >> *Comment >> XAxis >> *Comment >> YAxis
                >> *Constraint >> *Comment ) | ForEach )
          >> sp::as_lower_d[sp::str_p("</for-each")] >> '>';
      XAxis
        = *Comment
          >> sp::as_lower_d[sp::str_p("<xaxis")]
          >> *(AtName | AtLabel | AtType | AtError | AtIndex)
          >> ( "/>"
             | ( '>'
                 >> *Comment
                 >> sp::as_lower_d[sp::str_p("</xaxis")]
                 >> '>'
               )
             );
      YAxis
        = *Comment
          >> sp::as_lower_d[sp::str_p("<yaxis")]
          >> *(AtName | AtLabel | AtType | AtError | AtIndex)
          >> ( "/>"
             | ( '>'
                 >> *Comment
                 >> sp::as_lower_d[sp::str_p("</yaxis")]
                 >> '>'
               )
             );
      Constraint
        = *Comment
          >> sp::as_lower_d[sp::str_p("<constraint")]
          >> *(AtName | AtType | AtOperator | AtValue)
          >> ( "/>"
             | ( '>'
                 >> *Comment
                 >> sp::as_lower_d[sp::str_p("</yaxis")]
                 >> '>'
               )
             );
      Comment
        = sp::str_p("<!--")
          >> *( (sp::anychar_p - '-')
              | (sp::ch_p('-') >> (sp::anychar_p - '-'))
              )
          >> "-->";
      AtName
        = sp::as_lower_d[sp::str_p("name")]
          >> '='
          >> AttributeValue;
      AtOutput
        = sp::as_lower_d[sp::str_p("output")]
          >> '='
          >> (sp::ch_p('"') | '\'')
          >> sp::as_lower_d[(sp::str_p("text") | "html" | "gnuplot" | "xmgr" | "xml")]
          >> (sp::ch_p('"') | '\'');
      AtType
        = sp::as_lower_d[sp::str_p("type")]
          >> '='
          >> (sp::ch_p('"') | '\'')
          >> sp::as_lower_d[( sp::str_p("count")
                            | "mean"
                            | "error"
                            | "variance"
                            | "autocorr"
                            | "index"
                            | "parameter"
                            )
                           ]
          >> (sp::ch_p('"') | '\'');
      AtLabel
        = sp::as_lower_d[sp::str_p("label")]
          >> '='
          >> AttributeValue;
      AtError
        = sp::as_lower_d[sp::str_p("error")]
          >> '='
          >> (sp::ch_p('"') | '\'')
          >> sp::as_lower_d[sp::str_p("true")]
          >> (sp::ch_p('"') | '\'');
      AtIndex
        = sp::as_lower_d[sp::str_p("index")]
          >> '='
          >> AttributeValue;
      AtOperator
        = sp::as_lower_d[sp::str_p("operator")]
          >> '='
          >> (sp::ch_p('"') | '\'')
          >> sp::as_lower_d[( sp::str_p("lessthan")
                            | "lessorequalthan"
                            | "greaterthan"
                            | "greaterorequalthan"
                            | "notequal"
                            | "equal"
                            )
                           ]
          >> (sp::ch_p('"') | '\'');
      AtValue
        = sp::as_lower_d[sp::str_p("value")]
          >> '='
          >> AttributeValue;
      AttributeValue
        = ( sp::ch_p('"')
            >> ( *( sp::str_p("\"\"")
                  | (sp::anychar_p - sp::ch_p('"'))
                  )
               )
            >> '"'
          )
        | ( sp::ch_p('\'')
            >> ( *( sp::str_p("''")
                  | (sp::anychar_p - sp::ch_p('\''))
                  )
               )
            >> '\''
         );
    }

    sp::rule<ScannerT> const& start() const { return Constructor; }
  };
};

std::string XML::readFile(fs::path filename)
{
  std::ifstream fin(filename.string().c_str());
  std::string buffer;
  if (!fin.good())
    boost::throw_exception(std::runtime_error(std::string(
      "Could not open file: ") + filename.string()));
  fin.unsetf(std::ios::skipws);
  copy(std::istream_iterator<char>(fin), std::istream_iterator<char>(),
       std::back_inserter(buffer));
  fin.close();
  return buffer;
}

Node XML::operator()(fs::path inFileName, bool usePlotDTD)
{
  std::clock_t timer = std::clock();

  std::string fileContent = readFile(inFileName);
  if (usePlotDTD) {
    sp::parse_info<> plotInfo =
      sp::parse(fileContent.c_str(), PlotDTD(), sp::space_p);
    if (!plotInfo.hit)
      boost::throw_exception(std::runtime_error(std::string(
        "The Plotfile does not match the DTD! Fails parsing at: ") +
        plotInfo.stop));
  }

  sp::parse_info<> info =
    sp::parse(fileContent.c_str(), XMLEbnf(), sp::space_p);
  if (mVerbose) {
    if (info.full) {
      std::cout << "File parsed: " << inFileName.string() << " ("
                << std::setiosflags(std::ios_base::fixed)
                << std::setprecision(2)
                << std::difftime(std::clock(), timer) / 1000000 << "s)"
                << std::endl;
    } else
      std::cerr << std::endl << "Fails parsing file " << inFileName.string()
                << " at " << info.stop << std::endl;
  }
  return XMLActions::getRoot();
}
