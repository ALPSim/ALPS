/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2006 by Synge Todo <wistaria@comp-phys.org>
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

#include <alps/parser/xmlstream.h>

#include <cstdlib>
#include <stdexcept>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  double x = 3.14;

  alps::oxstream oxs;

  oxs << alps::header("MyEncoding");

  oxs << alps::stylesheet("URL to my stylesheet")
      << alps::processing_instruction("my_pi");

  oxs << alps::start_tag("tag0")
      << alps::attribute("name0", 1)

      << "this is a text"

      << alps::start_tag("tag1")
      << alps::start_tag("tag2")
      << alps::xml_namespace("MyNameSpace", "MyURL")
    
      << "text 2 "
      << "text 3 " << std::endl
      << alps::precision(3.14159265358979323846, 3) << ' '
      << alps::precision(3.14159265358979323846, 6) << '\n'
      << "text 4" << std::endl
      << alps::convert("text <&\">'")

      << alps::start_tag("tag3")
      << alps::end_tag

      << alps::precision(x, 6)

      << alps::start_tag("tag4") << alps::no_linebreak
      << "no linebreak"
      << alps::end_tag

      << alps::end_tag("tag2")
      << alps::end_tag("tag1")
      << alps::end_tag;

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exp) {
  std::cerr << exp.what() << std::endl;
  std::abort();
}
#endif
  return 0;
}
