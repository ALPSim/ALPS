/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2006 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
