/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2005 by Matthias Troyer <troyer@comp-phys.org>,
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

#ifndef ALPS_MODEL_SUBSTITUTE_H
#define ALPS_MODEL_SUBSTITUTE_H

#include <boost/algorithm/string/replace.hpp>

namespace alps {
  
inline std::string substitute(std::string const& text, unsigned int type)
{
  std::string n;
  for (unsigned int i=0;i<text.size();++i)
  if (text[i]=='#')
    n += boost::lexical_cast<std::string>(type);
  else
    n += text[i];
  // std::cerr << "Replaced " << text << " to " << n << " by substituting " << boost::lexical_cast<std::string>(type) << "\n";
  return n;
//  return boost::algorithm::replace_all_copy(text,"#",boost::lexical_cast<std::string>(type));
}
  
inline Parameters substitute(Parameters const& parms, unsigned int type)
{
  Parameters p;
  for (Parameters::const_iterator it = parms.begin() ; it != parms.end(); ++it)
    p[substitute(it->key(),type)] = substitute(it->value(),type);
  return p;
}

} // namespace alps

#endif
