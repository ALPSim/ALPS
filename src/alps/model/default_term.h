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

#ifndef ALPS_MODEL_DEFAULTTERM_H
#define ALPS_MODEL_DEFAULTTERM_H

#include <alps/model/siteterm.h>
#include <alps/model/bondterm.h>
#include <alps/model/substitute.h>

namespace alps {

template <class TERM>
class DefaultTermDescriptor : public TERM
{
public:
  typedef TERM term_type;
  DefaultTermDescriptor() {}
  DefaultTermDescriptor(const XMLTag& tag, std::istream& in) : term_type(tag,in) {}
  // operator term_type() const { return static_cast<term_type const&>(*this);}
  term_type get(unsigned int type) const;
  Parameters parms(unsigned int type) const { return substitute(TERM::parms(),type); }
};

template <class TERM>
TERM DefaultTermDescriptor<TERM>::get(unsigned int type) const
{
  return term_type(*this,substitute(this->term(),type),Parameters(),type);
}

typedef DefaultTermDescriptor<SiteTermDescriptor> DefaultSiteTermDescriptor;
typedef DefaultTermDescriptor<BondTermDescriptor> DefaultBondTermDescriptor;



} // namespace alps

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace alps {
#endif

template <class TERM>
inline alps::oxstream& operator<<(alps::oxstream& out, const alps::DefaultTermDescriptor<TERM>& q)
{
  q.write_xml(out);
  return out;
}

template <class TERM>
inline std::ostream& operator<<(std::ostream& out, const alps::DefaultTermDescriptor<TERM>& q)
{
  alps::oxstream xml(out);
  xml << q;
  return out;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // namespace alps
#endif

#endif
