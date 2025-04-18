/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Beat Ammon <ammon@ginnan.issp.u-tokyo.ac.jp>,
*                            Andreas Laeuchli <laeuchli@itp.phys.ethz.ch>,
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

#ifndef ALPS_ALEA_OUTPUT_HELPER_H
#define ALPS_ALEA_OUTPUT_HELPER_H

#include <alps/config.h>
#include <alps/xml.h>
#include <boost/filesystem/path.hpp>
#include <boost/mpl/bool.hpp>
#include <iostream>

namespace alps {
template <typename FLAG> struct output_helper {};


template <>
struct output_helper<boost::mpl::true_>
{
  template <class X, class L> static void output(const X& b, std::ostream& out, const L&)
  {
    b.output_scalar(out);
  }

  template <class X> static void output(const X& b, std::ostream& out)
  {
    b.output_scalar(out);
  }

  template <class X> static void write_xml(const X& b, oxstream& oxs, const boost::filesystem::path& fn_hdf5)
  {
    b.write_xml_scalar(oxs, fn_hdf5);
  }
  template <class X, class IT> static void write_more_xml(const X& b, oxstream& oxs, IT)
  {
    b.write_scalar_xml(oxs);
  }
};

template <>
struct output_helper<boost::mpl::false_>
{
  template <class T, class L> static void output(const T& b, std::ostream& out, const L& label)
  {
    b.output_vector(out,label);
  }

  template <class T> static void output(const T& b, std::ostream& out)
  {
    b.output_vector(out);
  }

  template <class T> static void write_xml(const T& b, oxstream& oxs, const boost::filesystem::path& fn_hdf5)
  {
    b.write_xml_vector(oxs, fn_hdf5);
  }
  
  template <class X, class IT> static void write_more_xml(const X& b, oxstream& oxs, IT i)
  {
    b.write_vector_xml(oxs, i);
  }
};

} // end namespace alps


#endif // ALPS_ALEA_OUTPUT_HELPER_H
