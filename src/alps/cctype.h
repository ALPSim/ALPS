/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

/// \file cctype.h
/// \brief A safe version of the standard cctype header
///
///  Some cctype headers do not undefine harmful macros, so undefine
///  them here.

#ifndef ALPS_CCTYPE_H
#define ALPS_CCTYPE_H

#include <cctype>

#ifdef isspace 
# undef isspace
#endif
#ifdef isprint
# undef isprint
#endif
#ifdef iscntrl
# undef iscntrl
#endif
#ifdef isupper
# undef isupper
#endif
#ifdef islower
# undef islower
#endif
#ifdef isalpha
# undef isalpha
#endif
#ifdef isdigit
# undef isdigit
#endif
#ifdef ispunct
# undef ispunct
#endif
#ifdef isxdigit
# undef isxdigit
#endif
#ifdef isalnum
# undef isalnum
#endif
#ifdef isgraph
# undef isgraph
#endif
#ifdef toupper
# undef toupper
#endif
#ifdef tolower
# undef tolower
#endif

#endif // ALPS_CCTYPE_H
