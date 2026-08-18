/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
