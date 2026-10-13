/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: convert2xml.C 3523 2009-12-12 05:52:24Z troyer $ */
#include <alps/config.h>

#include <string>

namespace alps {

  /// convert a file from XDR format to XML
  std::string ALPS_DECL convert2xml(std::string const& name);
} // end namespace
