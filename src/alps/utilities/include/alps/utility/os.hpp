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

#ifndef ALPS_UTILITY_OS_HPP
#define ALPS_UTILITY_OS_HPP

//=======================================================================
// This file includes low level functions which depend on the OS used
//=======================================================================

#include <alps/config.h>
#include <boost/filesystem/path.hpp>
#include <string>

namespace alps {

/// returns the hostname
ALPS_DECL std::string hostname();

/// returns the username
ALPS_DECL std::string username();

/// returns the username
ALPS_DECL boost::filesystem::path temp_directory_path();

/// returns the installation directory
ALPS_DECL boost::filesystem::path installation_directory();

/// returns the program directory
ALPS_DECL boost::filesystem::path bin_directory();

} // end namespace

#endif // ALPS_UTILITY_OS_HPP
