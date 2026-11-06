/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2008 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/utilities/os.hpp>
#include <alps/version.h>
#include <boost/throw_exception.hpp>
#include <filesystem>
#include <stdexcept>
#include <vector>

#include <iostream>

#include <boost/asio/ip/host_name.hpp>
#include <cstdlib>

namespace alps {

//=======================================================================
// hostname
//
// returns the host name
//-----------------------------------------------------------------------

std::string hostname() { return boost::asio::ip::host_name(); }

std::string username() {
  for (const char* key : {"LOGNAME", "USER", "USERNAME"})
    if (const char* value = std::getenv(key)) return value;
  return "unknown";
}

std::filesystem::path temp_directory_path() {
  return std::filesystem::temp_directory_path();
}

std::filesystem::path installation_directory()
{
  return std::filesystem::path(ALPS_PREFIX);
}

std::filesystem::path bin_directory()
{
  if (const char* bin = std::getenv("ALPS_BIN_PATH"))
    return std::filesystem::path(bin);
  return std::filesystem::path(ALPS_BIN_DIR);
}

  


} // end namespace alps
