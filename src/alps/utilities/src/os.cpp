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

#include <alps/utility/os.hpp>
#include <alps/version.h>
#include <boost/throw_exception.hpp>
#include <boost/filesystem/operations.hpp>
#include <stdexcept>
#include <vector>

#include <iostream>

#ifndef BOOST_POSIX_API
#include <windows.h>
#endif


#if defined(ALPS_HAVE_SYS_SYSTEMINFO_H)
# include <sys/systeminfo.h> // for sysinfo()
#elif defined(ALPS_HAVE_UNISTD_H)
# include <unistd.h>      // for gethostname() and getlogin()
# ifndef MAXHOSTNAMELEN
#   define MAXHOSTNAMELEN 256
# endif
#endif

#ifdef BOOST_NO_STDC_NAMESPACE
  namespace std { using ::system; }
#endif

namespace alps {

//=======================================================================
// hostname
//
// returns the host name
//-----------------------------------------------------------------------

std::string hostname()
{
#if defined(ALPS_HAVE_SYS_SYSTEMINFO_H)
  // the UNIX System V version
  char n[256];
  if(sysinfo(SI_HOSTNAME,n,256) < 0)
    boost::throw_exception(std::runtime_error("call to sysinfo failed in get_host_name"));
  return n;
#elif defined(ALPS_HAVE_UNISTD_H)
  // the BSD UNIX version
  char n[MAXHOSTNAMELEN];
  if(gethostname(n,MAXHOSTNAMELEN))
    boost::throw_exception(std::runtime_error("call to gethostname failed in get_host_name"));
  return n;
#else
  // the default version
  return "unnamed";
#endif
}

//=======================================================================
// username
//
// returns the username
//-----------------------------------------------------------------------

std::string username() {
#if defined(ALPS_HAVE_UNISTD_H)
  const char* login = getlogin();
  return (login ? std::string(getlogin()) : std::string("unknown"));
#else
  return std::string("unknown");
#endif
}


   // contributed by Jeff Flinn, from Boost 1.46
boost::filesystem::path temp_directory_path()
{
  using namespace boost::filesystem;
#   ifdef BOOST_POSIX_API
    const char* val = 0;
    
    (val = std::getenv("TMPDIR" )) ||
    (val = std::getenv("TMP"    )) ||
    (val = std::getenv("TEMP"   )) ||
    (val = std::getenv("TEMPDIR"));
    
    path p((val!=0) ? val : "/tmp");
    
    if (p.empty() || !is_directory(p))
      p=path(".");
      
    return p;
    
#   else  // Windows

    std::vector<TCHAR> buf(MAX_PATH);

    if (buf.empty() || GetTempPath(buf.size(), &buf[0])==0)
      return path(".");
        
    buf.pop_back();
    
    path p(buf.begin(), buf.end());
        
    if (!is_directory(p))
      p=path(".");

    
    return p;
#   endif
}

#ifdef ALPS_PREFIX
boost::filesystem::path installation_directory()
{
  return boost::filesystem::path(ALPS_PREFIX);
}
#else
boost::filesystem::path installation_directory()
{
  return boost::filesystem::path("");
}

#endif

boost::filesystem::path bin_directory()
{
  const char* val = 0;
  if (const char* bin = std::getenv("ALPS_BIN_PATH"))
    return boost::filesystem::path(bin);
  return installation_directory() / "bin";
}

  


} // end namespace alps
