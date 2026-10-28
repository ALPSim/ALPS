/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2009 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/parser/xslt_path.h>
#include <alps/config.h>
#include <alps/version.h>
#include <boost/throw_exception.hpp>
#include <filesystem>
#include <stdexcept>

std::string alps::xslt_path(const std::string& stylefile) {
  char* p = getenv("ALPS_XML_PATH");
  if (p==0)
    return stylefile == "job.xsl" || stylefile == "plot2html.xsl" ? "ALPS.xsl" : stylefile;
  
  std::string path = p;
  
  if (path != "http://xml.comp-phys.org" && path != "http://xml.comp-phys.org/")
    return path+"/"+stylefile;
  else if (stylefile == "ALPS.xsl")
    return "http://xml.comp-phys.org/2009/12/ALPS.xsl";
  else
    return "http://xml.comp-phys.org/"+stylefile;
}


std::string alps::search_xml_library_path(const std::string& file)
{
  // check whether this is a full path and the file exists 
  std::filesystem::path path(file);
  if (std::filesystem::exists(path))
    return path.string();
    
  char* p = getenv("ALPS_XML_PATH");
  if (p != 0)
    path = std::filesystem::path(p) /  file;
  else { 
    p = getenv("ALPS_ROOT");
    if (p != 0) 
      path = std::filesystem::path(p)
           / std::filesystem::path(ALPS_XML_DIR).lexically_relative(ALPS_PREFIX) / file;
    else {
      path = std::filesystem::path(ALPS_XML_DIR) / file;
    }
  }
  if (!std::filesystem::exists(path))
    boost::throw_exception(std::runtime_error("Cannot find file " + file + ", failed to find it at " + path.string()));
  return path.string();
}

void alps::copy_stylesheet(std::filesystem::path const& dir)
{
  std::filesystem::path dest = dir / "ALPS.xsl";
  if (! std::filesystem::exists(dest))
    std::filesystem::copy_file(search_xml_library_path("ALPS.xsl"),dest);
}
