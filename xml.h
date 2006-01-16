/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2006 by Matthias Troyer <troyer@comp-phys.org>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

/// @file xml.h
/// @brief defines the parser for XML schemas for the Green's functions and Weiss fields
///
/// @TODO this will need to be adapted as we go to mult-orbital or cluster models

#ifndef ALPS_DMFT_XML_H
#define ALPS_DMFT_XML_H

#include "types.h"
#include <alps/xml.h>
#include <alps/vectorio.h>

template <class T>
void write_freq(alps::oxstream& xml, std::vector<T> const& v)
{
  
  xml << alps::start_tag("VECTOR") 
      << alps::attribute("index","Matsubara frequency")
      << alps::attribute("size",v.size());
  xml << alps::write_vector(v," ");
  xml << alps::end_tag();
}

template <class T>
void write_freq(alps::oxstream& xml, std::pair<T,T> const& p)
{
  write_freq(xml,p.first);
  write_freq(xml,p.second);
}

template <class T>
void write_itime(alps::oxstream& xml, std::vector<T> const& v)
{
  xml << alps::start_tag("VECTOR") 
      << alps::attribute("index","imaginary time")
      << alps::attribute("size",v.size());
  xml << alps::write_vector(v," ");
  xml << alps::end_tag();
}


template <class T>
void write_itime(alps::oxstream& xml, std::pair<T,T> const& p)
{
  write_itime(xml,p.first);
  write_itime(xml,p.second);
}

template <class T>
void read_itime(std::istream& in, std::vector<T>& v)
{
  alps::XMLTag tag = alps::parse_tag(in);
  if (tag.attributes["index"] !="imaginary time")
    boost::throw_exception(std::runtime_error("Vector not in imaginary time"));
  alps::read_vector(alps::parse_content(in),v);
  tag = alps::parse_tag(in);
  if (tag.type !=alps::XMLTag::CLOSING)
    boost::throw_exception(std::runtime_error("Missing closing tag after reading vector"));
}


template <class T>
void read_itime(std::istream& in, std::pair<T,T>& p)
{
  read_itime(in,p.first);
  read_itime(in,p.second);
}

template <class T>
void read_freq(std::istream& in, std::vector<T>& v)
{
  alps::XMLTag tag = alps::parse_tag(in);
  if (tag.attributes["index"] !="Matsubara frequency")
    boost::throw_exception(std::runtime_error("Vector not in Matsubara frequency"));
  alps::read_vector(alps::parse_content(in),v);
  tag = alps::parse_tag(in);
  if (tag.type !=alps::XMLTag::CLOSING)
    boost::throw_exception(std::runtime_error("Missing closing tag after reading vector"));

}


template <class T>
void read_freq(std::istream& in, std::pair<T,T>& p)
{
  read_freq(in,p.first);
  read_freq(in,p.second);
}


#endif // ALPS_DMFT_XML_H
