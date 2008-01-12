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

#include <alps/xml.h>
#include <alps/vectorio.h>
#include "types.h"
#include "green_function.h"

template <class T>
void write_freq(alps::oxstream& xml, std::vector<T> const& v, int flavor)
{
  
  xml << alps::start_tag("VECTOR") 
      << alps::attribute("index","Matsubara frequency")
      << alps::attribute("flavor",flavor)
      << alps::attribute("size",v.size());
  xml << alps::write_vector(v," ");
  xml << alps::end_tag();
}

template <class T>
void write_freq(alps::oxstream& xml, std::pair<T,T> const& p)
{
  write_freq(xml,p.first, 0);
  write_freq(xml,p.second, 1);
}

template <class T>
void write_freq(alps::oxstream& xml, green_function<T> const& g)
{
    for(uint s1=0;s1<g.nsite();++s1){
      for(uint s2=0;s2<g.nsite();++s2){
	for(uint f=0;f<g.nflavor();++f){
        std::vector<std::complex<double> > g_f_s1_s2;
        for(uint w=0;w<g.nfreq();++w){
          g_f_s1_s2.push_back(g(w,s1,s2,f));
        }
        write_freq(xml,g_f_s1_s2, f);
      }
    }
  }
}

template <class T>
void write_itime(alps::oxstream& xml, std::vector<T> const& v, int flavor)
{
  xml << alps::start_tag("VECTOR") 
      << alps::attribute("index","imaginary time")
      << alps::attribute("flavor",flavor)
      << alps::attribute("size",v.size());
  xml << alps::write_vector(v," ");
  xml << alps::end_tag();
}


template <class T>
void write_itime(alps::oxstream& xml, std::pair<T,T> const& p)
{
  write_itime(xml,p.first, 0);
  write_itime(xml,p.second, 1);
}
template <class T>
void write_itime(alps::oxstream& xml, green_function<T> const& g)
{
    for(uint s1=0;s1<g.nsite();++s1){
      for(uint s2=0;s2<g.nsite();++s2){
	for(uint f=0;f<g.nflavor();++f){
        std::vector<double> g_f_s1_s2;
        for(uint t=0;t<g.ntime();++t){
          g_f_s1_s2.push_back(g(t,s1,s2,f));
        }
        write_itime(xml,g_f_s1_s2, f);
      }
    }
  }
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
void read_itime(std::istream& in, green_function<T>& g)
{
  for(uint s1=0;s1<g.nsite();++s1){
    for(uint s2=0;s2<g.nsite();++s2){
      for(uint f=0;f<g.nflavor();++f){
        std::vector<double> g_f_s1_s2;
        g_f_s1_s2.resize(g.ntime());
        read_itime(in,g_f_s1_s2);
        for(uint t=0;t<g.ntime();++t){
          g(t,s1,s2,f)=g_f_s1_s2[t];
        }
      }
    }
  }
}

template <class T>
void read_freq(std::istream& in, std::vector<T>& v)
{
  alps::XMLTag tag = alps::parse_tag(in);
  if (tag.attributes["index"] !="Matsubara frequency")
    boost::throw_exception(std::runtime_error("Vector not in Matsubara frequency"));
  v.resize(atoi(tag.attributes["size"].c_str()));
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
template <class T>
void read_freq(std::istream& in, green_function<T>& g)
{
    for(uint s1=0;s1<g.nsite();++s1){
      for(uint s2=0;s2<g.nsite();++s2){
	for(uint f=0;f<g.nflavor();++f){
        std::vector<std::complex<double> > g_f_s1_s2;
        g_f_s1_s2.resize(g.nfreq());
        read_freq(in,g_f_s1_s2);
        for(uint w=0;w<g.nfreq();++w){
          g(w,s1,s2,f)=g_f_s1_s2[w];
        }
      }
    }
  }
}


#endif // ALPS_DMFT_XML_H
