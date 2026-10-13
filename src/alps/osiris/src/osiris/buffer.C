/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2004 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/osiris/buffer.h>
#include <alps/osiris/dump.h>
#include <alps/osiris/std/vector.h>

#include <boost/throw_exception.hpp>
#include <stdexcept>

namespace alps {
namespace detail {

// write a few bytes
void Buffer::write_buffer(const void* p, size_type n)
{
  size_type write_pos=size();
  resize(write_pos+n);
  memcpy(&((*this)[write_pos]),p,n);
}


// read a few bytes and update position
void Buffer::read_buffer(void* p, size_type n)
{
  if(read_pos+n>size())
    boost::throw_exception(std::runtime_error("read past buffer"));
  memcpy(p,&((*this)[read_pos]),n);
  read_pos+=n;
}

} // end namespace detail
} // end namespace alps
