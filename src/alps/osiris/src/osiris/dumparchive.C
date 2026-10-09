/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2009 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                       Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <boost/version.hpp>

#if BOOST_VERSION >= 104000

#define BOOST_ARCHIVE_SOURCE
#include <boost/archive/detail/archive_serializer_map.hpp>
#include <boost/archive/impl/archive_serializer_map.ipp>
#include <alps/osiris/dumparchive.h>

namespace boost {
namespace archive {

template class detail::archive_serializer_map<alps::odump_archive>;
template class detail::archive_serializer_map<alps::idump_archive>;

} // namespace archive
} // namespace boost

#endif //  BOOST_VERSION >= 104000
