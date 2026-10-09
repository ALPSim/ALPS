/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_HDF5_STD_MAP
#define ALPS_NGS_HDF5_STD_MAP

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/cast.hpp>

#include <map>

namespace alps {
    namespace hdf5 {

        template <typename K, typename T, typename C, typename A> void save(
              archive & ar
            , std::string const & path
            , std::map<K, T, C, A> const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            if (ar.is_group(path))
                ar.delete_group(path);
            for(typename std::map<K, T, C, A>::const_iterator it = value.begin(); it != value.end(); ++it)
                save(ar, ar.complete_path(path) + "/" + ar.encode_segment(cast<std::string>(it->first)), it->second);
        }

        template <typename K, typename T, typename C, typename A> void load(
              archive & ar
            , std::string const & path
            , std::map<K, T, C, A> & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            std::vector<std::string> children = ar.list_children(path);
            for (typename std::vector<std::string>::const_iterator it = children.begin(); it != children.end(); ++it)
                load(ar, path + "/" +  *it, value[cast<K>(ar.decode_segment(*it))]);
        }
    }
}

#endif
