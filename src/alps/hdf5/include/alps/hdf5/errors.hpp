/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_HDF5_ERROR_HPP
#define ALPS_NGS_HDF5_ERROR_HPP

#include <string>
#include <stdexcept>

namespace alps {
    namespace hdf5 {

        class archive_error : public std::runtime_error {
            public:
                archive_error(std::string const & what)
                    : std::runtime_error(what) 
                {}
        };

        #define DEFINE_ALPS_HDF5_EXCEPTION(name)                                    \
            class name : public archive_error {                                     \
                public:                                                             \
                    name (std::string const & what)                                 \
                        : archive_error(what)                                       \
                    {}                                                              \
            };
        DEFINE_ALPS_HDF5_EXCEPTION(archive_not_found)
        DEFINE_ALPS_HDF5_EXCEPTION(archive_closed)
        DEFINE_ALPS_HDF5_EXCEPTION(invalid_path)
        DEFINE_ALPS_HDF5_EXCEPTION(path_not_found)
        DEFINE_ALPS_HDF5_EXCEPTION(wrong_type)
        #undef DEFINE_ALPS_HDF5_EXCEPTION
    }
};

#endif
