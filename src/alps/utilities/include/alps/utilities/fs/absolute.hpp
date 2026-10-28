/*
 * Copyright (C) 1998-2026 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */

#ifndef ALPS_UTILITY_FS_ABSOLUTE_HPP
#define ALPS_UTILITY_FS_ABSOLUTE_HPP

#include <filesystem>
#include <alps/utilities/fs/absolute.hpp>

namespace alps { namespace fs {
  /** @brief Resolves a path against a base directory, as Boost.Filesystem's
             two-argument absolute() did
      @param p : path to resolve; returned unchanged if it is absolute
      @param base : directory a relative @p p is taken from; itself made
             absolute against the current directory if it is relative
      @returns @p p if absolute, else absolute(@p base) / @p p
  */
  inline std::filesystem::path absolute(const std::filesystem::path& p,
                                        const std::filesystem::path& base) {
    return p.is_absolute() ? p : std::filesystem::absolute(base) / p;
  }
} }

#endif /* ALPS_UTILITY_FS_ABSOLUTE_HPP */
