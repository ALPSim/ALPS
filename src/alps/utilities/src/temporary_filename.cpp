// Copyright (C) 2008 - 2010 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#include <alps/utility/temporary_filename.hpp>
#include <boost/filesystem/operations.hpp>
#include <cerrno>
#include <cstdio>
#include <stdexcept>
#include <system_error>

namespace alps {
std::string temporary_filename(std::string prefix) {
    // C11 exclusive creation reserves the name, unlike mktemp. Close the
    // stream here: callers own the file, not an otherwise leaked descriptor.
    for (int attempt = 0; attempt < 128; ++attempt) {
        auto name = prefix + boost::filesystem::unique_path("%%%%%%%%%%%%").string();
        if (auto* file = std::fopen(name.c_str(), "wbx")) {
            if (std::fclose(file) != 0)
                throw std::system_error(errno, std::generic_category(), "close temporary file");
            return name;
        }
        if (errno != EEXIST)
            throw std::system_error(errno, std::generic_category(), "create temporary file");
    }
    throw std::runtime_error("Could not reserve a unique temporary filename");
}
}
