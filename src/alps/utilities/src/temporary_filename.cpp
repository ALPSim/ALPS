// Copyright (C) 2008 - 2010 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#include <alps/utilities/temporary_filename.hpp>
#include <cerrno>
#include <cstdio>
#include <stdexcept>
#include <system_error>
#ifndef _WIN32
#include <cstdlib>
#include <unistd.h>
#else
#include <random>
#endif

namespace alps {
std::string temporary_filename(std::string prefix) {
#ifndef _WIN32
    // mkstemp reserves the name atomically with owner-only permissions,
    // independently of the caller's umask.
    prefix += "XXXXXX";
    int descriptor = mkstemp(prefix.data());
    if (descriptor < 0)
        throw std::system_error(errno, std::generic_category(), "create temporary file");
    if (close(descriptor) != 0)
        throw std::system_error(errno, std::generic_category(), "close temporary file");
    return prefix;
#else
    // C11 exclusive creation reserves the name, unlike mktemp. Close the
    // stream here: callers own the file, not an otherwise leaked descriptor.
    std::random_device seed;
    std::mt19937 engine(seed());
    std::uniform_int_distribution<int> digit(0, 15);
    for (int attempt = 0; attempt < 128; ++attempt) {
        std::string name = prefix;
        for (int i = 0; i < 12; ++i)
            name += "0123456789abcdef"[digit(engine)];
        if (auto* file = std::fopen(name.c_str(), "wbx")) {
            if (std::fclose(file) != 0)
                throw std::system_error(errno, std::generic_category(), "close temporary file");
            return name;
        }
        if (errno != EEXIST)
            throw std::system_error(errno, std::generic_category(), "create temporary file");
    }
    throw std::runtime_error("Could not reserve a unique temporary filename");
#endif
}
}
