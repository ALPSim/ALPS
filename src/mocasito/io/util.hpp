// Copyright (C) 2008 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include <string>
#include <vector>
#include <sstream>
#ifndef IO_UTIL
#define IO_UTIL
namespace mocasito {
    namespace io {
        #define MOCASITO_IO_FOREACH_SCALAR(callback)                                        \
            callback(char)                                                                    \
            callback(signed char)                                                            \
            callback(unsigned char)                                                            \
            callback(short)                                                                    \
            callback(unsigned short)                                                        \
            callback(long)                                                                    \
            callback(int)                                                                    \
            callback(unsigned)                                                                \
            callback(unsigned long)                                                            \
            callback(long long)                                                                \
            callback(unsigned long long)                                                    \
            callback(float)                                                                    \
            callback(double)                                                                \
            callback(long double)                                                            \
            callback(bool)
        namespace detail {
            typedef enum node_t { IO_UNKNOWN = -1, IO_GROUP, IO_DATA, IO_ATTR } node_t;
        }
    }
}
#endif
