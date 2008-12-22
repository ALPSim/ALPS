// Copyright (C) 2008 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#ifndef IO_UTIL
#define IO_UTIL
namespace mocasito {
	namespace io {
		namespace detail {
			typedef enum node_t { IO_UNKNOWN = -1, IO_GROUP, IO_DATA, IO_GROUP_ATTR, IO_DATA_ATTR } node_t;
		}
	}
}
#endif
