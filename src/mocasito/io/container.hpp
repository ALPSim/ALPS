// Copyright (C) 2008 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include "iterator.hpp"
#include "context.hpp"
#include "free_fun.hpp"
#include <map>
#include <string>
#ifndef IO_CONTAINER
#define IO_CONTAINER
namespace mocasito {
	namespace io {
		template<typename E> class container {
			public:
				typedef detail::iterator<E, context<E> > iterator;
				container(std::string const & pin, std::string const & pout = std::string(""))
					: _engine(pin, pout)
				{
					MOCASITO_TRACE
				}
				context<E> & operator[](std::string const & p) {
					MOCASITO_TRACE
					if (_map.find(p) == _map.end())
						_map.insert(std::make_pair(p, context<E>(&_engine, p)));
					return _map[p];
				}
				void flush() const {
					MOCASITO_TRACE
					_engine.flush();
				}
			private:
				std::map<std::string, context<E> > _map;
				E _engine;
		};
	}
}
#endif
