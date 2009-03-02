// Copyright (C) 2008 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include <string>
#include <vector>
#include <stdexcept>
#include "util.hpp"
#include "../trace.hpp"
#ifndef IO_DETAIL_ITERATOR
#define IO_DETAIL_ITERATOR
namespace mocasito {
	namespace io {
		namespace detail {
			template<typename E, typename C> class iterator {
				public:
					iterator()
						: _closed(true)
					{
						MOCASITO_TRACE
					}
					iterator(iterator const & it)
						: _closed(it._closed), _index(_index), _list(it._list), _path(it._path), _engine(it._engine)
					{
						MOCASITO_TRACE					
					}
					iterator(E * e, std::string p, std::vector<std::string> l)
						: _closed(l.size() == 0), _index(0), _list(l), _path(p), _engine(e)
					{
						MOCASITO_TRACE					
					}
					iterator<E, C> operator++(int) {
						MOCASITO_TRACE
						iterator<E, C> it(*this);
						_closed = (++_index == _list.size());
						return it;
					}
					iterator<E, C> & operator++() {
						MOCASITO_TRACE
						_closed = (++_index == _list.size());
						return *this;
					}
					C operator*() {
						MOCASITO_TRACE
						return C(_engine, _path + "/" + _list[_index]);
					}
					bool operator==(iterator const & it) const {
						MOCASITO_TRACE
						return (_closed && it._closed) || (_closed == it._closed && _path == it._path && _index == it._index);
					}
					bool operator!=(iterator const & it) const {
						MOCASITO_TRACE
						return !(*this == it);
					}
				private:
					bool _closed;
					std::size_t _index;
					std::vector<std::string> _list;
					std::string _path;
					E * _engine;
			};
		}
	}
}
#endif

