// Copyright (C) 2008 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include "util.hpp"
#include "traits.hpp"
#include <string>
#include <cstring>
#include <vector>
#include <deque>
#include <stdexcept>
#include <hdf5.h>
#ifndef IO_CACHE
#define IO_CACHE
namespace mocasito {
	namespace io {
		namespace detail {
			struct cache_index {
				cache_index(bool is_g = false, bool is_s = false, bool is_n = false, bool is_a = false)
					: is_group(is_g), is_scalar(is_s), is_null(is_n), is_attr(is_a), type(0), size(1), offset(0), children(0), attr(0) 
				{}
				bool is_group;
				bool is_scalar;
				bool is_null;
				bool is_attr;
				std::size_t type;
				std::size_t size;
				std::size_t offset;
				std::vector<std::string> children;
				std::vector<std::string> attr;
			};
		}
		template<typename InEngine, typename OutEngine> class cache {
			public:
				cache(std::string const & p, std::string const & q)
					: _out(q), _mem(0) 
				{
					MOCASITO_TRACE
					InEngine input(p, "");
					std::deque<std::string> stack(1, "/");
					std::string path;
					std::vector<std::string> list;
					while (stack.size()) {
						_index.insert(std::make_pair(path = stack.front(), detail::cache_index()));
						stack.pop_front();
						if (input.is_group(path)) {
							list = input.list_children(path);
							for	(std::vector<std::string>::iterator it = list.begin(); it != list.end(); ++it)
								stack.push_back(path + (path == "/" ? "" : "/") + *it);
						}
/*						list = input.list_attr(path);
						for	(std::vector<std::string>::iterator it = list.begin(); it != list.end(); ++it)
							_index.insert(std::make_pair(path + (path == "/" ? "@" : "/@") + *it, detail::cache_index(false, false, false, true)));
*/					}
					for (std::map<std::string, detail::cache_index>::iterator it = _index.begin(); it != _index.end(); ++it) {
						it->second.attr = input.list_attr(it->first);
						if (input.is_group(it->first)) {
							it->second.is_group = true;
							it->second.children = input.list_children(it->first);
						} else {
							it->second.type = input.datatype(it->first);
							std::size_t size_of;
							switch (it->second.type) {
								#define MOCASITO_IO_CACHE_GET_SIZE_OF(T)					\
									case type_traits<T>::value: size_of = sizeof(T); break;
								MOCASITO_IO_FOREACH_SCALAR(MOCASITO_IO_CACHE_GET_SIZE_OF)
								#undef MOCASITO_IO_CACHE_GET_SIZE_OF
								default: 
									MOCASITO_IO_THROW("unknown type")
							}
							it->second.offset = _mem.size();
							if (input.is_scalar(it->first)) {
								it->second.is_scalar = true;
								_mem.resize(_mem.size() + size_of);
							} else {
								it->second.is_null = input.is_null(it->first);
								if(input.dimensions(it->first) > 1)
									MOCASITO_IO_THROW("the path has more than one dimension: " + it->first)
								it->second.size = input.extent(it->first)[0];
								_mem.resize(_mem.size() + it->second.size * size_of);
							}
							switch (it->second.type) {
								#define MOCASITO_IO_CACHE_READ_DATA(T)						\
									case type_traits<T>::value: input.get_data(it->first, reinterpret_cast<T *>(&_mem[it->second.offset])); break;
								MOCASITO_IO_FOREACH_SCALAR(MOCASITO_IO_CACHE_READ_DATA)
								#undef MOCASITO_IO_CACHE_READ_DATA
							}
						}
					}
				}
				void flush() const {
					MOCASITO_TRACE
					OutEngine input(_out, _out);
					MOCASITO_IO_THROW("not implemented")
				}
				bool is_group(std::string const & p) const {
					MOCASITO_TRACE
					if (_index.find(p) == _index.end())
						MOCASITO_IO_THROW("the paht does not exists: " + p)
					return _index.find(p)->second.is_group;
				}
				bool is_data(std::string const & p) const {
					if (_index.find(p) == _index.end())
						MOCASITO_IO_THROW("the paht does not exists: " + p)
					return !_index.find(p)->second.is_group;
				}
				std::vector<std::size_t> extent(std::string const & p) const {
					MOCASITO_TRACE
					if (_index.find(p) == _index.end())
						MOCASITO_IO_THROW("the paht does not exists: " + p)
					if (_index.find(p)->second.is_scalar)
						return std::vector<std::size_t>(0);
					else if (_index.find(p)->second.is_null)
						return std::vector<std::size_t>(1, 0);
					else
						return std::vector<std::size_t>(1, _index.find(p)->second.size);
				}
				std::size_t dimensions(std::string const & p) const {
					MOCASITO_TRACE
					if (_index.find(p) == _index.end())
						MOCASITO_IO_THROW("the paht does not exists: " + p)
					if (_index.find(p)->second.is_scalar)
						return 0;
					else
						return 1;
				}
				type_traits<>::type attrtype(detail::node_t t, std::string const & p, std::string const & s) const {
					MOCASITO_IO_THROW("not implemented")
					return type_traits<>::value;
				}
				type_traits<>::type datatype(std::string const & p) const {
					MOCASITO_TRACE
					if (_index.find(p) == _index.end())
						MOCASITO_IO_THROW("the paht does not exists: " + p)
					return _index.find(p)->second.type;
				}
				bool is_scalar(std::string const & p) const {
					MOCASITO_TRACE
					if (_index.find(p) == _index.end())
						MOCASITO_IO_THROW("the paht does not exists: " + p)
					return _index.find(p)->second.is_scalar;
				}
				bool is_null(std::string const & p) const {
					MOCASITO_TRACE
					if (_index.find(p) == _index.end())
						MOCASITO_IO_THROW("the paht does not exists: " + p)
					return _index.find(p)->second.is_null;
				}
				std::vector<std::string> list_children(std::string const & p) const {
					MOCASITO_TRACE
					if (!is_group(p))
						MOCASITO_IO_THROW("the paht is not a group: " + p)
					return _index.find(p)->second.children;
				}
				template<typename T> void get_data(std::string const & p, T * v) const {
					MOCASITO_TRACE
					if (!is_data(p))
						MOCASITO_IO_THROW("the paht does not contains data: " + p)
					switch (_index.find(p)->second.type) {
						#define MOCASITO_IO_CACHE_GET_DATA(U)								\
							case type_traits<U>::value:										\
								{															\
									std::deque<U> buffer(_index.find(p)->second.size);		\
									std::memcpy(&buffer.front(), &_mem[_index.find(p)->second.offset], _index.find(p)->second.size * sizeof(U));\
									std::copy(buffer.begin(), buffer.end(), v);				\
								}															\
							break;
						MOCASITO_IO_FOREACH_SCALAR(MOCASITO_IO_CACHE_GET_DATA)
						#undef MOCASITO_IO_CACHE_GET_DATA
					}
				}
				template<typename T> void get_group_attr(std::string const & p, std::string const & s, T & v) const {
					MOCASITO_TRACE
					MOCASITO_IO_THROW("not implemented")
				}
				template<typename T> void get_data_attr(std::string const & p, std::string const & s, T & v) const {
					MOCASITO_TRACE
					MOCASITO_IO_THROW("not implemented")
				}
				template<typename T> void set_data(std::string const & p, T const & v) {
					MOCASITO_TRACE
					MOCASITO_IO_THROW("not implemented")
				}
				template<typename T> void set_data(std::string const & p, T const * v, hsize_t s) {
					MOCASITO_TRACE
					MOCASITO_IO_THROW("not implemented")
				}
				template<typename T> void append_data(std::string const & p, T const * v, hsize_t s) {
					MOCASITO_TRACE
					MOCASITO_IO_THROW("not implemented")
				}
				void delete_data(std::string const & p, std::string const & s) {
					MOCASITO_TRACE
					MOCASITO_IO_THROW("not implemented")
				}
				template<typename T> void set_group_attr(std::string const & p, std::string const & s, T const & v) {
					MOCASITO_TRACE
					MOCASITO_IO_THROW("not implemented")
				}
				template<typename T> void set_data_attr(std::string const & p, std::string const & s, T const & v) {
					MOCASITO_TRACE
					MOCASITO_IO_THROW("not implemented")
				}
			private:
				std::string _out;
				std::map<std::string, detail::cache_index> _index;
				std::vector<char> _mem;
		};
	}
}
#endif
