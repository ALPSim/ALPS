// Copyright (C) 2008 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include "traits.hpp"
#include "util.hpp"
#include "../trace.hpp"
#include <string>
#include <cstring>
#include <vector>
#include <deque>
#include <map>
#include <stdexcept>
#include <algorithm>
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
					: _in(p), _out(q), _mem(0)
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
						list = input.list_attr(path);
						for	(std::vector<std::string>::iterator it = list.begin(); it != list.end(); ++it)
							_index.insert(std::make_pair(path + (path == "/" ? "@" : "/@") + *it, detail::cache_index(false, false, false, true)));
					}
					for (std::map<std::string, detail::cache_index>::iterator it = _index.begin(); it != _index.end(); ++it) {
						if (it->second.is_attr) {
							std::string p = it->first.substr(0, it->first.find_last_of('/'));
							std::string s = it->first.substr(it->first.find_last_of('@') + 1);
							it->second.type = input.attrtype(p, s);
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
							_mem.resize(_mem.size() + size_of);
							switch (it->second.type) {
								#define MOCASITO_IO_CACHE_READ_DATA(T)						\
									case type_traits<T>::value: input.get_attr(p, s, *reinterpret_cast<T *>(&_mem[it->second.offset])); break;
								MOCASITO_IO_FOREACH_SCALAR(MOCASITO_IO_CACHE_READ_DATA)
								#undef MOCASITO_IO_CACHE_READ_DATA
							}
						} else if (input.is_group(it->first)) {
							it->second.attr = input.list_attr(it->first);
							it->second.is_group = true;
							it->second.children = input.list_children(it->first);
						} else {
							it->second.attr = input.list_attr(it->first);
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
									case type_traits<T>::value:								\
										input.get_data(it->first, reinterpret_cast<T *>(&_mem[it->second.offset]));\
										break;
								MOCASITO_IO_FOREACH_SCALAR(MOCASITO_IO_CACHE_READ_DATA)
								#undef MOCASITO_IO_CACHE_READ_DATA
							}
						}
					}
				}
				void flush() const {
					MOCASITO_TRACE
					OutEngine output(_in, _out);
					for (std::map<std::string, detail::cache_index>::const_iterator it = _index.begin(); it != _index.end(); ++it)
						if (is_data(it->first) || it->second.is_attr)
							switch (_index.find(it->first)->second.type) {
								#define MOCASITO_IO_CACHE_FLUSH_SCALAR(T)				\
									case type_traits<T>::value:							\
										if (it->second.is_scalar || it->second.is_attr)	\
											output.set_data(it->first, *reinterpret_cast<T const *>(&_mem[it->second.offset]));\
										else											\
											output.set_data(it->first, reinterpret_cast<T const *>(&_mem[it->second.offset]), it->second.size);\
										break;
								MOCASITO_IO_FOREACH_SCALAR(MOCASITO_IO_CACHE_FLUSH_SCALAR)
								#undef MOCASITO_IO_CACHE_FLUSH_SCALAR
							}
					output.flush();	
				}
				bool is_group(std::string const & p) const {
					MOCASITO_TRACE
					return _index.find(p) != _index.end() && _index.find(p)->second.is_group;
				}
				bool is_data(std::string const & p) const {
					MOCASITO_TRACE
					return _index.find(p) != _index.end() && !_index.find(p)->second.is_group && !_index.find(p)->second.is_attr;
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
				type_traits<>::type attrtype(std::string const & p, std::string const & s) const {
					MOCASITO_TRACE
					if (_index.find(p + "/@" + s) == _index.end())
						MOCASITO_IO_THROW("the paht does not exists: " + p + "/@" + s)
					return _index.find(p + "/@" + s)->second.type;
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
				template<typename T> void get_attr(std::string const & p, std::string const & s, T & v) const {
					MOCASITO_TRACE
					if (_index.find(p + "/@" + s) == _index.end() || !_index.find(p + "/@" + s)->second.is_attr)
						MOCASITO_IO_THROW("the paht does not contains attribute data: " + p + "/@" + s)
					switch (_index.find(p)->second.type) {
						#define MOCASITO_IO_CACHE_GET_DATA(U)								\
							case type_traits<U>::value:										\
								std::memcpy(&v, &_mem[_index.find(p + "/@" + s)->second.offset], sizeof(U));\
							break;
						MOCASITO_IO_FOREACH_SCALAR(MOCASITO_IO_CACHE_GET_DATA)
						#undef MOCASITO_IO_CACHE_GET_DATA
					}
				}
				template<typename T> void set_data(std::string const & p, T const & v) {
					MOCASITO_TRACE
					if (_index.find(p) == _index.end())
						create_path(p);
					_index.find(p)->second.is_scalar = true;
					_index.find(p)->second.type = type_traits<T>::value;
					_index.find(p)->second.size = 0;
					_index.find(p)->second.offset = _mem.size();
					_mem.resize(_mem.size() + sizeof(T));
					std::memcpy(&_mem[_index.find(p)->second.offset], &v, sizeof(T));
				}
				template<typename T> void set_data(std::string const & p, T const * v, std::size_t s) {
					MOCASITO_TRACE
					if (_index.find(p) == _index.end())
						create_path(p);
					_index.find(p)->second.type = type_traits<T>::value;
					_index.find(p)->second.is_null = ((_index.find(p)->second.size = s) == 0);
					if (!_index.find(p)->second.is_null) {
						_index.find(p)->second.offset = _mem.size();
						_mem.resize(_mem.size() + s * sizeof(T));
						std::memcpy(&_mem[_index.find(p)->second.offset], v, s * sizeof(T));
					}
				}
				template<typename T> void append_data(std::string const & p, T const * v, std::size_t s) {
					MOCASITO_TRACE
					std::vector<T> buffer(extent(p)[0] + s);
					get_data(p, &buffer.front());
					std::memcpy(&buffer[extent(p)[0]], v, s * sizeof(T));
					set_data(p, &buffer.front(), buffer.size());
				}
				void delete_data(std::string const & p, std::string const & s) {
					MOCASITO_TRACE
					MOCASITO_IO_THROW("not implemented")
				}
				template<typename T> void set_attr(std::string const & p, std::string const & s, T const & v) {
					MOCASITO_TRACE
					if (_index.find(p) == _index.end())
						MOCASITO_IO_THROW("the paht does not exists: " + p)
					std::vector<std::string>::iterator it = std::find(_index.find(p)->second.attr.begin(), _index.find(p)->second.attr.end(), s);
					if (it == _index.find(p)->second.attr.end()) {
						_index.find(p)->second.attr.push_back(s);
						_index.insert(std::make_pair(p + "/@" + s, detail::cache_index(false, false, false, true)));
					}
					_index.find(p + "/@" + s)->second.type = type_traits<T>::value;
					_index.find(p + "/@" + s)->second.offset = _mem.size();
					_mem.resize(_mem.size() + sizeof(T));
					std::memcpy(&_mem[_index.find(p + "/@" + s)->second.offset], &v, sizeof(T));
				}
			private:
				void create_path(std::string const & p) {
					MOCASITO_TRACE
					std::size_t pos;
					for (pos = p.find_last_of('/'); _index.find(p.substr(0, pos)) == _index.end() && pos > 0 && pos < std::string::npos; pos = p.find_last_of('/', pos - 1));
					do {
						_index.find(p.substr(0, pos))->second.children.push_back(p.substr(pos, p.find_first_of('/', pos)));
						if ((pos = p.find_first_of('/', pos + 1)) != std::string::npos)
							_index.insert(std::make_pair(p.substr(0, pos), detail::cache_index(true, false, false, false)));
					} while(pos != std::string::npos);
					_index.find(p.substr(0, p.find_last_of('/')))->second.children.push_back(p.substr(p.find_last_of('/') + 1));
					_index.insert(std::make_pair(p, detail::cache_index(false, false, false, false)));
				}
				std::string _in, _out;
				std::map<std::string, detail::cache_index> _index;
				std::vector<char> _mem;
		};
	}
}
#endif
