// Copyright (C) 2008 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include "../trace.hpp"
#include "util.hpp"
#include "traits.hpp"
#include "iterator.hpp"
#include <string>
#include <vector>
#include <cstring>
#include <stdexcept>
#ifndef IO_CONTEXT
#define IO_CONTEXT
namespace mocasito {
	namespace io {
		namespace detail {
			template<typename C, typename T> C & append_helper(C & c, T const & v) {
				return append(c, v);
			}
		}
		template<typename E> class context {
			public:
				context() {}
				context(context<E> const & c)
					: _type(c._type), _segment(c._segment), _path(c._path), _engine(c._engine) 
				{
					MOCASITO_TRACE
				}
				context(E * e, std::string const & p)
					: _type(detail::IO_UNKNOWN), _path(p), _engine(e) 
				{
					MOCASITO_TRACE
					create(); 
				}
				context(context<E> const & c, std::string const & p)
					: _type(detail::IO_UNKNOWN), _path(p), _engine(c._engine) 
				{
					MOCASITO_TRACE
					create(); 
				}
				std::vector<std::size_t> extent() const {
					MOCASITO_TRACE
					if (_type == detail::IO_UNKNOWN)
						MOCASITO_IO_THROW("unknown path type: " + _path)
					else if (_type != detail::IO_DATA)
						MOCASITO_IO_THROW("only data nodes have extent: " + _path)
					return _engine->extent(_path); 
				}
				std::size_t dimensions() const {
					MOCASITO_TRACE
					if (_type == detail::IO_UNKNOWN)
						MOCASITO_IO_THROW("unknown path type: " + _path)
					else if (_type != detail::IO_DATA)
						MOCASITO_IO_THROW("only data nodes have dimensions: " + _path)
					return _engine->dimensions(_path);
				}
				type_traits<>::type datatype() const {
					MOCASITO_TRACE
					switch (_type) {
						case detail::IO_UNKNOWN:
							MOCASITO_IO_THROW("unknown path type: " + _path)
						case detail::IO_GROUP:
							MOCASITO_IO_THROW("groups have no datatype: " + _path)
						case detail::IO_ATTR:
							return _engine->attrtype(_path, _segment);
						case detail::IO_DATA:
						default:
							return _engine->datatype(_path);
					}
				}
				std::string path() const { 
					MOCASITO_TRACE
					return _path + (_segment.size() ? "@" : "") + _segment; 
				}
				bool is_attribute() const {
					MOCASITO_TRACE
					if (_type == detail::IO_UNKNOWN)
						MOCASITO_IO_THROW("unknown path type: " + _path)
					return _type == detail::IO_ATTR;
				};
				bool is_leaf() const {
					MOCASITO_TRACE
					if (_type == detail::IO_UNKNOWN)
						MOCASITO_IO_THROW("unknown path type: " + _path)
					return _type != detail::IO_GROUP;
				};
				bool is_scalar() const {
					MOCASITO_TRACE
					if (_type == detail::IO_UNKNOWN)
						MOCASITO_IO_THROW("unknown path type: " + _path)
					else if (_type != detail::IO_DATA)
						MOCASITO_IO_THROW("only data nodes can be check to be scalar: " + _path)
					return _engine->is_scalar(_path);
				};
				bool exists() const  {
					MOCASITO_TRACE
					return (_type != detail::IO_UNKNOWN);
				};
				detail::iterator<E, context<E> > begin() { 
					MOCASITO_TRACE
					if (_type != detail::IO_GROUP)
						MOCASITO_IO_THROW("iterators can only loop over groups: " + _path)
					return detail::iterator<E, context<E> >(_engine, _path, _engine->list_children(_path)); 
				}
				detail::iterator<E, context<E> > end() { 
					MOCASITO_TRACE
					return detail::iterator<E, context<E> >(); 
				}
				template<typename T> operator T() const {
					MOCASITO_TRACE
					if (_type == detail::IO_UNKNOWN)
						MOCASITO_IO_THROW("unknown path type: " + _path)
					else if (_type == detail::IO_GROUP)
						MOCASITO_IO_THROW("groups have no data: " + _path)
					T v;
					return assign(v, *this);
				}
				template<typename T> context<E> & operator=(T const & v) {
					MOCASITO_TRACE
					return assign(*this, v);
				}
				template<typename T> context<E> & operator<<(T const & v) {
					MOCASITO_TRACE
					return detail::append_helper(*this, v);
				}
				context<E> operator+(char const * p) const {
					MOCASITO_TRACE
					this->operator+(std::string(p));
				}
				context<E> operator+(std::string const & p) const {
					MOCASITO_TRACE
					if (p[0] == '/')
						return context<E>(_engine, p);
					else if (p.substr(0, 2) != "..")
						return context<E>(_engine, _path + "/" + p);
					else {
						std::string q = _path;
						for (std::size_t i = 0; p.substr(i, 2) == ".."; i += 3)
							q = q.substr(0, q.find_last_of('/'));
						return context<E>(_engine, q + "/" + p);
					}
				}
				template<typename T> void get(T * v) const { 
					MOCASITO_TRACE
					switch (_type) {
						case detail::IO_GROUP:
							MOCASITO_IO_THROW("groups have no data " + _path)
						case detail::IO_DATA:
							_engine->get_data(_path, v); 
							break;
						case detail::IO_ATTR:
							_engine->get_attr(_path, _segment, *v); 
							break;
						case detail::IO_UNKNOWN:
						default:
							MOCASITO_IO_THROW("unknown path type " + _path)
					}
				}
				template<typename T> void set(T const & v) {
					MOCASITO_TRACE
					switch (_type) {
						case detail::IO_GROUP:
							MOCASITO_IO_THROW("groups have no data " + _path)
						case detail::IO_UNKNOWN:
						case detail::IO_DATA:
							_engine->set_data(_path, v); 
							_type = detail::IO_DATA;
							break;
						case detail::IO_ATTR:
							_engine->set_attr(_path, _segment, v); 
							break;
					}
				}
				template<typename T> void set(T const * v, std::size_t s) {
					MOCASITO_TRACE
					switch (_type) {
						case detail::IO_GROUP:
							MOCASITO_IO_THROW("groups have no data " + _path)
						case detail::IO_DATA:
						case detail::IO_UNKNOWN:
							_engine->set_data(_path, v, s);
							_type = detail::IO_DATA;
							break;
						case detail::IO_ATTR:
						default:
							MOCASITO_IO_THROW("attrbuts are scalar data " + _path)
					}
				}
				template<typename T> void append(T const * v, std::size_t s) { 
					MOCASITO_TRACE
					if (_type != detail::IO_DATA)
						MOCASITO_IO_THROW("append can only be used for data: " + _path)
					_engine->append_data(_path, v, s);
				}
			private:
				void create() {
					MOCASITO_TRACE
					if (_path.find_last_of('@') != std::string::npos) {
						_segment = _path.substr(_path.find_last_of('@') + 1);
						_path = _path.substr(0, _path.find_last_of('@') - 1);
						_type = detail::IO_ATTR;
					} else if (_engine->is_group(_path))
						_type = detail::IO_GROUP;
					else if (_engine->is_data(_path))
						_type = detail::IO_DATA;
					if (*_path.rbegin() == '/')
						_path = _path.substr(0, _path.size() - 1);
				}
				detail::node_t _type;
				std::string _segment;
				std::string _path;
				E * _engine;
		};
		template<typename T, typename E> bool operator==(context<E> const & c, T const & v) { MOCASITO_TRACE T w; assign(w, c); return w == v; }
		template<typename E> bool operator==(context<E> const & c, char const * v) { MOCASITO_TRACE return c == std::string(v); }
		template<typename T, typename E> bool operator==(T const & v, context<E> const & c) { MOCASITO_TRACE return c == v; }
		template<typename T, typename E> bool operator!=(context<E> const & c, T const & v) { MOCASITO_TRACE return !(c == v); }
		template<typename T, typename E> bool operator!=(T const & v, context<E> const & c) { MOCASITO_TRACE return !(v == c); }
	}
}
#endif
