// Copyright (C) 2008 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
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
				{}
				context(E * e, std::string const & p)
					: _type(detail::IO_UNKNOWN), _path(p), _engine(e) 
				{
					create(); 
				}
				context(context<E> const & c, std::string const & p)
					: _type(detail::IO_UNKNOWN), _path(p), _engine(c._engine) 
				{
					create(); 
				}
				std::vector<std::size_t> extent() const {
					if (_type == detail::IO_UNKNOWN)
						throw(std::runtime_error("unknown path type: " + _path));
					else if (_type != detail::IO_DATA)
						throw(std::runtime_error("only data nodes have extent: " + _path));
					return _engine->extent(_path); 
				}
				std::size_t dimensions() const {
					if (_type == detail::IO_UNKNOWN)
						throw(std::runtime_error("unknown path type: " + _path));
					else if (_type != detail::IO_DATA)
						throw(std::runtime_error("only data nodes have dimensions: " + _path));
					return _engine->dimensions(_path);
				}
				type_traits<>::type datatype() const {
					switch (_type) {
						case detail::IO_UNKNOWN:
							throw(std::runtime_error("unknown path type: " + _path));
						case detail::IO_GROUP:
							throw(std::runtime_error("groups have no datatype: " + _path));
						case detail::IO_GROUP_ATTR:
						case detail::IO_DATA_ATTR:
							return _engine->attrtype(_type, _path, _segment);
						case detail::IO_DATA:
						default:
							return _engine->datatype(_path);
					}
				}
				std::string path() const { 
					return _path + (_segment.size() ? "@" : "") + _segment; 
				}
				bool is_attribute() const {
					if (_type == detail::IO_UNKNOWN)
						throw(std::runtime_error("unknown path type: " + _path));
					return _type == detail::IO_GROUP_ATTR || _type == detail::IO_DATA_ATTR;
				};
				bool is_leaf() const {
					if (_type == detail::IO_UNKNOWN)
						throw(std::runtime_error("unknown path type: " + _path));
					return _type != detail::IO_GROUP;
				};
				bool is_scalar() const {
					if (_type == detail::IO_UNKNOWN)
						throw(std::runtime_error("unknown path type: " + _path));
					else if (_type != detail::IO_DATA)
						throw(std::runtime_error("only data nodes can be check to be scalar: " + _path));
					return _engine->is_scalar(_path);
				};
				bool exists() const  {
					return (_type != detail::IO_UNKNOWN);
				};
				detail::iterator<E, context<E> > begin() { 
					if (_type != detail::IO_GROUP)
						throw(std::runtime_error("iterators can only loop over groups: " + _path));
					return detail::iterator<E, context<E> >(_engine, _path, _engine->list_children(_path)); 
				}
				detail::iterator<E, context<E> > end() { 
					return detail::iterator<E, context<E> >(); 
				}
				template<typename T> operator T() {
					if (_type == detail::IO_UNKNOWN)
						throw(std::runtime_error("unknown path type: " + _path));
					else if (_type == detail::IO_GROUP)
						throw(std::runtime_error("groups have no data: " + _path));
					T v;
					return assign(v, *this);
				}
				template<typename T> context<E> & operator=(T const & v) {
					return assign(*this, v);
				}
				template<typename T> context<E> & operator<<(T const & v) { 
					return detail::append_helper(*this, v);
				}
				context<E> operator+(std::string const & p) const { 
					std::cerr<<"operator+: this is: "<<this<<std::endl;
					std::cerr<<"operator+: &p is: "<<&p<<std::endl;
					std::cerr<<"strings are: p: "<<p<<" path: "<<_path<<std::endl;
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
					switch (_type) {
						case detail::IO_GROUP:
							throw(std::runtime_error("groups have no data " + _path));
						case detail::IO_DATA:
							_engine->get_data(_path, v); 
							break;
						case detail::IO_GROUP_ATTR:
							_engine->get_group_attr(_path, _segment, *v); 
							break;
						case detail::IO_DATA_ATTR:
							_engine->get_data_attr(_path, _segment, *v); 
							break;
						case detail::IO_UNKNOWN:
						default:
							throw(std::runtime_error("unknown path type " + _path));
					}
				}
				template<typename T> void set(T const & v) {
					switch (_type) {
						case detail::IO_GROUP:
							throw(std::runtime_error("groups have no data " + _path));
						case detail::IO_UNKNOWN:
						case detail::IO_DATA:
							_engine->set_data(_path, v); 
							_type = detail::IO_DATA;
							break;
						case detail::IO_GROUP_ATTR:
							_engine->set_group_attr(_path, _segment, v); 
							break;
						case detail::IO_DATA_ATTR:
							_engine->set_data_attr(_path, _segment, v); 
							break;
					}
				}
				template<typename T> void set(T const * v, std::size_t s) {
					switch (_type) {
						case detail::IO_GROUP:
							throw(std::runtime_error("groups have no data " + _path));
						case detail::IO_DATA:
						case detail::IO_UNKNOWN:
							_engine->set_data(_path, v, s);
							_type = detail::IO_DATA;
							break;
						case detail::IO_GROUP_ATTR:
						case detail::IO_DATA_ATTR:
						default:
							throw(std::runtime_error("attrbuts are scalar data " + _path));
					}
				}
				template<typename T> void append(T const * v, std::size_t s) { 
					if (_type != detail::IO_DATA)
						throw(std::runtime_error("append can only be used for data: " + _path));
					_engine->append_data(_path, v, s);
				}
			private:
				void create() {
					if (_path.find_last_of('@') != std::string::npos) {
						_segment = _path.substr(_path.find_last_of('@'));
						_path = _path.substr(0, _path.find_last_of('@'));
						if (_engine->is_group(_path))
							_type = detail::IO_GROUP_ATTR;
						else if (_engine->is_data(_path))
							_type = detail::IO_DATA_ATTR;
						else
							throw(std::runtime_error("attributes can only be set on existing pathes: " + _path));
					} else if (_engine->is_group(_path))
						_type = detail::IO_GROUP;
					else if (_engine->is_data(_path))
						_type = detail::IO_DATA;
				}
				detail::node_t _type;
				std::string _segment;
				std::string _path;
				E * _engine;
		};
	}
}
#endif
