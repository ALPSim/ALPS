// Copyright (C) 2008 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include "util.hpp"
#include <string>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <hdf5.h>
#ifndef IO_CACHE
#define IO_CACHE
namespace mocasito {
	namespace io {
		namespace detail {
			template<typename T1, typename T2, typename T3> class triple {
				public:
					triple(T1 v1, T2 v2, T3 v3)
						: first(v1), second(v2), third(v3)
					{}
					T1 first;
					T2 second;
					T3 third;
			}
		}
		template<typename in_handler, typename out_handler> class cache {
			public:
				cache(std::string const & p, std::string const & q)
					: _in(p), _out(q)
				{
					_in.get(*this);
				}
				void flush() const {
					_out.write(*this);
				}
				bool is_group(std::string const & p) const {
					if (_tree->find(p) == _tree->end() && _index->find(p) == _index->end())
						throw(std::runtime_error("the path does not exists: " + p));
					return _tree->find(p) == _tree->end();
				}
				bool is_data(std::string const & p) const {
					if (_tree->find(p) == _tree->end() && _index->find(p) == _index->end())
						throw(std::runtime_error("the path does not exists: " + p));
					return _index->find(p) == _index->end();
				}
				std::vector<std::size_t> extends(std::string const & p) const {
					if (_index->find(p) == _index->end())
						throw(std::runtime_error("no data path: " + p));
					return _index->find(p)->third;
				}
				std::size_t dimensions(std::string const & p) const {
					if (_index->find(p) == _index->end())
						throw(std::runtime_error("no data path: " + p));
					return _index->find(p)->third.size();
				}
				type_traits<>::type attrtype(detail::node_t t, std::string const & p, std::string const & s) const {
// TODO: impl
throw(std::runtime_error("not Impl"));
					return -1;
				}
				type_traits<>::type datatype(detail::node_t t, std::string const & p) const {
					if (_index->find(p) == _index->end())
						throw(std::runtime_error("no data path: " + p));
					return _index->find(p)->first;
				}
				bool is_scalar(std::string const & p) const {
					if (_index->find(p) == _index->end())
						throw(std::runtime_error("no data path: " + p));
					return _index->find(p)->third.size() == 0;
				}
				std::vector<std::string> list_children(std::string const & p) const {
					if (_tree->find(p) == _tree->end())
						throw(std::runtime_error("no group path: " + p));
					return _tree->find(p)->second;
				}
				template<typename T> void get_data(std::string const & p, T * v) const {
					if (_index->find(p) == _index->end())
						throw(std::runtime_error("no data path: " + p));
					switch (_index->find(p)->first) {
						case type_traits<bool>::value:
							std::memcoyp
						
							template<> struct type_traits<char> { enum { value = 2 }; };
							template<> struct type_traits<signed char> { enum { value = 3 }; };
							template<> struct type_traits<unsigned char> { enum { value = 4 }; };
							template<> struct type_traits<short> { enum { value = 5 }; };
							template<> struct type_traits<unsigned short> { enum { value = 6 }; };
							template<> struct type_traits<int> { enum { value = 7 }; };
							template<> struct type_traits<unsigned> { enum { value = 8 }; };
							template<> struct type_traits<long> { enum { value = 9 }; };
							template<> struct type_traits<long long> { enum { value = 10 }; };
							template<> struct type_traits<unsigned long long> { enum { value = 11 }; };
							template<> struct type_traits<float> { enum { value = 12 }; };
							template<> struct type_traits<double> { enum { value = 13 }; };
							template<> struct type_traits<long double> { enum { value = 14 }; };
						default:
							throw(std::runtime_error("unknown type: " + _index->find(p)->first));
					}
					
					
					
					
						
						
						
					v = static_cast<T>(*reinterprete_cast<long long *>(&_data->front() + _index->find(_path)->second));
						
						
						
					return _index->find(p)->third.size() == 0;

detail::h5d_t data_id(H5Dopen(_file_id, p.c_str(), H5P_DEFAULT));
					detail::h5t_t type_id(get_native_type(v));
					detail::h5e_t(H5Dread(data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, v));
				}
				template<typename T> void get_group_attr(std::string const & p, std::string const & s, T & v) const {
// TODO: impl
throw(std::runtime_error("not Impl"));
				}
				template<typename T> void get_data_attr(std::string const & p, std::string const & s, T & v) const {
// TODO: impl
throw(std::runtime_error("not Impl"));
				}
				template<typename T> void set_data(std::string const & p, T const & v) {
// TODO: impl
throw(std::runtime_error("not Impl"));
				}
				template<typename T> void set_data(std::string const & p, T const * v, hsize_t s) {
// TODO: impl
throw(std::runtime_error("not Impl"));
				}
				template<typename T> void append_data(std::string const & p, T const * v, hsize_t s) {
// TODO: impl
throw(std::runtime_error("not Impl"));
				}
				template<typename T> void set_group_attr(std::string const & p, std::string const & s, T const & v) {
// TODO: impl
throw(std::runtime_error("not Impl"));
				}
				template<typename T> void set_data_attr(std::string const & p, std::string const & s, T const & v) {
// TODO: impl
throw(std::runtime_error("not Impl"));
				}
			private:
				in_handler _in;
				out_handler _out;
				std::map<std::string, std::vector<std::string> > _tree;
				std::map<std::string, detail::triple<std::size_t, std::size_t, std::vector<std::size_t> > _index;
				std::vector<char> _data;
	}
}
#endif