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
#ifndef IO_HDF5
#define IO_HDF5
namespace mocasito {
	namespace io {
		namespace detail {
			class h5_e {
				public:
					static herr_t noop(hid_t) { return 0; }
					static herr_t callback(unsigned n, const H5E_error2_t *desc, void * buffer) {
						*reinterpret_cast<std::ostringstream *>(buffer) << std::endl << "#" << n << " " << desc->file_name << " line " << desc->line << " in " << desc->func_name << "(): " << desc->desc;
						return 0;
					}
					static void invoke() { 
						std::ostringstream buffer;
						H5Ewalk(H5E_DEFAULT, H5E_WALK_DOWNWARD, callback, &buffer);
						throw(std::runtime_error(buffer.str()));
					}
			};
			template<herr_t(*F)(hid_t)> class h5_t {
				public:
					h5_t(): _id(-1) {}
					h5_t(hid_t id): _id(id) {  if (_id < 0) h5_e::invoke(); H5Eclear(H5E_DEFAULT); }
					~h5_t() { if (_id >= 0 && F(_id) < 0) h5_e::invoke(); H5Eclear(H5E_DEFAULT); }
					operator hid_t() const { return _id; }
					h5_t & operator=(hid_t id) { if ((_id = id) < 0) h5_e::invoke(); H5Eclear(H5E_DEFAULT); return *this; }
				private:
					hid_t _id;
			};
			typedef h5_t<H5Fclose> h5f_t;
			typedef h5_t<H5Gclose> h5g_t;
			typedef h5_t<H5Dclose> h5d_t;
			typedef h5_t<H5Aclose> h5a_t;
			typedef h5_t<H5Sclose> h5s_t;
			typedef h5_t<H5Tclose> h5t_t;
			typedef h5_t<h5_e::noop> h5e_t;
		}
		class hdf5 {
			public:
				hdf5(std::string const & p, std::string const & q) {
					if (q != "" && p != q)
						throw(std::runtime_error("input file needs to be the same as the output file"));
// TODO: pooling
					H5Eset_auto(H5E_DEFAULT, NULL, NULL);
					hid_t id = H5Fopen(p.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
					_file_id = id < 0 ? H5Fcreate(p.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT) : id;
				}
				void flush() const {
					H5Fflush(_file_id, H5F_SCOPE_GLOBAL);
				}
				bool is_group(std::string const & p) const {
					hid_t id = H5Gopen(_file_id, p.c_str(), H5P_DEFAULT);
					return id < 0 ? false : static_cast<bool>(detail::h5g_t(id));
				}
				bool is_data(std::string const & p) const {
					hid_t id = H5Dopen(_file_id, p.c_str(), H5P_DEFAULT);
					return id < 0 ? false : static_cast<bool>(detail::h5d_t(id));
				}
				std::vector<std::size_t> extent(std::string const & p) const {
					if (is_null(p))
						return std::vector<std::size_t>(1, 0);
					std::vector<hsize_t> buffer(dimensions(p));
					{
						detail::h5d_t data_id(H5Dopen(_file_id, p.c_str(), H5P_DEFAULT));
						detail::h5s_t space_id(H5Dget_space(data_id));
						detail::h5e_t(H5Sget_simple_extent_dims(space_id, &buffer.front(), NULL));
					}
					std::vector<std::size_t> extend(buffer.size());
					std::copy(buffer.begin(), buffer.end(), extend.begin());
					return extend;
				}
				std::size_t dimensions(std::string const & p) const {
					detail::h5d_t data_id(H5Dopen(_file_id, p.c_str(), H5P_DEFAULT));
					detail::h5s_t space_id(H5Dget_space(data_id));
					return static_cast<hid_t>(detail::h5e_t(H5Sget_simple_extent_dims(space_id, NULL, NULL)));
				}
				type_traits<>::type attrtype(detail::node_t t, std::string const & p, std::string const & s) const {
// TODO: impl
					return -1;
				}
				type_traits<>::type datatype(std::string const & p) const {
					detail::h5d_t data_id(H5Dopen(_file_id, p.c_str(), 0));
					detail::h5t_t type_id(H5Dget_type(data_id));
					detail::h5t_t native_id(H5Tget_native_type(type_id, H5T_DIR_ASCEND));
					if (detail::h5e_t(H5Tequal(detail::h5t_t(H5Tcopy(native_id)), detail::h5t_t(get_native_type(static_cast<char>(0))))) > 0) return type_traits<char>::value;
					else if (detail::h5e_t(H5Tequal(detail::h5t_t(H5Tcopy(native_id)), detail::h5t_t(get_native_type(static_cast<signed char>(0))))) > 0) return type_traits<signed char>::value;
					else if (detail::h5e_t(H5Tequal(detail::h5t_t(H5Tcopy(native_id)), detail::h5t_t(get_native_type(static_cast<unsigned char>(0))))) > 0) return type_traits<unsigned char>::value;
					else if (detail::h5e_t(H5Tequal(detail::h5t_t(H5Tcopy(native_id)), detail::h5t_t(get_native_type(static_cast<short>(0))))) > 0) return type_traits<short>::value;
					else if (detail::h5e_t(H5Tequal(detail::h5t_t(H5Tcopy(native_id)), detail::h5t_t(get_native_type(static_cast<unsigned char>(0))))) > 0) return type_traits<unsigned char>::value;
					else if (detail::h5e_t(H5Tequal(detail::h5t_t(H5Tcopy(native_id)), detail::h5t_t(get_native_type(static_cast<int>(0))))) > 0) return type_traits<int>::value;
					else if (detail::h5e_t(H5Tequal(detail::h5t_t(H5Tcopy(native_id)), detail::h5t_t(get_native_type(static_cast<unsigned>(0))))) > 0) return type_traits<unsigned>::value;
					else if (detail::h5e_t(H5Tequal(detail::h5t_t(H5Tcopy(native_id)), detail::h5t_t(get_native_type(static_cast<long>(0))))) > 0) return type_traits<long>::value;
					else if (detail::h5e_t(H5Tequal(detail::h5t_t(H5Tcopy(native_id)), detail::h5t_t(get_native_type(static_cast<unsigned long>(0))))) > 0) return type_traits<unsigned long>::value;
					else if (detail::h5e_t(H5Tequal(detail::h5t_t(H5Tcopy(native_id)), detail::h5t_t(get_native_type(static_cast<long long>(0))))) > 0) return type_traits<long long>::value;
					else if (detail::h5e_t(H5Tequal(detail::h5t_t(H5Tcopy(native_id)), detail::h5t_t(get_native_type(static_cast<unsigned long long>(0))))) > 0) return type_traits<unsigned long long>::value;
					else if (detail::h5e_t(H5Tequal(detail::h5t_t(H5Tcopy(native_id)), detail::h5t_t(get_native_type(static_cast<float>(0))))) > 0) return type_traits<float>::value;
					else if (detail::h5e_t(H5Tequal(detail::h5t_t(H5Tcopy(native_id)), detail::h5t_t(get_native_type(static_cast<double>(0))))) > 0) return type_traits<double>::value;
					else if (detail::h5e_t(H5Tequal(detail::h5t_t(H5Tcopy(native_id)), detail::h5t_t(get_native_type(static_cast<long double>(0))))) > 0) return type_traits<long double>::value;
					else if (detail::h5e_t(H5Tequal(detail::h5t_t(H5Tcopy(native_id)), detail::h5t_t(get_native_type(static_cast<bool>(0))))) > 0) return type_traits<bool>::value;
					else throw(std::runtime_error("error comparing types " + p));
				}
				bool is_scalar(std::string const & p) const {
					detail::h5d_t data_id(H5Dopen(_file_id, p.c_str(), H5P_DEFAULT));
					detail::h5s_t space_id(H5Dget_space(data_id));
					H5S_class_t type = H5Sget_simple_extent_type(space_id);
					if (type == H5S_NO_CLASS)
						throw(std::runtime_error("error reading class " + p));
					return type == H5S_SCALAR;
				}
				bool is_null(std::string const & p) const {
					detail::h5d_t data_id(H5Dopen(_file_id, p.c_str(), H5P_DEFAULT));
					detail::h5s_t space_id(H5Dget_space(data_id));
					H5S_class_t type = H5Sget_simple_extent_type(space_id);
					if (type == H5S_NO_CLASS)
						throw(std::runtime_error("error reading class " + p));
					return type == H5S_NULL;
				}
				std::vector<std::string> list_children(std::string const & p) const {
					std::vector<std::string> list;
					H5Giterate(_file_id, p.c_str(), NULL, group_visitor, reinterpret_cast<void *>(&list));
					return list;
				}
				template<typename T> void get_data(std::string const & p, T * v) const {
					detail::h5d_t data_id(H5Dopen(_file_id, p.c_str(), H5P_DEFAULT));
					if (!is_null(p)) {
						detail::h5t_t type_id(get_native_type(v));
						detail::h5e_t(H5Dread(data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, v));
					}
				}
				template<typename T> void get_group_attr(std::string const & p, std::string const & s, T & v) const {
					get_attr<detail::h5g_t, T>(H5Gopen(_file_id, p.c_str(), H5P_DEFAULT), s, &v);
				}
				template<typename T> void get_data_attr(std::string const & p, std::string const & s, T & v) const {
					get_attr<detail::h5d_t, T>(H5Dopen(_file_id, p.c_str(), H5P_DEFAULT), s, &v);
				}
				template<typename T> void set_data(std::string const & p, T const & v) {
					detail::h5t_t type_id(get_native_type(v));
					hid_t id = H5Dopen(_file_id, p.c_str(), H5P_DEFAULT);
					if (id < 0) {
						create_path(p);
						detail::h5s_t space_id(H5Screate(H5S_SCALAR));
						id = H5Dcreate(_file_id, p.c_str(), type_id, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
					}
					detail::h5d_t data_id(id);
					detail::h5e_t(H5Dwrite(data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v));
				}
				template<typename T> void set_data(std::string const & p, T const * v, hsize_t s) {
					detail::h5t_t type_id(get_native_type(v));
					hid_t id = H5Dopen(_file_id, p.c_str(), H5P_DEFAULT);
					if (id < 0) {
						create_path(p);
						if (s == 0) {
							detail::h5s_t space_id(H5Screate(H5S_NULL));
							id = H5Dcreate(_file_id, p.c_str(), type_id, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
						} else {
							detail::h5s_t space_id(H5Screate_simple(1, &s, NULL));
							id = H5Dcreate(_file_id, p.c_str(), type_id, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
						}
					} else 
						detail::h5e_t(H5Dset_extent(id, &s));
					detail::h5d_t data_id(id);
					if (s > 0)
						detail::h5e_t(H5Dwrite(data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, v));
				}
				template<typename T> void append_data(std::string const & p, T const * v, hsize_t s) {
					detail::h5t_t type_id(get_native_type(v));
					hid_t id = H5Dopen(_file_id, p.c_str(), H5P_DEFAULT);
					if (id < 0)
						return set_data(p, v, s);
					detail::h5d_t data_id(id);
					hsize_t start = extent(p)[0], count = start + s;
					detail::h5e_t(H5Dset_extent(data_id, &count));
					detail::h5s_t space_id(H5Dget_space(data_id));
					detail::h5e_t(H5Sselect_hyperslab(space_id, H5S_SELECT_SET, &start, NULL, &s, NULL));
					detail::h5e_t(H5Dwrite(data_id, type_id, H5S_ALL, space_id, H5P_DEFAULT, v));
				}
				void delete_data(std::string const & p, std::string const & s) {
					detail::h5g_t data_id(H5Dopen(_file_id, p.c_str(), H5P_DEFAULT));
					detail::h5e_t(H5Ldelete(_file_id, s.c_str(), data_id));
				}
				template<typename T> void set_group_attr(std::string const & p, std::string const & s, T const & v) {
					set_attr<detail::h5g_t, T>(H5Gopen(_file_id, p.c_str(), H5P_DEFAULT), s, v);
				}
				template<typename T> void set_data_attr(std::string const & p, std::string const & s, T const & v) {
					set_attr<detail::h5d_t, T>(H5Dopen(_file_id, p.c_str(), H5P_DEFAULT), s, v);
				}
			private:
				template<typename T> hid_t get_native_type(T &) const { throw(std::runtime_error("unknown type")); }
				hid_t get_native_type(char) const { return H5Tcopy(H5T_NATIVE_CHAR); }
				hid_t get_native_type(signed char) const { return H5Tcopy(H5T_NATIVE_SCHAR); }
				hid_t get_native_type(unsigned char) const { return H5Tcopy(H5T_NATIVE_UCHAR); }
				hid_t get_native_type(short) const { return H5Tcopy(H5T_NATIVE_SHORT); }
				hid_t get_native_type(unsigned short) const { return H5Tcopy(H5T_NATIVE_USHORT); }
				hid_t get_native_type(int) const { return H5Tcopy(H5T_NATIVE_INT); }
				hid_t get_native_type(unsigned) const { return H5Tcopy(H5T_NATIVE_UINT); }
				hid_t get_native_type(long) const { return H5Tcopy(H5T_NATIVE_LONG); }
				hid_t get_native_type(unsigned long) const { return H5Tcopy(H5T_NATIVE_ULONG); }
				hid_t get_native_type(long long) const { return H5Tcopy(H5T_NATIVE_LLONG); }
				hid_t get_native_type(unsigned long long) const { return H5Tcopy(H5T_NATIVE_ULLONG); }
				hid_t get_native_type(float) const { return H5Tcopy(H5T_NATIVE_FLOAT); }
				hid_t get_native_type(double) const { return H5Tcopy(H5T_NATIVE_DOUBLE); }
				hid_t get_native_type(long double) const { return H5Tcopy(H5T_NATIVE_LDOUBLE); }
				hid_t get_native_type(bool) const { return H5Tcopy(H5T_NATIVE_HBOOL); }
				template<typename T> hid_t get_native_type(T * v) const { return get_native_type(*v); }
				static herr_t group_visitor(hid_t id, char const * n, void * d) {
					reinterpret_cast<std::vector<std::string> *>(d)->push_back(n);
					return 0;
				}
				template<typename I, typename T> void get_attr(I const & data_id, std::string const & s, T * v) const {
					detail::h5t_t type_id(get_native_type(v));
					detail::h5a_t attr_id(H5Aopen(data_id, s.c_str(), H5P_DEFAULT));
					detail::h5e_t(H5Aread(attr_id, type_id, v));
				}
				template<typename I, typename T> void set_attr(I const & data_id, std::string const & s, T const & v) {
					detail::h5t_t type_id(get_native_type(v));
					hid_t id = H5Aopen(data_id, s.c_str(), H5P_DEFAULT);
					if (id < 0) {
						detail::h5s_t space_id(H5Screate(H5S_SCALAR));
						id = H5Acreate(data_id, s.c_str(), type_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
					}
					detail::h5a_t attr_id(id);
					detail::h5e_t(H5Awrite(attr_id, type_id, &v));
				}
				void create_path(std::string const & p) {
					std::size_t pos;
					hid_t data_id = -1;
					for (pos = p.find_last_of('/'); data_id < 0 && pos > 0 && pos < std::string::npos; pos = p.find_last_of('/', pos - 1))
						data_id = H5Gopen(_file_id, p.substr(0, pos).c_str(), H5P_DEFAULT);
					if (data_id < 0) {
						pos = p.find_first_of('/', 1);
						detail::h5g_t(H5Gcreate(_file_id, p.substr(0, pos).c_str(), 0, H5P_DEFAULT, H5P_DEFAULT));
					} else {
						pos = p.find_first_of('/', pos + 1);
						detail::h5g_t(data_id);
					}
					while ((pos = p.find_first_of('/', pos + 1)) != std::string::npos)
						detail::h5g_t(H5Gcreate(_file_id, p.substr(0, pos).c_str(), 0, H5P_DEFAULT, H5P_DEFAULT));
				}
// TODO: needs to be a pool!
				detail::h5f_t _file_id;
		};
	}
}
#endif
