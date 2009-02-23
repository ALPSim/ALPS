// Copyright (C) 2008 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include "util.hpp"
#include "traits.hpp"
#include <string>
#include <sstream>
#include <vector>
#include <map>
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
					static herr_t callback(unsigned n, H5E_error2_t const * desc, void * buffer) {
						*reinterpret_cast<std::ostringstream *>(buffer) << "    #" << n << " " << desc->file_name << " line " << desc->line << " in " << desc->func_name << "(): " << desc->desc << std::endl;
						return 0;
					}
					static void invoke(MOCASITO_TRACE) { 
						std::ostringstream buffer;
						buffer << "HDR5 trace:" << std::endl;
						H5Ewalk(H5E_DEFAULT, H5E_WALK_DOWNWARD, callback, &buffer);
						MOCASITO_IO_THROW(buffer.str())
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
			typedef h5_t<H5Gclose> h5g_t;
			typedef h5_t<H5Dclose> h5d_t;
			typedef h5_t<H5Aclose> h5a_t;
			typedef h5_t<H5Sclose> h5s_t;
			typedef h5_t<H5Tclose> h5t_t;
			typedef h5_t<H5Pclose> h5p_t;
			typedef h5_t<h5_e::noop> h5e_t;
			template<typename T> class h5_fptr {
				public:
					h5_fptr(std::string const & p) {
						H5Eset_auto(H5E_DEFAULT, NULL, NULL);
						if (_files.find(_name = p) == _files.end()) {
							T id = H5Fopen(p.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
							_files.insert(std::make_pair(p, std::make_pair(id < 0 ? H5Fcreate(p.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT) : id, 1)));
						} else
							_files[p].second++;
					}
					~h5_fptr() {
						if (!--_files[_name].second) {
							H5Fclose(_files[_name].first);
							_files.erase(_name);
						}
					}
					operator T() const {
						return _files[_name].first;
					}
				private:
					std::string _name;
					static std::map<std::string, std::pair<T, std::size_t> > _files;
			};
			template<typename T> std::map<std::string, std::pair<T, std::size_t> > h5_fptr<T>::_files;
		}
		class hdf5 {
			public:
				hdf5(std::string const & p, std::string const & q, MOCASITO_TRACE): _file(p) {
					if (q != "" && p != q)
						MOCASITO_IO_THROW("input file needs to be the same as the output file")
				}
				void flush(MOCASITO_TRACE) const {
					H5Fflush(_file, H5F_SCOPE_GLOBAL);
				}
				bool is_group(std::string const & p, MOCASITO_TRACE) const {
					hid_t id = H5Gopen(_file, p.c_str(), H5P_DEFAULT);
					return id < 0 ? false : static_cast<bool>(detail::h5g_t(id));
				}
				bool is_data(std::string const & p, MOCASITO_TRACE) const {
					hid_t id = H5Dopen(_file, p.c_str(), H5P_DEFAULT);
					return id < 0 ? false : static_cast<bool>(detail::h5d_t(id));
				}
				std::vector<std::size_t> extent(std::string const & p, MOCASITO_TRACE) const {
					if (is_null(p))
						return std::vector<std::size_t>(1, 0);
					std::vector<hsize_t> buffer(dimensions(p), 0);
					{
						detail::h5d_t data_id(H5Dopen(_file, p.c_str(), H5P_DEFAULT));
						detail::h5s_t space_id(H5Dget_space(data_id));
						detail::h5e_t(H5Sget_simple_extent_dims(space_id, &buffer.front(), NULL));
					}
					std::vector<std::size_t> extend(buffer.size(), 0);
					std::copy(buffer.begin(), buffer.end(), extend.begin());
					return extend;
				}
				std::size_t dimensions(std::string const & p, MOCASITO_TRACE) const {
					detail::h5d_t data_id(H5Dopen(_file, p.c_str(), H5P_DEFAULT));
					detail::h5s_t space_id(H5Dget_space(data_id));
					return static_cast<hid_t>(detail::h5e_t(H5Sget_simple_extent_dims(space_id, NULL, NULL)));
				}
				type_traits<>::type attrtype(detail::node_t t, std::string const & p, std::string const & s, MOCASITO_TRACE) const {
					detail::h5d_t data_id(H5Dopen(_file, p.c_str(), 0));
					detail::h5a_t attr_id(H5Aopen(data_id, s.c_str(), H5P_DEFAULT));
					return get_type_id(H5Aget_type(attr_id), p + "/@" + s);
				}
				type_traits<>::type datatype(std::string const & p, MOCASITO_TRACE) const {
					detail::h5d_t data_id(H5Dopen(_file, p.c_str(), 0));
					return get_type_id(H5Dget_type(data_id), p);
				}
				bool is_scalar(std::string const & p, MOCASITO_TRACE) const {
					detail::h5d_t data_id(H5Dopen(_file, p.c_str(), H5P_DEFAULT));
					detail::h5s_t space_id(H5Dget_space(data_id));
					H5S_class_t type = H5Sget_simple_extent_type(space_id);
					if (type == H5S_NO_CLASS)
						MOCASITO_IO_THROW("error reading class " + p)
					return type == H5S_SCALAR;
				}
				bool is_null(std::string const & p, MOCASITO_TRACE) const {
					detail::h5d_t data_id(H5Dopen(_file, p.c_str(), H5P_DEFAULT));
					detail::h5s_t space_id(H5Dget_space(data_id));
					H5S_class_t type = H5Sget_simple_extent_type(space_id);
					if (type == H5S_NO_CLASS)
						MOCASITO_IO_THROW("error reading class " + p)
					return type == H5S_NULL;
				}
				std::vector<std::string> list_children(std::string const & p, MOCASITO_TRACE) const {
					std::vector<std::string> list;
					H5Giterate(_file, p.c_str(), NULL, child_visitor, reinterpret_cast<void *>(&list));
					return list;
				}
				std::vector<std::string> list_attr(std::string const & p, MOCASITO_TRACE) const {
					std::vector<std::string> list;
					if (is_group(p)) {
						detail::h5g_t id(H5Gopen(_file, p.c_str(), H5P_DEFAULT));
						H5Aiterate(id, H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, NULL, attr_visitor, reinterpret_cast<void *>(&list));
					} else {
						detail::h5d_t id(H5Dopen(_file, p.c_str(), H5P_DEFAULT));
						H5Aiterate(id, H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, NULL, attr_visitor, reinterpret_cast<void *>(&list));
					}
					return list;
				}
				template<typename T> void get_data(std::string const & p, T * v, MOCASITO_TRACE) const {
					detail::h5d_t data_id(H5Dopen(_file, p.c_str(), H5P_DEFAULT));
					if (!is_null(p)) {
						detail::h5t_t type_id(get_native_type(v));
						detail::h5e_t(H5Dread(data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, v));
					}
				}
				template<typename T> void get_group_attr(std::string const & p, std::string const & s, T & v, MOCASITO_TRACE) const {
					get_attr<detail::h5g_t, T>(H5Gopen(_file, p.c_str(), H5P_DEFAULT), s, &v);
				}
				template<typename T> void get_data_attr(std::string const & p, std::string const & s, T & v, MOCASITO_TRACE) const {
					get_attr<detail::h5d_t, T>(H5Dopen(_file, p.c_str(), H5P_DEFAULT), s, &v);
				}
				template<typename T> void set_data(std::string const & p, T const & v, MOCASITO_TRACE) {
					detail::h5t_t type_id(get_native_type(v));
					hid_t id = H5Dopen(_file, p.c_str(), H5P_DEFAULT);
					if (id < 0)
						id = create_path(p, type_id, H5Screate(H5S_SCALAR), 0);
					detail::h5d_t data_id(id);
					detail::h5e_t(H5Dwrite(data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v));
				}
				template<typename T> void set_data(std::string const & p, T const * v, hsize_t s, MOCASITO_TRACE) {
					detail::h5t_t type_id(get_native_type(v));
					hid_t id = H5Dopen(_file, p.c_str(), H5P_DEFAULT);
					if (id < 0) {
						id = create_path(p, type_id, s ? H5Screate_simple(1, &s, NULL) : H5Screate(H5S_NULL), s);
					} else 
						detail::h5e_t(H5Dset_extent(id, &s));
					detail::h5d_t data_id(id);
					if (s > 0)
						detail::h5e_t(H5Dwrite(data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, v));
				}
				template<typename T> void append_data(std::string const & p, T const * v, hsize_t s, MOCASITO_TRACE) {
					detail::h5t_t type_id(get_native_type(v));
					hid_t id = H5Dopen(_file, p.c_str(), H5P_DEFAULT);
					if (id < 0)
						return set_data(p, v, s);
					detail::h5d_t data_id(id);
					hsize_t start = extent(p)[0], count = start + s;
					detail::h5e_t(H5Dset_extent(data_id, &count));
					detail::h5s_t space_id(H5Dget_space(data_id));
					detail::h5e_t(H5Sselect_hyperslab(space_id, H5S_SELECT_SET, &start, NULL, &s, NULL));
					detail::h5s_t mem_id(H5Screate_simple(1, &s, NULL));
					detail::h5e_t(H5Dwrite(data_id, type_id, mem_id, space_id, H5P_DEFAULT, v));
				}
				void delete_data(std::string const & p, std::string const & s, MOCASITO_TRACE) {
					detail::h5g_t data_id(H5Dopen(_file, p.c_str(), H5P_DEFAULT));
					detail::h5e_t(H5Ldelete(_file, s.c_str(), data_id));
				}
				template<typename T> void set_group_attr(std::string const & p, std::string const & s, T const & v, MOCASITO_TRACE) {
					set_attr<detail::h5g_t, T>(H5Gopen(_file, p.c_str(), H5P_DEFAULT), s, v);
				}
				template<typename T> void set_data_attr(std::string const & p, std::string const & s, T const & v, MOCASITO_TRACE) {
					set_attr<detail::h5d_t, T>(H5Dopen(_file, p.c_str(), H5P_DEFAULT), s, v);
				}
			private:
				template<typename T> hid_t get_native_type(T &) const { MOCASITO_IO_THROW("unknown type") }
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
				template<typename T> hid_t get_native_type(T * ) const { return get_native_type(T()); }
				static herr_t child_visitor(hid_t id, char const * n, void * d) {
					reinterpret_cast<std::vector<std::string> *>(d)->push_back(n);
					return 0;
				}
				static herr_t attr_visitor(hid_t id, char const * n, const H5A_info_t *, void * d) {
					reinterpret_cast<std::vector<std::string> *>(d)->push_back(n);
					return 0;
				}
				type_traits<>::type get_type_id(detail::h5t_t const & type_id, std::string const & p, MOCASITO_TRACE) const {
					detail::h5t_t native_id(H5Tget_native_type(type_id, H5T_DIR_ASCEND));
					if (false);
					#define MOCASITO_IO_HDF5_GET_TYPE_ID(T)									\
						else if (detail::h5e_t(H5Tequal(detail::h5t_t(H5Tcopy(native_id)), detail::h5t_t(get_native_type(static_cast<T>(0))))) > 0) return type_traits<T>::value;
					MOCASITO_IO_FOREACH_SCALAR(MOCASITO_IO_HDF5_GET_TYPE_ID)
					#undef MOCASITO_IO_HDF5_GET_TYPE_ID
					else MOCASITO_IO_THROW("error comparing types " + p)
				}
				template<typename I, typename T> void get_attr(I const & data_id, std::string const & s, T * v, MOCASITO_TRACE) const {
					detail::h5t_t type_id(get_native_type(v));
					detail::h5a_t attr_id(H5Aopen(data_id, s.c_str(), H5P_DEFAULT));
					detail::h5e_t(H5Aread(attr_id, type_id, v));
				}
				template<typename I, typename T> void set_attr(I const & data_id, std::string const & s, T const & v, MOCASITO_TRACE) {
					detail::h5t_t type_id(get_native_type(v));
					hid_t id = H5Aopen(data_id, s.c_str(), H5P_DEFAULT);
					if (id < 0) {
						detail::h5s_t space_id(H5Screate(H5S_SCALAR));
						id = H5Acreate(data_id, s.c_str(), type_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
					}
					detail::h5a_t attr_id(id);
					detail::h5e_t(H5Awrite(attr_id, type_id, &v));
				}
				hid_t create_path(std::string const & p, hid_t type_id, hid_t space_id, hsize_t s, MOCASITO_TRACE) {
					std::size_t pos;
					hid_t data_id = -1;
					for (pos = p.find_last_of('/'); data_id < 0 && pos > 0 && pos < std::string::npos; pos = p.find_last_of('/', pos - 1))
						data_id = H5Gopen(_file, p.substr(0, pos).c_str(), H5P_DEFAULT);
					if (data_id < 0) {
						pos = p.find_first_of('/', 1);
						detail::h5g_t(H5Gcreate(_file, p.substr(0, pos).c_str(), 0, H5P_DEFAULT, H5P_DEFAULT));
					} else {
						pos = p.find_first_of('/', pos + 1);
						detail::h5g_t(data_id);
					}
					while ((pos = p.find_first_of('/', pos + 1)) != std::string::npos && pos > 0)
						detail::h5g_t(H5Gcreate(_file, p.substr(0, pos).c_str(), 0, H5P_DEFAULT, H5P_DEFAULT));
					detail::h5p_t prop_id(H5Pcreate(H5P_DATASET_CREATE)); 
					detail::h5e_t(H5Pset_fill_time(prop_id, H5D_FILL_TIME_NEVER));
					if (s > 0)
						detail::h5e_t(H5Pset_chunk (prop_id, 1, &s));
					return H5Dcreate(_file, p.c_str(), type_id, detail::h5s_t(space_id), H5P_DEFAULT, prop_id, H5P_DEFAULT);
				}
				detail::h5_fptr<hid_t> _file;
		};
	}
}
#endif
