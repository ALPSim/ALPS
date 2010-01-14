// Copyright (C) 2008 - 2010 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef ALPS_HDF5_HPP
#define ALPS_HDF5_HPP

#ifndef _HDF5USEDLL_
# define _HDF5USEDLL_
#endif
#ifndef _HDF5USEHLDLL_
# define _HDF5USEHLDLL_
#endif

#include <string>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <complex>
#include <stdexcept>
#include <valarray>
#include <iostream>

#include <boost/config.hpp>
#include <boost/utility.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/if.hpp>

#include <hdf5.h>

namespace alps {
	namespace hdf5 {
		namespace detail {
			#define HDF5_ADD_CV(callback, T)														\
				callback(T)																			\
				callback(T &)																		\
				callback(T const)																	\
				callback(T const &)																	\
				callback(T volatile)																\
				callback(T volatile &)																\
				callback(T const volatile)															\
				callback(T const volatile &)
			#define HDF5_FOREACH_SCALAR(callback)													\
				callback(char)																		\
				callback(signed char)																\
				callback(unsigned char)																\
				callback(short)																		\
				callback(unsigned short)															\
				callback(int)																		\
				callback(unsigned int)																\
				callback(long)																		\
				callback(unsigned long)																\
				callback(long long)																	\
				callback(unsigned long long)														\
				callback(float)																		\
				callback(double)																	\
				callback(long double)
			struct scalar_tag {};
			struct stl_complex_tag {};
			struct stl_pair_tag {};
			struct stl_container_of_unknown_tag {};
			struct stl_container_of_string_tag {};
			struct stl_container_of_scalar_tag {};
			struct stl_container_of_container_tag {};
			struct stl_string_tag {};
			struct c_string_tag {};
			struct enum_tag {};
			template<typename T> struct is_writable : boost::mpl::or_<typename boost::is_scalar<T>::type, typename boost::is_enum<T>::type>::type { 
				typedef typename boost::mpl::if_<typename boost::is_enum<T>::type, enum_tag, scalar_tag>::type category; 
			};
			template<typename T> struct is_writable<std::complex<T> > : boost::mpl::true_ { typedef stl_complex_tag category; };
			template<typename T, typename U> struct is_writable<std::pair<T, U> > : boost::mpl::true_ { typedef stl_pair_tag category; };
			template<typename T> struct is_writable<std::vector<T> > : boost::mpl::true_ { typedef stl_container_of_unknown_tag category; };
			template<typename T> struct is_writable<std::valarray<T> > : boost::mpl::true_ { typedef stl_container_of_unknown_tag category; };
			template<typename T> struct is_writable<std::deque<T> > : boost::mpl::true_ { typedef stl_container_of_unknown_tag category; };
			template<typename T> struct is_writable<std::list<T> > : boost::mpl::true_ { typedef stl_container_of_unknown_tag category; };
			template<> struct is_writable<std::vector<std::string> > : boost::mpl::true_ { typedef stl_container_of_string_tag category; };
			template<> struct is_writable<std::valarray<std::string> > : boost::mpl::true_ { typedef stl_container_of_string_tag category; };
			#define HDF5_CONTAINER_OF_SCALAR_CV(T)																						\
				template<> struct is_writable<std::vector<T> > : boost::mpl::true_ { typedef stl_container_of_scalar_tag category; };	\
				template<> struct is_writable<std::valarray<T> > : boost::mpl::true_ { typedef stl_container_of_scalar_tag category; };
			#define HDF5_CONTAINER_OF_SCALAR(T)																							\
				HDF5_ADD_CV(HDF5_CONTAINER_OF_SCALAR_CV, T)
			HDF5_FOREACH_SCALAR(HDF5_CONTAINER_OF_SCALAR)
			#undef HDF5_CONTAINER_OF_SCALAR
			template<typename T> struct is_writable<std::vector<std::valarray<T> > > : boost::mpl::true_ { typedef stl_container_of_container_tag category; };
			template<typename T, typename C, typename A> struct is_writable<std::set<T, C, A> > : boost::mpl::true_ { typedef stl_container_of_unknown_tag category; };
			template<typename K, typename D, typename C, typename A> struct is_writable<std::map<K, D, C, A> > : boost::mpl::true_ { typedef stl_container_of_unknown_tag category; };
			template<> struct is_writable<std::string> : boost::mpl::true_ { typedef stl_string_tag category; };
			template<> struct is_writable<char const *> : boost::mpl::true_ { typedef c_string_tag category; };
			template<std::size_t N> struct is_writable<const char [N]> : boost::mpl::true_ { typedef c_string_tag category; };
			template<std::size_t N> struct is_writable<char [N]> : boost::mpl::true_ { typedef c_string_tag category; };
			class error {
				public:
					static herr_t noop(hid_t) { return 0; }
					static herr_t callback(unsigned n, H5E_error2_t const * desc, void * buffer) {
						*reinterpret_cast<std::ostringstream *>(buffer) << "    #" << n << " " << desc->file_name << " line " << desc->line << " in " << desc->func_name << "(): " << desc->desc << std::endl;
						return 0;
					}
					static std::string invoke() {
						std::ostringstream buffer;
						buffer << "HDF5 error:" << std::endl;
						H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, callback, &buffer);
						return buffer.str();
					}
			};
			template<herr_t(*F)(hid_t)> class ressource {
				public:
					ressource(): _id(-1) {}
					ressource(hid_t id): _id(id) {
						if (_id < 0)
							throw std::runtime_error(error::invoke()); 
						H5Eclear2(H5E_DEFAULT);
					}
					~ressource() {
						if(_id < 0 || F(_id) < 0) {
							std::cerr << error::invoke() << std::endl;
							std::abort();
						}
						H5Eclear2(H5E_DEFAULT); 
					}
					operator hid_t() const { 
						return _id; 
					}
					ressource<F> & operator=(hid_t id) { 
						if ((_id = id) < 0) 
							throw std::runtime_error(error::invoke()); 
						H5Eclear2(H5E_DEFAULT); 
						return *this; 
					}
				private:
					hid_t _id;
			};
			typedef ressource<H5Fclose> file_type;
			typedef ressource<H5Gclose> group_type;
			typedef ressource<H5Dclose> data_type;
			typedef ressource<H5Aclose> attribute_type;
			typedef ressource<H5Sclose> space_type;
			typedef ressource<H5Tclose> type_type;
			typedef ressource<H5Pclose> property_type;
			typedef ressource<error::noop> error_type;
			template <typename T> T check_file(T id) { file_type unused(id); return unused; }
			template <typename T> T check_group(T id) { group_type unused(id); return unused; }
			template <typename T> T check_data(T id) { data_type unused(id); return unused; }
			template <typename T> T check_attribute(T id) { attribute_type unused(id); return unused; }
			template <typename T> T check_space(T id) { space_type unused(id); return unused; }
			template <typename T> T check_type(T id) { type_type unused(id); return unused; }
			template <typename T> T check_property(T id) { property_type unused(id); return unused; }
			template <typename T> T check_error(T id) { error_type unused(id); return unused; }
			struct write {};
			struct read {};
			template <typename Tag> class archive: boost::noncopyable {
				public:
					archive(std::string const & file) {
						H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
						if (boost::is_same<Tag, write>::value) {
							if (H5Fis_hdf5(file.c_str()) == 0)
								throw std::runtime_error("no valid hdf5 file " + file);
							hid_t id = H5Fopen(file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
							_file = (id < 0 ? H5Fcreate(file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT) : id);
						} else {
							if (check_error(H5Fis_hdf5(file.c_str())) == 0)
								throw std::runtime_error("no valid hdf5 file " + file);
							_file = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
						}
					}
					~archive() {
						H5Fflush(_file, H5F_SCOPE_GLOBAL);
						if (H5Fget_obj_count(_file, H5F_OBJ_ALL) - H5Fget_obj_count(_file, H5F_OBJ_FILE) > 0) {
							std::cerr << "Not all ressources closed" << std::endl;
							std::abort();
						}
					}
					std::string get_context() const {
						return _context;
					}
					std::string compute_path(std::string const & p) const {
						if (p.size() && p[0] == '/')
							return p;
						else if (p.size() < 2 || p.substr(0, 2) != "..")
							return _context + (_context.size() == 1 || !p.size() ? "" : "/") + p;
						else {
							std::string s = _context;
							std::size_t i = 0;
							for (; s.size() && p.substr(i, 2) == ".."; i += 3)
								s = s.substr(0, s.find_last_of('/'));
							return s + (s.size() == 1 || !p.substr(i).size() ? "" : "/") + p.substr(i);
						}
					}
					template<typename T> typename boost::enable_if<
						typename boost::mpl::and_<is_writable<T>, typename boost::is_same<Tag, write>::type >
					>::type serialize(std::string const & p, T const & v) {
						if (p.find_last_of('@') != std::string::npos)
							set_attr(compute_path(p).substr(0, compute_path(p).find_last_of('@') - 1), p.substr(p.find_last_of('@') + 1), v, typename is_writable<T>::category());
						else
							set_data(compute_path(p), v, typename is_writable<T>::category());
					}
					template<typename T> typename boost::enable_if<
						typename boost::mpl::and_<is_writable<T>, typename boost::is_same<Tag, read>::type >
					>::type serialize(std::string const & p, T & v) {
						if (p.find_last_of('@') != std::string::npos)
							get_attr(compute_path(p.substr(0, p.find_last_of('@') - 1)), p.substr(p.find_last_of('@') + 1), v, typename is_writable<T>::category());
						else
							get_data(compute_path(p), v, typename is_writable<T>::category());
					}
					template<typename T> typename boost::disable_if<is_writable<T> >::type serialize(std::string const & p, T & v) {
						std::string c = _context;
						_context = compute_path(p);
						v.serialize(*this);
						_context = c;
					}
					bool is_group(std::string const & p) const {
						hid_t id = H5Gopen2(_file, compute_path(p).c_str(), H5P_DEFAULT);
						return id < 0 ? false : check_group(id) != 0;
					}
					bool is_data(std::string const & p) const {
						hid_t id = H5Dopen2(_file, compute_path(p).c_str(), H5P_DEFAULT);
						return id < 0 ? false : check_data(id) != 0;
					}
					std::vector<std::size_t> extent(std::string const & p) const {
						if (is_null(p))
							return std::vector<std::size_t>(1, 0);
						std::vector<hsize_t> buffer(dimensions(p), 0);
						{
							data_type data_id(H5Dopen2(_file, compute_path(p).c_str(), H5P_DEFAULT));
							space_type space_id(H5Dget_space(data_id));
							check_error(H5Sget_simple_extent_dims(space_id, &buffer.front(), NULL));
						}
						std::vector<std::size_t> extend(buffer.size(), 0);
						std::copy(buffer.begin(), buffer.end(), extend.begin());
						return extend;
					}
					std::size_t dimensions(std::string const & p) const {
						data_type data_id(H5Dopen2(_file, compute_path(p).c_str(), H5P_DEFAULT));
						space_type space_id(H5Dget_space(data_id));
						return static_cast<hid_t>(check_error(H5Sget_simple_extent_dims(space_id, NULL, NULL)));
					}
					bool is_scalar(std::string const & p) const {
						data_type data_id(H5Dopen2(_file, compute_path(p).c_str(), H5P_DEFAULT));
						space_type space_id(H5Dget_space(data_id));
						H5S_class_t type = H5Sget_simple_extent_type(space_id);
						if (type == H5S_NO_CLASS)
							throw std::runtime_error("error reading class " + compute_path(p));
						return type == H5S_SCALAR;
					}
					bool is_null(std::string const & p) const {
						data_type data_id(H5Dopen2(_file, compute_path(p).c_str(), H5P_DEFAULT));
						space_type space_id(H5Dget_space(data_id));
						H5S_class_t type = H5Sget_simple_extent_type(space_id);
						if (type == H5S_NO_CLASS)
							throw std::runtime_error("error reading class " + compute_path(p));
						return type == H5S_NULL;
					}
					std::vector<std::string> list_children(std::string const & p) const {
						std::vector<std::string> list;
						group_type group_id(H5Gopen2(_file, compute_path(p).c_str(), H5P_DEFAULT));
						check_error(H5Literate(group_id, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, child_visitor, reinterpret_cast<void *>(&list)));
						return list;
					}
					std::vector<std::string> list_attr(std::string const & p) const {
						std::vector<std::string> list;
						if (is_group(p)) {
							group_type id(H5Gopen2(_file, compute_path(p).c_str(), H5P_DEFAULT));
							check_error(H5Aiterate2(id, H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, NULL, attr_visitor, reinterpret_cast<void *>(&list)));
						} else {
							data_type id(H5Dopen2(_file, compute_path(p).c_str(), H5P_DEFAULT));
							check_error(H5Aiterate2(id, H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, NULL, attr_visitor, reinterpret_cast<void *>(&list)));
						}
						return list;
					}
					void set_group(std::string const & p) const {
						std::string q = compute_path(p);
						if (!is_group(q)) {
							std::size_t pos;
							hid_t group_id = -1;
							for (pos = q.find_last_of('/'); group_id < 0 && pos > 0 && pos < std::string::npos; pos = q.find_last_of('/', pos - 1))
								group_id = H5Gopen2(_file, q.substr(0, pos).c_str(), H5P_DEFAULT);
							if (group_id < 0) {
								if ((pos = q.find_first_of('/', 1)) != std::string::npos)
									check_group(H5Gcreate2(_file, q.substr(0, pos).c_str(), 0, H5P_DEFAULT, H5P_DEFAULT));
							} else {
								pos = q.find_first_of('/', pos + 1);
								check_group(group_id);
							}
							while (pos != std::string::npos && (pos = q.find_first_of('/', pos + 1)) != std::string::npos && pos > 0)
								check_group(H5Gcreate2(_file, q.substr(0, pos).c_str(), 0, H5P_DEFAULT, H5P_DEFAULT));
							check_group(H5Gcreate2(_file, q.c_str(), 0, H5P_DEFAULT, H5P_DEFAULT));
						}
					}
					template<typename T> typename boost::enable_if<boost::is_scalar<T> >::type get_data(std::string const & p, T * v) const {
						if (!is_null(p)) {
							data_type data_id(H5Dopen2(_file, compute_path(p).c_str(), H5P_DEFAULT));
							type_type type_id(get_native_type(v));
							check_error(H5Dread(data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, v));
						}
					}
					template<typename T> typename boost::enable_if<boost::is_scalar<T> >::type get_data(std::string const & p, std::complex<T> * v) const {
						get_data(p, reinterpret_cast<T *>(v));
					}
					template<typename T> typename boost::enable_if<boost::is_scalar<T> >::type set_data(std::string const & p, T const * v, hsize_t s) const {
						type_type type_id(get_native_type(v));
						data_type data_id(create_path(compute_path(p), type_id, s ? H5Screate_simple(1, &s, NULL) : H5Screate(H5S_NULL), s ? 1 : 0, &s));
						if (s > 0)
							check_error(H5Dwrite(data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, v));
					}
					template<typename T> typename boost::enable_if<boost::is_scalar<T> >::type set_data(std::string const & p, std::complex<T> const * v, hsize_t s) const {
						set_data(p, reinterpret_cast<T const *>(v), 2 * s);
					}
				private:
					template<typename T> hid_t get_native_type(T &) const { throw std::runtime_error("unknown type"); }
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
					static herr_t child_visitor(hid_t, char const * n, const H5L_info_t *, void * d) {
						reinterpret_cast<std::vector<std::string> *>(d)->push_back(n);
						return 0;
					}
					static herr_t attr_visitor(hid_t, char const * n, const H5A_info_t *, void * d) {
						reinterpret_cast<std::vector<std::string> *>(d)->push_back(n);
						return 0;
					}
					hid_t create_path(std::string const & p, hid_t type_id, hid_t space_id, hsize_t d, hsize_t const * s = NULL, bool set_prop = true) const {
						hid_t id = H5Dopen2(_file, p.c_str(), H5P_DEFAULT);
						if (d > 0 && s[0] > 0 && id >= 0 && is_null(p)) {
							check_data(id);
							check_error(H5Ldelete(_file, p.c_str(), H5P_DEFAULT));
						} else if (id < 0) {
							if (p.find_last_of('/') < std::string::npos && p.find_last_of('/') > 0)
								set_group(p.substr(0, p.find_last_of('/')));
						} else if (d > 0 && s[0] > 0) {
							check_error(H5Dset_extent(id, s));
							return id;
						} else
							return id;
						if (set_prop) {
							property_type prop_id(H5Pcreate(H5P_DATASET_CREATE));
							check_error(H5Pset_fill_time(prop_id, H5D_FILL_TIME_NEVER));
							if (d > 0)
								check_error(H5Pset_chunk(prop_id, d, s));
							return H5Dcreate2(_file, p.c_str(), type_id, space_type(space_id), H5P_DEFAULT, prop_id, H5P_DEFAULT);
						} else
							return H5Dcreate2(_file, p.c_str(), type_id, space_type(space_id), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
					}
					template<typename T> void get_data(std::string const & p, T & v, scalar_tag) const {
						if (is_null(p))
							throw std::runtime_error("the path '" + p + "' is null");
						if (!is_scalar(p))
							throw std::runtime_error("the path '" + p + "' is not a scalar");
						get_data(p, &v);
					}
					template<typename T> void get_data(std::string const & p, T & v, stl_pair_tag) const {
						const_cast<archive<Tag> &>(*this) >> make_pvp(p + "/first", v.first);
						const_cast<archive<Tag> &>(*this) >> make_pvp(p + "/second", v.second);
					}
					template<typename T> void get_data(std::string const & p, T & v, stl_complex_tag) const {
						get_data(p, static_cast<typename T::value_type *>(&v));
					}
					template<typename T> void get_data(std::string const & p, T & v, stl_container_of_scalar_tag) const {
						if (is_null(p))
							v.resize(0);
						else {
							if (dimensions(p) != 1)
								throw std::runtime_error("the path " + p + " has not dimension 1");
							v.resize(extent(p)[0]);
							get_data(p, &(v[0]));
						}
					}
					template<typename T> void get_data(std::string const & p, T & v, stl_string_tag) const {
						if (is_null(p))
							v = "";
						else if (is_scalar(p)) {
							data_type data_id(H5Dopen2(_file, p.c_str(), H5P_DEFAULT));
							type_type type_id(H5Dget_type(data_id));
							type_type native_id(H5Tget_native_type(type_id, H5T_DIR_ASCEND));
							if (H5Tget_class(native_id) == H5T_STRING) {
								v.resize(H5Tget_size(native_id));
								check_error(H5Dread(data_id, type_type(H5Tcopy(type_id)), H5S_ALL, H5S_ALL, H5P_DEFAULT, &v[0]));
							}
							#define HDF5_GET_STRING(T)																		\
								else if (check_error(H5Tequal(																\
									check_type(H5Tcopy(native_id)), check_type(get_native_type<T>(0))						\
								)) > 0) {																					\
									T t;																					\
									check_error(H5Dread(																	\
										data_id, check_type(H5Tcopy(type_id)), H5S_ALL, H5S_ALL, H5P_DEFAULT, &t			\
									));																						\
									v = boost::lexical_cast<std::string>(t);												\
								}
							HDF5_FOREACH_SCALAR(HDF5_GET_STRING)
							#undef HDF5_GET_STRING
							else throw std::runtime_error("error in types: " + p);
						} else
							throw std::runtime_error("vectors cannot be read into strings: " + p);
					}
					template<typename T> void get_data(std::string const & p, T & v, stl_container_of_unknown_tag) const {
						v.clear();
						if (is_group(p)) {
							std::size_t size;
							const_cast<archive<Tag> *>(this)->serialize(p + "/@length", size);
							for (std::size_t i = 0; i <  size; ++i) {
								typename T::value_type t;
								const_cast<archive<Tag> &>(*this) >> make_pvp(p + "/" + boost::lexical_cast<std::string>(i), t);
								v.push_back(t);
							}
						}
					}
					template<typename K, typename D, typename C, typename A> void get_data(std::string const & p, std::map<K, D, C, A> & v, stl_container_of_unknown_tag) const {
						v.clear();
						if (is_group(p)) {
							std::size_t size;
							const_cast<archive<Tag> *>(this)->serialize(p + "/@length", size);
							for (std::size_t i = 0; i <  size; ++i) {
								std::pair<K, D> t;
								const_cast<archive<Tag> &>(*this) >> make_pvp(p + "/" + boost::lexical_cast<std::string>(i), t);
								v.insert(t);
							}
						}
					}
					template<typename T> void get_data(std::string const & p, T & v, stl_container_of_string_tag) const {
						if (is_null(p))
							v.clear();
						else {
							if (dimensions(p) != 1)
								throw std::runtime_error("the path " + p + " has not dimension 1");
							v.resize(extent(p)[0]);
							data_type data_id(H5Dopen2(_file, p.c_str(), 0));
							type_type type_id(H5Dget_type(data_id));
							type_type native_id(H5Tget_native_type(type_id, H5T_DIR_ASCEND));
							if (H5Tget_class(native_id) != H5T_STRING)
								throw std::runtime_error("the path " + p + " does not contain strings");
							char **data = new char*[v.size()];
							check_error(H5Dread(data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data));
							for (std::size_t i = 0; i < v.size(); ++i) {
								v[i] = data[i];
								free(data[i]);
							}
							delete[] data;
						}
					}
					template<typename T> void get_data(std::string const & p, T & v, stl_container_of_container_tag) const {
						if (is_null(p))
							v.clear();
						else {
							if (dimensions(p) != 2)
								throw std::runtime_error("the path " + p + " has not dimension 2");
							v.resize(extent(p)[0]);
							v[0].resize(extent(p)[1]);
							for (typename T::iterator it = v.begin() + 1; it != v.end(); ++it)
								it->resize(v[0].size());
							data_type data_id(H5Dopen2(_file, p.c_str(), H5P_DEFAULT));
							type_type type_id(get_native_type(const_cast<T &>(v)[0][0]));
							for (std::size_t i = 0; i < v.size(); ++i) {
								hsize_t start[2] = { i, 0 }, count[2] = { 1, v[i].size() };
								space_type space_id(H5Dget_space(data_id));
								check_error(H5Sselect_hyperslab(space_id, H5S_SELECT_SET, start, NULL, count, NULL));
								space_type mem_id(H5Screate_simple(2, count, NULL));
								check_error(H5Dread(data_id, type_id, mem_id, space_id, H5P_DEFAULT, &(const_cast<T &>(v)[i][0])));
							}
						}
					}
					template<typename T> void get_data(std::string const & p, T & v, c_string_tag) const {
						std::string s;
						get_data(p, s, stl_string_tag());
						std::strcpy(v, s.c_str());
					}
					template<typename T> void get_attr(std::string const &, std::string const & p, T &, stl_container_of_scalar_tag) const {
						throw std::runtime_error("attributes needs to be a scalar type or a string" + p);
					}
					template<typename T> void get_attr(std::string const &, std::string const & p, T &, stl_container_of_unknown_tag) const {
						throw std::runtime_error("attributes needs to be a scalar type or a string" + p);
					}
					template<typename T> void get_attr(std::string const &, std::string const & p, T &, stl_container_of_string_tag) const {
						throw std::runtime_error("attributes needs to be a scalar type or a string" + p);
					}
					template<typename T> void get_attr(std::string const &, std::string const & p, T &, stl_container_of_container_tag) const {
						throw std::runtime_error("attributes needs to be a scalar type or a string" + p);
					}
					template<typename T> void get_attr(std::string const &, std::string const & p, T &, stl_pair_tag) const {
						throw std::runtime_error("attributes needs to be a scalar type or a string" + p);
					}
					template<typename T> void get_attr(std::string const &, std::string const & p, T &, stl_complex_tag) const {
						throw std::runtime_error("attributes needs to be a scalar type or a string" + p);
					}
					template<typename T> void get_attr(std::string const & p, std::string const & s, T & v, scalar_tag) const {
						hid_t parent_id;
						if (is_group(p))
							parent_id = H5Gopen2(_file, p.c_str(), H5P_DEFAULT);
						else if (is_data(p))
							parent_id = H5Dopen2(_file, p.c_str(), H5P_DEFAULT);
						else
							throw std::runtime_error("unknown path: " + p);
						type_type type_id(get_native_type(v));
						attribute_type attr_id(H5Aopen(parent_id, s.c_str(), H5P_DEFAULT));
						check_error(H5Aread(attr_id, type_id, &v));
						if (is_group(p))
							check_group(parent_id);
						else
							check_data(parent_id);
					}
					template<typename T> void get_attr(std::string const & p, std::string const & s, T & v, stl_string_tag) const {
						hid_t parent_id;
						if (is_group(p))
							parent_id = H5Gopen2(_file, p.c_str(), H5P_DEFAULT);
						else if (is_data(p))
							parent_id = H5Dopen2(_file, p.c_str(), H5P_DEFAULT);
						else
							throw std::runtime_error("unknown path: " + p);
						attribute_type attr_id(H5Aopen(parent_id, s.c_str(), H5P_DEFAULT));
						type_type type_id(H5Aget_type(attr_id));
						type_type native_id(H5Tget_native_type(type_id, H5T_DIR_ASCEND));
						if (H5Tget_class(type_id) == H5T_STRING) {
							v.resize(H5Tget_size(native_id));
							check_error(H5Aread(attr_id, type_id, &v[0]));
						}
						#define HDF5_GET_STRING(T)																		\
							else if (check_error(H5Tequal(																\
								check_type(H5Tcopy(native_id)), check_type(get_native_type<T>(0))						\
							)) > 0) {																					\
								T t;																					\
								get_data(p, &t);																		\
								check_error(H5Aread(attr_id, type_id, &t));												\
								v = boost::lexical_cast<std::string>(t);												\
							}
						HDF5_FOREACH_SCALAR(HDF5_GET_STRING)
						#undef HDF5_GET_STRING
						else throw std::runtime_error("error in types: " + p);
						if (is_group(p))
							check_group(parent_id);
						else
							check_data(parent_id);
					}
					template<typename T> void set_data(std::string const & p, T v, scalar_tag) const {
						type_type type_id(get_native_type(v));
						data_type data_id(create_path(p, type_id, H5Screate(H5S_SCALAR), 0));
						check_error(H5Dwrite(data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v));
					}
					template<typename T> void set_data(std::string const & p, T const & v, stl_container_of_scalar_tag) const {
						if (!v.size())
							set_data(p, static_cast<typename T::value_type const *>(NULL), 0);
						else
							set_data(p, &(const_cast<T &>(v)[0]), v.size());
					}
					template<typename T> void set_data(std::string const & p, T const & v, stl_string_tag) const {
						if (!v.size())
							set_data(p, static_cast<char const *>(NULL), 0);
						else {
							type_type type_id(H5Tcopy(H5T_C_S1));
							check_error(H5Tset_size(type_id, v.size()));
							data_type data_id(create_path(p, type_id, H5Screate(H5S_SCALAR), 0, false));
							check_error(H5Dwrite(data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(v[0])));
						}
					}
					template<typename T> void set_data(std::string const & p, T const & v, c_string_tag) const {
						set_data(p, std::string(v), stl_string_tag());
					}
					template<typename T> void set_data(std::string const & p, T const & v, stl_pair_tag) const {
						const_cast<archive<Tag> &>(*this) << make_pvp(p + "/first", v.first);
						const_cast<archive<Tag> &>(*this) << make_pvp(p + "/second", v.second);
					}
					template<typename T> void set_data(std::string const & p, T const & v, stl_complex_tag) const {
						set_data(p, static_cast<typename T::value_type const *>(&v), 2);
					}
					template<typename T> void set_data(std::string const & p, T const & v, stl_container_of_unknown_tag) const {
						if (!v.size())
							set_data(p, static_cast<int const *>(NULL), 0);
						else {
							std::size_t pos = 0;
							for (typename T::const_iterator it = v.begin(); it != v.end(); ++it)
								const_cast<archive<Tag> &>(*this) << make_pvp(p + "/" + boost::lexical_cast<std::string>(pos++), *it);
							const_cast<archive<Tag> *>(this)->serialize(p + "/@length", pos);
						}
					}
					template<typename T> void set_data(std::string const & p, T const & v, stl_container_of_string_tag) const {
						if (!v.size())
							set_data(p, static_cast<typename T::value_type::value_type const *>(NULL), 0);
						else {
							type_type type_id(H5Tcopy(H5T_C_S1));
							check_error(H5Tset_size(type_id, H5T_VARIABLE));
							hsize_t s = v.size();
							data_type data_id(create_path(p, type_id, H5Screate_simple(1, &s, NULL), 1, &s, false));
							for (std::size_t i = 0; i < v.size(); ++i) {
								hsize_t start = i, count = 1;
								space_type space_id(H5Dget_space(data_id));
								check_error(H5Sselect_hyperslab(space_id, H5S_SELECT_SET, &start, NULL, &count, NULL));
								space_type mem_id(H5Screate_simple(1, &count, NULL));
								char const * c = const_cast<T &>(v)[i].c_str();
								check_error(H5Dwrite(data_id, type_id, mem_id, space_id, H5P_DEFAULT, &c));
							}
						}
					}
					template<typename T> void set_data(std::string const & p, T const & v, stl_container_of_container_tag) const {
						if (!v.size() || !v[0].size())
							set_data(p, static_cast<typename T::value_type::value_type const *>(NULL), 0);
						else {
							type_type type_id(get_native_type(const_cast<T &>(v)[0][0]));
							hsize_t s[2] = { v.size(), v[0].size() };
							data_type data_id(create_path(p, type_id, H5Screate_simple(2, s, NULL), 2, s));
							for (std::size_t i = 0; i < v.size(); ++i)
								if (v[i].size() != v[0].size())
									throw std::runtime_error(p + " is not a rectengual matrix");
								else {
									hsize_t start[2] = { i, 0 }, count[2] = { 1, v[i].size() };
									space_type space_id(H5Dget_space(data_id));
									check_error(H5Sselect_hyperslab(space_id, H5S_SELECT_SET, start, NULL, count, NULL));
									space_type mem_id(H5Screate_simple(2, count, NULL));
									check_error(H5Dwrite(data_id, type_id, mem_id, space_id, H5P_DEFAULT, &(const_cast<T &>(v)[i][0])));
								}
						}
					}
					template<typename T> void set_attr(std::string const &, std::string const & p, T const &, stl_container_of_scalar_tag) const {
						throw std::runtime_error("attributes needs to be a scalar type or a string" + p);
					}
					template<typename T> void set_attr(std::string const &, std::string const & p, T const &, stl_container_of_container_tag) const {
						throw std::runtime_error("attributes needs to be a scalar type or a string" + p);
					}
					template<typename T> void set_attr(std::string const &, std::string const & p, T const &, stl_container_of_unknown_tag) const {
						throw std::runtime_error("attributes needs to be a scalar type or a string" + p);
					}
					template<typename T> void set_attr(std::string const &, std::string const & p, T const &, stl_container_of_string_tag) const {
						throw std::runtime_error("attributes needs to be a scalar type or a string" + p);
					}
					template<typename T> void set_attr(std::string const &, std::string const & p, T const &, stl_pair_tag) const {
						throw std::runtime_error("attributes needs to be a scalar type or a string" + p);
					}
					template<typename T> void set_attr(std::string const &, std::string const & p, T const &, stl_complex_tag) const {
						throw std::runtime_error("attributes needs to be a scalar type or a string" + p);
					}
					template<typename T> void set_attr(std::string const & p, std::string const & s, T const & v, scalar_tag) const {
						hid_t parent_id;
						if (is_group(p))
							parent_id = H5Gopen2(_file, p.c_str(), H5P_DEFAULT);
						else if (is_data(p))
							parent_id = H5Dopen2(_file, p.c_str(), H5P_DEFAULT);
						else
							throw std::runtime_error("unknown path: " + p);
						type_type type_id(get_native_type(v));
						hid_t id = H5Aopen(parent_id, s.c_str(), H5P_DEFAULT);
						if (id < 0) {
							space_type space_id(H5Screate(H5S_SCALAR));
							id = H5Acreate(parent_id, s.c_str(), type_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
						}
						attribute_type attr_id(id);
						check_error(H5Awrite(attr_id, type_id, &v));
						if (is_group(p))
							check_group(parent_id);
						else
							check_data(parent_id);
					}
					template<typename T> void set_attr(std::string const & p, std::string const & s, T const & v, stl_string_tag) {
						hid_t parent_id;
						if (is_group(p))
							parent_id = H5Gopen2(_file, p.c_str(), H5P_DEFAULT);
						else if (is_data(p))
							parent_id = H5Dopen2(_file, p.c_str(), H5P_DEFAULT);
						else
							throw std::runtime_error("unknown path: " + p);
						type_type type_id(H5Tcopy(H5T_C_S1));
						check_error(H5Tset_size(type_id, v.size()));
						hid_t id = H5Aopen(parent_id, s.c_str(), H5P_DEFAULT);
						if (id < 0) {
							space_type space_id(H5Screate(H5S_SCALAR));
							id = H5Acreate(parent_id, s.c_str(), type_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
						}
						attribute_type attr_id(id);
						check_error(H5Awrite(attr_id, type_id, &v[0]));
						if (is_group(p))
							check_group(parent_id);
						else
							check_data(parent_id);
					}
					template<typename T> void set_attr(std::string const & p, std::string const & s, T const & v, c_string_tag) {
						set_attr(p, s, std::string(v), stl_string_tag());
					}
					std::string _context;
					file_type _file;
			};
		}
		typedef detail::archive<detail::read> iarchive;
		typedef detail::archive<detail::write> oarchive;
		template <typename T> iarchive & serialize(iarchive & ar, std::string const & p, T & v) {
			ar.serialize(p, v);
			return ar;
		}
		template <typename T> oarchive & serialize(oarchive & ar, std::string const & p, T const & v) {
			ar.serialize(p, v);
			return ar;
		}
		namespace detail {
			typedef enum {plain, ptr, array} pvp_type;
			template <typename T, pvp_type U = plain> class pvp {
				public:
					pvp(std::string const & p, T v): _p(p), _v(v) {}
					pvp(pvp<T, U> const & c): _p(c._p), _v(c._v) {}
					template<typename Tag> archive<Tag> & serialize(archive<Tag> & ar) const { return ::alps::hdf5::serialize(ar, _p, _v); }
				private:
					std::string _p;
					mutable T _v;
			};
			template <typename T> class pvp<T, ptr> {
				public:
					pvp(std::string const & p, T * v): _p(p), _v(v) {}
					pvp(pvp<T, ptr> const & c): _p(c._p), _v(c._v) {}
					template<typename Tag> archive<Tag> & serialize(archive<Tag> & ar) const { return ::alps::hdf5::serialize(ar, _p, *_v); }
				private:
					std::string _p;
					mutable T * _v;
			};
			template <> class pvp<char const, ptr> : public pvp<char const *, plain> {
				public:
					pvp(std::string const & p, char const * v): pvp<char const *, plain>(p, v) {}
					pvp(pvp<char const, ptr> const & c): pvp<char const *, plain>(*dynamic_cast<pvp<char const *, plain> const * >(&c)) {}
			};
			template <typename T> class pvp<T, array> {
				public:
					pvp(std::string const & p, T * v, std::size_t s): _p(p), _v(v), _s(s) {}
					pvp(pvp<T, ptr> const & c): _p(c._p), _v(c._v), _s(c._s) {}
					archive<write> & serialize(archive<write> & ar) const { 
						ar.set_data(ar.compute_path(_p), _v, _s);
						return ar;
					}
					archive<read> & serialize(archive<read> & ar) const { 
						ar.get_data(ar.compute_path(_p), _v); 
						return ar;
					}
				private:
					std::string _p;
					mutable T * _v;
					std::size_t _s;
			};
			template <typename T, pvp_type U> archive<write> & operator<< (archive<write> & ar, pvp<T, U> const & v) { 
				return v.serialize(ar);
			}
			template <typename T, pvp_type U> archive<read> & operator>> (archive<read> & ar, pvp<T, U> const & v) { 
				return v.serialize(ar);
			}
			#undef HDF5_ADD_CV
			#undef HDF5_FOREACH_SCALAR
		}
	}
	template <typename T> hdf5::detail::pvp<T &, hdf5::detail::plain> make_pvp(std::string const & p, T & v) {
		return hdf5::detail::pvp<T &, hdf5::detail::plain>(p, v);
	}
	template <typename T> hdf5::detail::pvp<T const &, hdf5::detail::plain> make_pvp(std::string const & p, T const & v) {
		return hdf5::detail::pvp<T const &, hdf5::detail::plain>(p, v);
	}
	template <typename T> hdf5::detail::pvp<T, hdf5::detail::ptr> make_pvp(std::string const & p, T * & v) {
		return hdf5::detail::pvp<T, hdf5::detail::ptr>(p, v);
	}
	template <typename T> hdf5::detail::pvp<T, hdf5::detail::ptr> make_pvp(std::string const & p, T * const & v) {
		return hdf5::detail::pvp<T, hdf5::detail::ptr>(p, v);
	}
	template <typename T> hdf5::detail::pvp<T, hdf5::detail::array> make_pvp(std::string const & p, T * v, std::size_t s) {
		return hdf5::detail::pvp<T, hdf5::detail::array>(p, v, s);
	}
}
#endif
