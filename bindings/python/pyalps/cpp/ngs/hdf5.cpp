// Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>
//               Matthias Troyer <troyer@comp-phys.org>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
// Save path: extract_from_pyobject_py11 dispatches a nb::handle to a
// visitor that writes a concrete C++ value. Load path: dispatches on
// the archive's inspected type and reads into a concrete C++ type
// before wrapping it back as a nb::object.
//
// Exception translation: pyalps/hdf5.py creates ArchiveError etc. and
// calls register_archive_exception_type(id, type); the translators
// below fire PyErr_SetString against whichever Python type was handed
// in.
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/pair.hpp>
#include <alps/hdf5/pointer.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/ngs/stacktrace.hpp>
#include "extract_from_pyobject.hpp"
#include "../numpy_compat.hpp"
#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
namespace nb = nanobind;
namespace alps {
    namespace detail {
        // Save-side visitor: receives a concrete C++ value (or a
        // nb::list / nb::dict) from extract_from_pyobject_py11 and
        // writes it to the archive at `path`.
        struct hdf5_save_py11_visitor {
            alps::hdf5::archive & ar;
            std::string path;
            template <typename U>
            void operator()(U const & v) const {
                ar[path] << v;
            }
            template <typename U>
            void operator()(U const * ptr, std::vector<std::size_t> const & sizes) const {
                // Use make_pvp(path, ptr, size-vector) to preserve the
                // dimensional shape — a plain vector<U> flatten would
                // round-trip the data but lose the rank.
                ar << alps::make_pvp(path, ptr, sizes);
            }
            void operator()(nb::list const & l) const {
                // Order: flat numeric first, then nested numeric, then
                // strings. Heterogeneous / deeply-nested / mixed-type
                // lists fall through to the descent branch below which
                // stores each entry under a numeric child path.
                try { ar[path] << nb::cast<std::vector<double>>(l); return; }
                catch (nb::cast_error const &) {}
                try { ar[path] << nb::cast<std::vector<int>>(l); return; }
                catch (nb::cast_error const &) {}
                try { ar[path] << nb::cast<std::vector<std::complex<double>>>(l); return; }
                catch (nb::cast_error const &) {}
                try { ar[path] << nb::cast<std::vector<std::vector<double>>>(l); return; }
                catch (nb::cast_error const &) {}
                try { ar[path] << nb::cast<std::vector<std::vector<int>>>(l); return; }
                catch (nb::cast_error const &) {}
                try { ar[path] << nb::cast<std::vector<std::string>>(l); return; }
                catch (nb::cast_error const &) {}
                // Inhomogeneous — recurse per-element into
                // <path>/<index>, letting each entry be stored as its
                // own native type.
                ar.create_group(path);
                Py_ssize_t i = 0;
                for (auto item : l) {
                    std::string child = path + "/" + std::to_string(static_cast<long long>(i++));
                    hdf5_save_py11_visitor child_visitor{ar, child};
                    extract_from_pyobject_py11(child_visitor, item);
                }
            }
            void operator()(nb::dict const & d) const {
                // Store a dict as a group with one child per key. Keys
                // are stringified (HDF5 paths are strings), values go
                // through the same save dispatch recursively.
                ar.create_group(path);
                for (auto item : d) {
                    std::string key = nb::cast<std::string>(nb::str(item.first));
                    std::string child = path + "/" + key;
                    hdf5_save_py11_visitor child_visitor{ar, child};
                    extract_from_pyobject_py11(child_visitor, item.second);
                }
            }
        };
        std::string python_hdf5_get_filename(alps::hdf5::archive & ar) {
            return ar.get_filename();
        }
        void python_hdf5_save(alps::hdf5::archive & ar,
                              std::string const & path,
                              nb::handle data) {
            hdf5_save_py11_visitor visitor{ar, path};
            extract_from_pyobject_py11(visitor, data);
        }
        // Helper: load a multi-dim HDF5 dataset of element type T
        // into a flat std::vector, then wrap as a numpy array with
        // the original shape so that Python sees a 2-D np.array for
        // rank-2 writes etc. Preserves the dimensionality encoded on
        // the save path (alps::make_pvp(path, ptr, size-vector)).
        template <typename T>
        nb::object load_nd_array(alps::hdf5::archive & ar,
                                 std::string const & path,
                                 std::vector<std::size_t> const & shape) {
            std::size_t total = 1;
            for (auto s : shape) total *= s;
            std::vector<T> flat(total);
            if (shape.size() <= 1) {
                // vector<T> overload works directly.
                ar[path] >> flat;
            } else {
                // make_pvp with explicit size-vector to read a
                // multi-dim dataset into a flat buffer.
                ar >> alps::make_pvp(path, flat.data(), shape);
            }
            return alps::python::make_numpy_array<T>(flat.data(), shape);
        }
        nb::object python_hdf5_load_impl(alps::hdf5::archive & ar,
                                         std::string const & path);
        nb::object python_hdf5_load(alps::hdf5::archive & ar,
                                    std::string const & path) {
            return python_hdf5_load_impl(ar, path);
        }
        nb::object python_hdf5_load_impl(alps::hdf5::archive & ar,
                                         std::string const & path) {
            // Groups (not datasets) get loaded recursively. Children
            // whose names are consecutive decimal integers starting at 0
            // are recovered as a Python list (preserving round-trip for
            // list-saved-as-group); otherwise a dict.
            if (ar.is_group(path)) {
                auto children = ar.list_children(path);
                bool list_shaped = true;
                for (std::size_t i = 0; list_shaped && i < children.size(); ++i) {
                    if (children[i] != std::to_string(i))
                        list_shaped = false;
                }
                if (list_shaped) {
                    nb::list result;
                    for (auto const & child : children)
                        result.append(
                            python_hdf5_load_impl(ar, path + "/" + child));
                    return nb::object(std::move(result));
                } else {
                    nb::dict result;
                    for (auto const & child : children)
                        result[nb::str(child.c_str())] =
                            python_hdf5_load_impl(ar, path + "/" + child);
                    return nb::object(std::move(result));
                }
            }
            // Complex values have a quirky HDF5 representation: a
            // single complex is stored as rank-1 dims=[2] (real,imag)
            // and a 2x2 array of complex as rank-3 dims=[2,2,2]. So
            // is_scalar returns false for a scalar complex — branch
            // on is_complex first and use the rank minus 1 (stripping
            // the trailing complex-pair dim) to tell scalar from
            // array.
            if (ar.is_complex(path)) {
                auto ext = ar.extent(path);
                if (ext.size() == 1) {
                    std::complex<double> v; ar[path] >> v; return nb::cast(v);
                }
                std::vector<std::size_t> shape(ext.begin(), ext.end() - 1);
                return load_nd_array<std::complex<double>>(ar, path, shape);
            }
            // Convenience macros for the scalar path: check each
            // candidate integer width in turn (numpy's default int is
            // platform-dependent — int64 on macOS/Linux, int32 on
            // Windows — so we can't rely on just `int` matching).
            #define TRY_SCALAR(T)                                                                \
                if (ar.is_datatype<T>(path)) { T v; ar[path] >> v; return nb::cast(v); }
            if (ar.is_scalar(path)) {
                TRY_SCALAR(std::string)
                TRY_SCALAR(double)
                TRY_SCALAR(float)
                TRY_SCALAR(bool)
                TRY_SCALAR(std::int64_t)
                TRY_SCALAR(std::int32_t)
                TRY_SCALAR(std::int16_t)
                TRY_SCALAR(std::int8_t)
                TRY_SCALAR(std::uint64_t)
                TRY_SCALAR(std::uint32_t)
                TRY_SCALAR(std::uint16_t)
                TRY_SCALAR(std::uint8_t)
                throw std::runtime_error(
                    "Unknown HDF5 scalar type at " + path + ALPS_STACKTRACE);
            } else {
                // String datasets don't map to nb::ndarray<T> the way
                // numeric types do; keep the flat vector path.
                if (ar.is_datatype<std::string>(path)) {
                    std::vector<std::string> v; ar[path] >> v; return nb::cast(v);
                }
                auto shape = ar.extent(path);
                if (ar.is_datatype<double>(path))        return load_nd_array<double>(ar, path, shape);
                if (ar.is_datatype<float>(path))         return load_nd_array<float>(ar, path, shape);
                if (ar.is_datatype<std::int64_t>(path))  return load_nd_array<std::int64_t>(ar, path, shape);
                if (ar.is_datatype<std::int32_t>(path))  return load_nd_array<std::int32_t>(ar, path, shape);
                if (ar.is_datatype<std::int16_t>(path))  return load_nd_array<std::int16_t>(ar, path, shape);
                if (ar.is_datatype<std::int8_t>(path))   return load_nd_array<std::int8_t>(ar, path, shape);
                if (ar.is_datatype<std::uint64_t>(path)) return load_nd_array<std::uint64_t>(ar, path, shape);
                if (ar.is_datatype<std::uint32_t>(path)) return load_nd_array<std::uint32_t>(ar, path, shape);
                if (ar.is_datatype<std::uint16_t>(path)) return load_nd_array<std::uint16_t>(ar, path, shape);
                if (ar.is_datatype<std::uint8_t>(path))  return load_nd_array<std::uint8_t>(ar, path, shape);
                throw std::runtime_error(
                    "Unknown HDF5 vector type at " + path + ALPS_STACKTRACE);
            }
            #undef TRY_SCALAR
        }
        nb::list python_hdf5_extent(alps::hdf5::archive & ar,
                                    std::string const & path) {
            nb::list result;
            std::vector<std::size_t> ext = ar.extent(path);
            if (ar.is_complex(path)) {
                if (ext.size() > 1)
                    ext.pop_back();
                else
                    ext.back() = 1;
            }
            for (auto const & s : ext)
                result.append(s);
            return result;
        }
        // Python exception types registered by pyalps.hdf5 at import
        // time. Translators below fire PyErr_SetString against these
        // pre-registered PyObject*'s so the Python side sees its own
        // subclasses (ArchiveError / ArchiveNotFound / ...).
        std::array<PyObject *, 6> exception_type = {};
        #define TRANSLATE_CPP_ERROR_TO_PYTHON(T, ID)                            \
            static void translate_ ## T (hdf5:: T const & e) {                  \
                std::string message =                                           \
                    std::string(e.what()).substr(                               \
                        0, std::string(e.what()).find_first_of('\n'));          \
                PyErr_SetString(exception_type[ID] ? exception_type[ID]         \
                                                   : PyExc_RuntimeError,       \
                                message.c_str());                               \
            }
        TRANSLATE_CPP_ERROR_TO_PYTHON(archive_error, 0)
        TRANSLATE_CPP_ERROR_TO_PYTHON(archive_not_found, 1)
        TRANSLATE_CPP_ERROR_TO_PYTHON(archive_closed, 2)
        TRANSLATE_CPP_ERROR_TO_PYTHON(invalid_path, 3)
        TRANSLATE_CPP_ERROR_TO_PYTHON(path_not_found, 4)
        TRANSLATE_CPP_ERROR_TO_PYTHON(wrong_type, 5)
        #undef TRANSLATE_CPP_ERROR_TO_PYTHON
        void register_exception_type(int id, nb::object type) {
            if (id < 0 || id >= static_cast<int>(exception_type.size()))
                throw std::out_of_range(
                    "register_archive_exception_type: id out of range");
            // Py_INCREF the incoming type so it survives past this call
            // (we're keeping a raw PyObject* in a static array).
            Py_INCREF(type.ptr());
            exception_type[id] = type.ptr();
        }
    }
}
NB_MODULE(pyngshdf5_c, m) {
    // Install the six C++→Python exception translators. Each calls the
    // matching translate_* above, which forwards to whichever Python
    // class was registered via register_archive_exception_type. If
    // pyalps/hdf5.py hasn't run yet, the translator falls back to
    // RuntimeError so the module is safely loadable on its own.
    nb::register_exception_translator(
        [](const std::exception_ptr &p, void * /*payload*/) {
            try { std::rethrow_exception(p); }
            catch (alps::hdf5::archive_not_found const & e) {
                alps::detail::translate_archive_not_found(e);
            } catch (alps::hdf5::archive_closed const & e) {
                alps::detail::translate_archive_closed(e);
            } catch (alps::hdf5::invalid_path const & e) {
                alps::detail::translate_invalid_path(e);
            } catch (alps::hdf5::path_not_found const & e) {
                alps::detail::translate_path_not_found(e);
            } catch (alps::hdf5::wrong_type const & e) {
                alps::detail::translate_wrong_type(e);
            } catch (alps::hdf5::archive_error const & e) {
                // Base class — must be caught LAST since the specialized
                // types above inherit from it.
                alps::detail::translate_archive_error(e);
            }
        });
    m.def("register_archive_exception_type",
          &alps::detail::register_exception_type);
    nb::class_<alps::hdf5::archive>(m, "hdf5_archive_impl")
        .def(nb::init<std::string, std::string>())
        .def("__deepcopy__",
             // copy.deepcopy() hands us (self, memo); memo unused.
             [](alps::hdf5::archive const & self, nb::handle /*memo*/) {
                 return alps::hdf5::archive(self);
             })
        .def_prop_ro("filename", &alps::detail::python_hdf5_get_filename)
        .def_prop_ro("context",  &alps::hdf5::archive::get_context)
        .def_prop_ro("is_open",  &alps::hdf5::archive::is_open)
        .def("set_context",     &alps::hdf5::archive::set_context)
        .def("is_group",        &alps::hdf5::archive::is_group)
        .def("is_data",         &alps::hdf5::archive::is_data)
        .def("is_attribute",    &alps::hdf5::archive::is_attribute)
        .def("close",           &alps::hdf5::archive::close)
        .def("extent",          &alps::detail::python_hdf5_extent)
        .def("dimensions",      &alps::hdf5::archive::dimensions)
        .def("is_scalar",       &alps::hdf5::archive::is_scalar)
        .def("is_complex",      &alps::hdf5::archive::is_complex)
        .def("is_null",         &alps::hdf5::archive::is_null)
        .def("list_children",   &alps::hdf5::archive::list_children)
        .def("list_attributes", &alps::hdf5::archive::list_attributes)
        .def("__setitem__",     &alps::detail::python_hdf5_save)
        .def("__getitem__",     &alps::detail::python_hdf5_load)
        .def("create_group",    &alps::hdf5::archive::create_group)
        .def("delete_data",     &alps::hdf5::archive::delete_data)
        .def("delete_group",    &alps::hdf5::archive::delete_group)
        .def("delete_attribute",&alps::hdf5::archive::delete_attribute);
}
