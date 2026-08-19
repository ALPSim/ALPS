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
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
namespace nb = nanobind;
namespace alps {
    namespace detail {
        // Decode one layer of the two entities produced by
        // archive::encode_segment.  Do this locally instead of calling
        // archive::decode_segment unconditionally: HDF5 files created by
        // other tools may legitimately contain a raw '&' in a child name.
        // A literal name containing exactly "&#38;" or "&#47;" remains
        // ambiguous because the existing ALPS format has no type marker.
        static std::string decode_dict_key(std::string const & segment) {
            std::string result;
            result.reserve(segment.size());
            for (std::size_t pos = 0; pos < segment.size();) {
                if (segment.compare(pos, 5, "&#38;") == 0) {
                    result.push_back('&');
                    pos += 5;
                } else if (segment.compare(pos, 5, "&#47;") == 0) {
                    result.push_back('/');
                    pos += 5;
                } else {
                    result.push_back(segment[pos++]);
                }
            }
            return result;
        }
        // Analysis of a Python list/tuple tree against the legacy
        // Boost.Python vectorization rules (src/alps/hdf5/python.cpp,
        // is_vectorizable_generic): a list is written as one dataset
        // only when every leaf has the SAME exact scalar type — plain
        // bool was never a vectorizable dtype — and all nested extents
        // are rectangular. Everything else becomes a group with one
        // child per index. The checked-in pyhdf5io fixture documents
        // this contract ([1, 2, 3] must stay int32 on disk).
        struct list_vectorizer {
            enum class leaf_kind { none, integral, floating, cplx, text };
            std::vector<std::size_t> extent;   // rectangular extents per depth
            leaf_kind kind = leaf_kind::none;  // != none also means "a leaf was seen"
            std::vector<long long> ints;
            std::vector<double> reals;
            std::vector<std::complex<double>> cplxs;
            std::vector<std::string> texts;
            bool fits_int = true;
            bool analyze(nb::handle node, std::size_t depth) {
                std::size_t const n = nb::len(node);
                if (depth == extent.size())
                    extent.push_back(n);
                else if (extent[depth] != n)
                    return false;                                   // ragged
                for (std::size_t i = 0; i < n; ++i) {
                    nb::object item = node[i];
                    PyObject * raw = item.ptr();
                    if (PyBool_Check(raw))
                        return false;         // legacy: bool never vectorizes
                    if (PyList_Check(raw) || PyTuple_Check(raw)) {
                        // once a leaf has fixed the depth (extent can no
                        // longer grow), sequences may not appear at or
                        // below the leaf level
                        if (kind != leaf_kind::none && depth + 1 >= extent.size())
                            return false;
                        if (!analyze(item, depth + 1))
                            return false;
                        continue;
                    }
                    // scalar leaf: all leaves sit at one depth — the
                    // deepest extent recorded so far
                    if (depth + 1 != extent.size())
                        return false;
                    // Exact numeric types only, mirroring the legacy
                    // tp_name dispatch: numpy scalars (np.float64 and
                    // np.complex128 included, although they subclass the
                    // builtins) take the numpy-stacking path below,
                    // which preserves their dtype like the legacy
                    // scalar_types table did. str accepts subclasses —
                    // np.str_ was a legacy string dtype too.
                    if (PyLong_Check(raw)) {   // np ints don't subclass int
                        int overflow = 0;
                        long long v = PyLong_AsLongLongAndOverflow(raw, &overflow);
                        if (overflow)
                            return false;   // → descent; the per-element save raises, like legacy
                        if (!accept(leaf_kind::integral))
                            return false;
                        if (v < std::numeric_limits<int>::min() || v > std::numeric_limits<int>::max())
                            fits_int = false;
                        ints.push_back(v);
                    } else if (PyFloat_CheckExact(raw)) {
                        if (!accept(leaf_kind::floating))
                            return false;
                        reals.push_back(PyFloat_AsDouble(raw));
                    } else if (PyComplex_CheckExact(raw)) {
                        if (!accept(leaf_kind::cplx))
                            return false;
                        Py_complex c = PyComplex_AsCComplex(raw);
                        cplxs.emplace_back(c.real, c.imag);
                    } else if (PyUnicode_Check(raw)) {
                        if (!accept(leaf_kind::text))
                            return false;
                        texts.push_back(nb::cast<std::string>(item));
                    } else
                        return false;   // numpy scalars/arrays, other objects
                }
                return true;
            }
            bool accept(leaf_kind k) {
                if (kind == leaf_kind::none)
                    kind = k;
                return kind == k;       // legacy: mixed scalar kinds → group
            }
        };
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
                // NumPy uses rank zero for a 0-D array.  Passing an empty
                // size vector to archive::write creates an HDF5 NULL
                // dataspace, which silently turns the scalar into an empty
                // array.  Store the pointed-to value as a scalar instead.
                if (sizes.empty()) {
                    ar[path] << *ptr;
                    return;
                }
                // Use make_pvp(path, ptr, size-vector) to preserve the
                // dimensional shape — a plain vector<U> flatten would
                // round-trip the data but lose the rank.
                ar << alps::make_pvp(path, ptr, sizes);
            }
            void operator()(nb::list const & l) const {
                // Reproduce the legacy vectorization rules (see
                // list_vectorizer above): exact-type homogeneous
                // rectangular list trees become one N-D dataset that
                // keeps the element type — [1, 2, 3] stays int32 on
                // disk, floats stay float64 — while bool-containing,
                // mixed-type and ragged lists become a group with one
                // child per index.
                if (nb::len(l) == 0) {
                    // legacy wrote an empty integer dataset
                    ar[path] << std::vector<int>();
                    return;
                }
                list_vectorizer v;
                if (v.analyze(l, 0)) {
                    switch (v.kind) {
                        case list_vectorizer::leaf_kind::integral:
                            if (v.fits_int) {
                                std::vector<int> buf(v.ints.begin(), v.ints.end());
                                (*this)(buf.data(), v.extent);
                            } else
                                (*this)(v.ints.data(), v.extent);
                            return;
                        case list_vectorizer::leaf_kind::floating:
                            (*this)(v.reals.data(), v.extent);
                            return;
                        case list_vectorizer::leaf_kind::cplx:
                            (*this)(v.cplxs.data(), v.extent);
                            return;
                        case list_vectorizer::leaf_kind::text:
                            if (v.extent.size() == 1) {
                                ar[path] << v.texts;
                                return;
                            }
                            break;   // nested string lists → group descent
                        case list_vectorizer::leaf_kind::none:
                            break;   // e.g. [[], []] → group descent
                    }
                } else if (numpy_stackable(l)) {
                    // Legacy vectorized numpy content too: homogeneous
                    // numpy-scalar lists (numpy.int64 etc. were
                    // scalar_types entries) and rectangular trees
                    // mixing ndarrays with nested sequences all became
                    // one dataset. Delegate to numpy so shape checking
                    // and dtype handling match numpy's rules, then feed
                    // the stacked array through the ndarray save path.
                    // Ragged shapes (numpy raises) and non-numeric
                    // dtypes fall through to the group descent below.
                    nb::object arr;
                    try {
                        arr = nb::borrow<nb::object>(alps::python::numpy_module())
                                  .attr("asarray")(l);
                    } catch (nb::python_error &) {
                        arr = nb::object();
                    }
                    if (arr.is_valid()) {
                        std::string const dtype_kind =
                            nb::cast<std::string>(arr.attr("dtype").attr("kind"));
                        if (dtype_kind.size() == 1
                            && std::strchr("biufc", dtype_kind[0])) {
                            hdf5_save_py11_visitor child_visitor{ar, path};
                            extract_from_pyobject_py11(child_visitor, arr);
                            return;
                        }
                    }
                }
                // Heterogeneous / ragged / bool-containing — recurse
                // per-element into <path>/<index>, letting each entry
                // be stored as its own native type (legacy behaviour).
                // Legacy wiped any existing group before a list save;
                // create_group alone would keep stale children around.
                if (ar.is_group(path))
                    ar.delete_group(path);
                ar.create_group(path);
                Py_ssize_t i = 0;
                for (auto item : l) {
                    std::string child = path + "/" + std::to_string(static_cast<long long>(i++));
                    hdf5_save_py11_visitor child_visitor{ar, child};
                    extract_from_pyobject_py11(child_visitor, item);
                }
            }
            static bool is_ndarray(PyObject * raw) {
                return std::strcmp(Py_TYPE(raw)->tp_name, "numpy.ndarray") == 0;
            }
            static bool is_numpy_scalar(PyObject * raw) {
                static std::array<char const *, 16> const scalar_types{{
                    "numpy.str_", "numpy.str", "numpy.bool_", "numpy.bool",
                    "numpy.int8", "numpy.int16", "numpy.int32", "numpy.int64",
                    "numpy.uint8", "numpy.uint16", "numpy.uint32", "numpy.uint64",
                    "numpy.float32", "numpy.float64",
                    "numpy.complex64", "numpy.complex128",
                }};
                for (char const * scalar_type : scalar_types)
                    if (std::strcmp(Py_TYPE(raw)->tp_name, scalar_type) == 0)
                        return true;
                return false;
            }
            struct tree_scan {
                bool has_ndarray = false;
                bool has_numpy_scalar = false;
                bool has_other_scalar = false;
                bool has_bool_leaf = false;
                bool homogeneous_numpy_scalars = true;
                PyTypeObject * numpy_scalar_type = nullptr;
            };
            static void scan_tree(nb::handle node, tree_scan & scan) {
                std::size_t const n = nb::len(node);
                for (std::size_t i = 0; i < n; ++i) {
                    nb::object item = node[i];
                    PyObject * raw = item.ptr();
                    if (is_ndarray(raw)) {
                        scan.has_ndarray = true;
                    } else if (PyList_Check(raw) || PyTuple_Check(raw)) {
                        scan_tree(item, scan);
                    } else if (is_numpy_scalar(raw)) {
                        scan.has_numpy_scalar = true;
                        PyTypeObject * scalar_type = Py_TYPE(raw);
                        if (!scan.numpy_scalar_type)
                            scan.numpy_scalar_type = scalar_type;
                        else if (scan.numpy_scalar_type != scalar_type)
                            scan.homogeneous_numpy_scalars = false;
                        if (std::strncmp(Py_TYPE(raw)->tp_name, "numpy.bool", 10) == 0)
                            scan.has_bool_leaf = true;
                    } else {
                        scan.has_other_scalar = true;
                        if (PyBool_Check(raw))
                            scan.has_bool_leaf = true;
                    }
                }
            }
            // The list shapes the legacy build stacked into one
            // dataset beyond plain scalars: (a) numpy scalars of ONE
            // type (exact tp_name match, like legacy scalar_types), or
            // (b) sequences/ndarrays only, with an ndarray somewhere in
            // the tree (legacy vectorized extent-matched mixes of
            // list/tuple/ndarray nodes) — but never when a plain bool
            // sits among the leaves, which numpy would silently promote
            // to 0/1. Pure-list trees never reach (b) — their
            // exact-type handling stays with list_vectorizer.
            static bool numpy_stackable(nb::list const & l) {
                char const * first_scalar = nullptr;
                bool scalars_only = true;
                bool sequences_only = true;
                for (auto item : l) {
                    PyObject * raw = item.ptr();
                    char const * tp = Py_TYPE(raw)->tp_name;
                    if (is_ndarray(raw) || PyList_Check(raw) || PyTuple_Check(raw)) {
                        scalars_only = false;
                        continue;
                    }
                    sequences_only = false;
                    if (std::strncmp(tp, "numpy.", 6) != 0)
                        return false;
                    if (!first_scalar)
                        first_scalar = tp;
                    else if (std::strcmp(tp, first_scalar) != 0)
                        return false;
                }
                if (scalars_only && first_scalar)
                    return true;
                if (!sequences_only)
                    return false;
                tree_scan scan;
                scan_tree(l, scan);
                // A rectangular tree made solely from one exact NumPy
                // scalar type is vectorizable at any nesting depth.
                // np.asarray below performs the final rectangularity
                // check and preserves the scalar dtype.
                if (scan.has_numpy_scalar && !scan.has_ndarray
                    && !scan.has_other_scalar)
                    return scan.homogeneous_numpy_scalars;
                // Preserve the legacy ndarray/list stacking path.  Bool
                // leaves remain a veto because NumPy would silently turn
                // them into 0/1 when combined with a numeric ndarray.
                return scan.has_ndarray && !scan.has_bool_leaf;
            }
            void operator()(nb::dict const & d) const {
                // Store a dict as a group with one child per key. Keys
                // are stringified (HDF5 paths are strings), values go
                // through the same save dispatch recursively. Like the
                // list descent above (and the legacy build), wipe an
                // existing group first so stale keys don't survive.
                if (ar.is_group(path))
                    ar.delete_group(path);
                ar.create_group(path);
                for (auto item : d) {
                    std::string key = nb::cast<std::string>(nb::str(item.first));
                    std::string child = path + "/" + ar.encode_segment(key);
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
            // archive::read rejects a zero-sized chunk.  The HDF5 dataset
            // already carries the complete extent, so for arrays such as
            // (0, 2) and (2, 0) there is no payload to read: construct the
            // correctly shaped NumPy array directly.
            if (total == 0) {
                return alps::python::make_numpy_array<T>(nullptr, shape);
            } else if (shape.size() <= 1) {
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
            // whose names are exactly the decimal integers 0..n-1 are
            // recovered as a Python list (preserving round-trip for
            // list-saved-as-group); otherwise a dict. The backend
            // yields child names in lexicographic order ("0", "1",
            // "10", "2", ...), so the check is on the name SET and the
            // list is loaded in numeric order — the legacy loader was
            // order-insensitive the same way, indexing
            // value[cast<size_t>(name)].
            if (ar.is_group(path)) {
                auto children = ar.list_children(path);
                // Match the legacy dynamic loader: an empty group is a
                // dict. Empty lists use the dataset representation.
                bool list_shaped = !children.empty();
                std::vector<bool> seen(children.size(), false);
                for (auto const & child : children) {
                    bool numeric = !child.empty() && child.size() < 20;
                    for (char c : child)
                        if (c < '0' || c > '9') {
                            numeric = false;
                            break;
                        }
                    std::size_t index = numeric
                        ? static_cast<std::size_t>(std::strtoull(child.c_str(), nullptr, 10))
                        : 0;
                    if (!numeric || std::to_string(index) != child
                        || index >= children.size() || seen[index]) {
                        list_shaped = false;
                        break;
                    }
                    seen[index] = true;
                }
                if (list_shaped) {
                    nb::list result;
                    for (std::size_t i = 0; i < children.size(); ++i)
                        result.append(
                            python_hdf5_load_impl(ar, path + "/" + std::to_string(i)));
                    return nb::object(std::move(result));
                } else {
                    nb::dict result;
                    for (auto const & child : children) {
                        std::string const key = decode_dict_key(child);
                        result[nb::str(key.c_str())] =
                            python_hdf5_load_impl(ar, path + "/" + child);
                    }
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
                bool const single_value = ext.size() == 1;
                // Preserve the component precision.  The legacy loader
                // returned complex64 datasets as NumPy complex64 rather than
                // widening them to complex128; only a complex128 scalar used
                // the ordinary Python ``complex`` shortcut.
                if (ar.is_datatype<float>(path)) {
                    if (single_value) {
                        std::complex<float> value;
                        ar[path] >> value;
                        return alps::python::make_numpy_array(
                            &value, std::vector<std::size_t>());
                    }
                    std::vector<std::size_t> shape(ext.begin(), ext.end() - 1);
                    return load_nd_array<std::complex<float>>(ar, path, shape);
                }
                if (single_value) {
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
            // Keep a strong reference in the static table — the entry
            // is deliberately pinned until process exit because the
            // translators can fire at any time — but release any
            // previous entry so re-registration doesn't leak it.
            PyObject * previous = exception_type[id];
            Py_INCREF(type.ptr());
            exception_type[id] = type.ptr();
            Py_XDECREF(previous);
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
