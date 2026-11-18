// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
//
// How pyalps holds the observables and results of <alps/alea/observables.hpp>,
// shared by the modules that bind them (pyngsobservable_c, pyngsobservables_c,
// pyngsresult_c, pyngsresults_c, pyngsapi_c), so that each passes the same
// C++ types.
//
// - pyalps.ngs.observable is alps::alea::accumulator, a handle: the object a
//   set returns for a name refers to the observable in the set, and outlives
//   its removal from the set or the set itself.
// - pyalps.ngs.result is pyalps::result, which holds a result through a
//   pointer that is null for a result without measurements. In-place
//   arithmetic replaces what the pointer refers to.
#ifndef PYALPS_ACCUMULATOR_HANDLES_HPP
#define PYALPS_ACCUMULATOR_HANDLES_HPP

#include <alps/alea/observables.hpp>
#include <alps/alea/observables_io.hpp>
#include <alps/hdf5/archive.hpp>
#include <nanobind/nanobind.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace pyalps {

namespace alea = alps::alea;

/// A result, or none ("No Measurements").
struct result {
    std::shared_ptr<alea::observable_result> value;

    alea::observable_result & get() const {
        if (!value)
            throw nanobind::value_error("result has no measurements");
        return *value;
    }
    bool empty() const { return !value || value->count() == 0; }
};

inline result wrap(alea::observable_result const & value) {
    return result{std::make_shared<alea::observable_result>(value)};
}

/// An observable with binning analysis and bins for jackknife, as ALEA's
/// RealObservable and RealVectorObservable had. A bin number of 0 asks for
/// the default, 128; bins come in pairs, so an odd number is rounded up.
inline alea::observable declare(std::string const & name, std::uint32_t binnum) {
    return alea::observable(name).bins(binnum ? binnum + binnum % 2 : 128);
}

/// Takes numbers only.
inline alea::accumulator real_observable(std::string const & name, std::uint32_t binnum = 0) {
    return alea::accumulator(declare(name, binnum).components(1));
}

/// Takes arrays, all of the size of the first.
inline alea::accumulator real_vector_observable(std::string const & name, std::uint32_t binnum = 0) {
    return alea::accumulator(declare(name, binnum));
}

/// A set that holds `observable` alone, to merge, save or load it.
inline alea::observable_set alone(alea::accumulator const & observable) {
    alea::observable_set set;
    set.insert(observable.name(), observable);
    return set;
}

/// The set or the results stored at `path`. Loading into a new object first
/// leaves the caller's and the archive's context as they were on failure.
template <typename Set> Set load_at(alps::hdf5::archive & ar, std::string const & path) {
    alps::hdf5::archive reader(ar);
    reader.set_context(ar.complete_path(path));
    Set loaded;
    loaded.load(reader);
    return loaded;
}

template <typename Set> nanobind::object names_iterator(Set const & set) {
    nanobind::list keys;
    for (std::string const & name : set.names())
        keys.append(nanobind::cast(name));
    return keys.attr("__iter__")();
}

template <typename Set> void need(Set const & set, std::string const & name) {
    if (!set.has(name))
        throw nanobind::key_error(name.c_str());
}

} // namespace pyalps

#endif
