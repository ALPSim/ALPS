// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
//
// How pyalps holds ALPSCore's accumulators and results, shared by the modules
// that bind them (pyngsobservable_c, pyngsobservables_c, pyngsresult_c,
// pyngsresults_c, pyngsapi_c), so that each passes the same C++ types.
//
// - ALPSCore's accumulator_set and result_set keep their items through
//   std::shared_ptr. Lookups return that pointer, so an item a Python name
//   refers to outlives its removal from the set or the set itself.
// - ALPSCore's sets have no erase, and copying one gives an empty set. Sets
//   are therefore built in place and never returned by value.
// - ALPSCore's result_wrapper has no checked empty state: its operations
//   dereference a null pointer. pyalps::result holds a result through a
//   pointer that is null for a result without measurements, and checks it.
#ifndef PYALPS_ACCUMULATOR_HANDLES_HPP
#define PYALPS_ACCUMULATOR_HANDLES_HPP

#include <alps/accumulators.hpp>
#include <alps/hdf5/archive.hpp>
#include <nanobind/nanobind.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace pyalps {

namespace acc = alps::accumulators;

/// A result, or none ("No Measurements").
struct result {
    std::shared_ptr<acc::result_wrapper> value;

    acc::result_wrapper & get() const {
        if (!value)
            throw nanobind::value_error("result has no measurements");
        return *value;
    }
    bool empty() const { return !value || value->count() == 0; }
};

/// Whether a result or accumulator holds values of type T.
template <typename T, typename W> bool holds(W & wrapper) {
    try {
        wrapper.template get<T>();
        return true;
    } catch (std::runtime_error const &) {
        return false;
    }
}

/// The pointer a set keeps for `name`.
template <typename Set>
std::shared_ptr<typename Set::value_type> item(Set & set, std::string const & name) {
    for (auto & entry : set)
        if (entry.first == name)
            return entry.second;
    throw nanobind::key_error(name.c_str());
}

/// Puts `value` into `set` under `name`, replacing an item of that name.
template <typename Set>
void assign(Set & set, std::string const & name, std::shared_ptr<typename Set::value_type> value) {
    std::vector<std::pair<std::string, std::shared_ptr<typename Set::value_type>>> kept;
    for (auto & entry : set)
        if (entry.first != name)
            kept.emplace_back(entry.first, entry.second);
    set.clear();
    for (auto & entry : kept)
        set.insert(entry.first, entry.second);
    set.insert(name, std::move(value));
}

/// Removes `name` from `set`.
template <typename Set> void erase(Set & set, std::string const & name) {
    item(set, name);
    std::vector<std::pair<std::string, std::shared_ptr<typename Set::value_type>>> kept;
    for (auto & entry : set)
        if (entry.first != name)
            kept.emplace_back(entry.first, entry.second);
    set.clear();
    for (auto & entry : kept)
        set.insert(entry.first, entry.second);
}

/// Replaces the items of `set` by those of `loaded`.
template <typename Set> void replace(Set & set, Set & loaded) {
    set.clear();
    for (auto & entry : loaded)
        set.insert(entry.first, entry.second);
}

/// An accumulator with full binning analysis, as ALEA's RealObservable
/// and RealVectorObservable had. A bin number of 0 asks for the default, 128,
/// as it did there; ALPSCore's accumulators would keep no bins at all.
template <typename T>
std::shared_ptr<acc::accumulator_wrapper> full_binning(std::string const & name, std::uint32_t binnum = 0) {
    return acc::FullBinningAccumulator<T>(name, std::size_t(binnum ? binnum : 128)).wrapper;
}

template <typename R>
void load_as(alps::hdf5::archive & ar, std::shared_ptr<acc::result_wrapper> & loaded) {
    if (!loaded && R::can_load(ar))
        loaded = std::make_shared<acc::result_wrapper>(R());
}

/// Loads the result stored at the archive's context. ALPSCore's result_set
/// finds a stored result's type in a registry that is private to the
/// library; this tries the same predefined types, in the same order (most
/// features first), for the value types pyalps handles.
inline std::shared_ptr<acc::result_wrapper> load_result(alps::hdf5::archive & ar) {
    std::shared_ptr<acc::result_wrapper> loaded;
    load_as<acc::FullBinningAccumulator<double>::result_type>(ar, loaded);
    load_as<acc::FullBinningAccumulator<std::vector<double>>::result_type>(ar, loaded);
    load_as<acc::LogBinningAccumulator<double>::result_type>(ar, loaded);
    load_as<acc::LogBinningAccumulator<std::vector<double>>::result_type>(ar, loaded);
    load_as<acc::NoBinningAccumulator<double>::result_type>(ar, loaded);
    load_as<acc::NoBinningAccumulator<std::vector<double>>::result_type>(ar, loaded);
    load_as<acc::MeanAccumulator<double>::result_type>(ar, loaded);
    load_as<acc::MeanAccumulator<std::vector<double>>::result_type>(ar, loaded);
    if (!loaded)
        throw std::runtime_error("no result of double or vector of double values is stored at "
                                 + ar.get_context());
    loaded->load(ar);
    return loaded;
}

} // namespace pyalps

#endif
