// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#ifndef PYALPS_MAPPING_LIFETIME_HPP
#define PYALPS_MAPPING_LIFETIME_HPP

#include <nanobind/nanobind.h>
#include <memory>
#include <string>

namespace pyalps {

// A reference_internal lookup keeps its map alive, but an erased map node
// still dies immediately. Preserve that node for any Python references to
// its value, just as Boost.Python's indexing proxies did. Keeping the actual
// node (rather than returning copies on lookup) also preserves operations
// such as observable.merge() that can replace the value's internal payload.
template <typename Mapping>
void erase_map_item(Mapping & mapping, std::string const & key) {
    namespace nb = nanobind;
    if (!mapping.has(key))
        throw nb::key_error(key.c_str());

    nb::object value = nb::cast(&mapping[key], nb::rv_policy::reference);
    using node_type = typename Mapping::node_type;
    auto pending = std::make_unique<node_type>();
    nb::capsule owner(pending.get(), [](void * pointer) noexcept {
        delete static_cast<node_type *>(pointer);
    });
    node_type * node = pending.release();
    // Establish ownership before mutating the map, so even a Python
    // allocation failure leaves every existing reference valid. Use the
    // public call policy rather than nanobind's private keep-alive API.
    nb::cpp_function([](nb::handle, nb::handle) {}, nb::keep_alive<1, 2>())(
        value, owner);
    *node = mapping.extract(key);
}

} // namespace pyalps
#endif
