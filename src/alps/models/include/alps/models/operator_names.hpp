// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The XML catalog's spellings (lib/xml/models.xml) of the kit's built-in
// site operators. SiteBasis::op accepts the catalog spelling of a kit
// operator, and CatalogModel::compile_operator the kit spelling of a catalog
// one, so one name reaches the same operator on either side. Header-only, so
// the catalog importer needs no ALPS::models link for it.
#pragma once

#include <array>
#include <string_view>
#include <utility>

namespace alps {
namespace models {

// (kit name, catalog name). A catalog name listed twice resolves to the first
// kit name its site defines; Sz, n, n_up, c, c_up and b agree already.
inline constexpr std::array<std::pair<std::string_view, std::string_view>, 11> catalog_operator_names{{
    {"Sp", "Splus"},
    {"Sm", "Sminus"},
    {"c_dag", "cdag"},
    {"c_dag_up", "cdag_up"},
    {"c_dn", "c_down"},
    {"c_dag_dn", "cdag_down"},
    {"n_dn", "n_down"},
    {"n_updn", "double_occupancy"},
    {"b_dag", "bdag"},
    {"a", "b"},
    {"a_dag", "bdag"},
}};

} // namespace models
} // namespace alps
