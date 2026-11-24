// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT

#ifndef ALPS_DMFT_CTINT_MAX_ORDER_H
#define ALPS_DMFT_CTINT_MAX_ORDER_H

/// @file ctint_max_order.h
/// @brief the warning of both CT-INT solvers when the expansion reached MAX_ORDER
///
/// Both solvers refuse to add a vertex at MAX_ORDER vertices, so they sample
/// the series cut at that order. That can be what is wanted, to compare with
/// low-order methods, but it changes the results, and it was silent (#389).
/// The fraction of the updates spent at MAX_ORDER is measured as
/// "MaxOrderReached"; a nonzero one is reported with this warning.

#include <ostream>

/// Warns on `os` if `fraction` (of the updates at order MAX_ORDER) is nonzero.
inline void warn_if_max_order_reached(std::ostream& os, double fraction, unsigned int max_order)
{
  if (!(fraction > 0))
    return;
  os << "WARNING: CT-INT: the expansion order was at MAX_ORDER = " << max_order << " in " << 100 * fraction
     << "% of the updates. The series is cut at that order: the results are those of the expansion up to order "
     << max_order << ", not of the full model. Increase MAX_ORDER for the full model." << std::endl;
}

#endif
