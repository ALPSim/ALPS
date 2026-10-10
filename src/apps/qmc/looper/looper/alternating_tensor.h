/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef LOOPER_ALTERNATING_TENSOR_H
#define LOOPER_ALTERNATING_TENSOR_H

#include <boost/tuple/tuple.hpp>

namespace looper {

//
// function alternating_tensor
//

inline int alternating_tensor(int i, int j, int k) {
  switch (i) {
  case 0 :
    switch (j) {
    case 1:
      return (k == 2) ? 1 : 0;
    case 2:
      return (k == 1) ? -1 : 0;
    default:
      return 0;
    }
  case 1 :
    switch (j) {
    case 2:
      return (k == 0) ? 1 : 0;
    case 0:
      return (k == 2) ? -1 : 0;
    default:
      return 0;
    }
  case 2 :
    switch (j) {
    case 0:
      return (k == 1) ? 1 : 0;
    case 1:
      return (k == 0) ? -1 : 0;
    default:
      return 0;
    }
  default :
    break;
  }
  return 0;
}

inline int alternating_tensor(boost::tuple<int, int, int> const& x) {
  return alternating_tensor(x.get<0>(), x.get<1>(), x.get<2>());
}

} // end namespace looper

#endif // LOOPER_ALTERNATING_TENSOR_H
