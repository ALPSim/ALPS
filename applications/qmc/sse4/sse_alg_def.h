/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2010 by Sergei Isakov <isakov@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef __SSE_ALG_DEF_H__
#define __SSE_ALG_DEF_H__

#include <alps/osiris/dump.h>

#include "lattice.h"

struct Operator {
    unsigned vertex_index;
    unsigned unit_ref;
    unsigned linked[2 * UNIT_SIZE];
};

inline alps::ODump& operator<<(alps::ODump& dump, Operator const& op)
{
    return dump << op.vertex_index << op.unit_ref;
}
 
inline alps::IDump& operator>>(alps::IDump& dump, Operator& op)
{
    return dump >> op.vertex_index >> op.unit_ref;
}

typedef std::vector<Operator>::iterator op_iterator;
typedef std::vector<Operator>::const_iterator op_c_iterator;

const unsigned IDENTITY = std::numeric_limits<unsigned>::max();
const unsigned MAX_NUMBER = std::numeric_limits<unsigned>::max();

#endif
