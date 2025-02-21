/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2010 by Sergei Isakov <isakov@itp.phys.ethz.ch>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
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
