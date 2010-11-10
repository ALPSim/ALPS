/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2010 by Sergei Isakov <isakov@itp.phys.ethz.ch>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
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
