/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
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

#ifndef __DMTK_H__
#define __DMTK_H__

namespace dmtk
{

enum
{
  LEFT,
  RIGHT,
};

} //namespace dmtk

#include <math.h>
#include <valarray>
#include <dmtk/conj.h>
#include <dmtk/ctimer.h>
#include <dmtk/bits.h>
#include <dmtk/enums.h>
#include <dmtk/basis.h>
#include <dmtk/vector.h>
#include <dmtk/matrix.h>
#include <dmtk/qn.h>
#include <dmtk/subspace.h>
#include <dmtk/operators.h>
#include <dmtk/block_matrix.h>
#include <dmtk/block.h>
#include <dmtk/state.h>
#include <dmtk/state_slice.h>
#include <dmtk/lattice.h>
#include <dmtk/system.h>
#include <dmtk/util.h>
#include <dmtk/hami.h>

#endif // __DMTK_H__
