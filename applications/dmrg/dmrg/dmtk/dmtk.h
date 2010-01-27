/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
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
