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

#ifndef __DMTK_ENUMS_H__
#define __DMTK_ENUMS_H__

namespace dmtk
{

#define DMTK_SQRT2 1.41421356237309504880168872420969807856967
#define DMTK_EULER 0.5772156649015328606065120900824024310422
#define DMTK_PI 3.141592653589793238462643383279502884197
#define DMTK_PIO2 1.57079632679489661923132169163975144209858
#define DMTK_TWOPI 6.283185307179586476925286766559005768394
#define DMTK_ERROR 2147483647 // largest positive 32 bit int 

enum
{
  LEFT2RIGHT,
  RIGHT2LEFT,
};

enum
{
  MASK_BLOCK1 = 1 << 0,
  MASK_BLOCK2 = 1 << 1,
  MASK_BLOCK3 = 1 << 2,
  MASK_BLOCK4 = 1 << 3,
};

enum
{
  BLOCK_NONE,
  BLOCK1,
  BLOCK2,
  BLOCK3,
  BLOCK4,
};

} // namespace dmtk

#endif // __DMTK_ENUMS_H__
