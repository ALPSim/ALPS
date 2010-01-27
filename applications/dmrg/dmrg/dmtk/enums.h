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
