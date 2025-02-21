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

#ifndef __DMTK_BITS_H__
#define __DMTK_BITS_H__

namespace dmtk
{

#define IS_EVEN(n) ((n & 1) == 1 ? false : true)
#define SGN(n) (1 - (((n) & 1) << 1))
#define IEOR(n1,n2) ((n1) ^ (n2))
#define IAND(n1,n2) ((n1) & (1 << (n2)))
#define IOR(n1,n2) ((n1) | (n2))
#define IBITS(n,i) (((n) & (1 << i)) >> i)
#define IBSET(n,i) ((n) | (1 << i))
#define IBCLR(n,i) ((n) ^ (IBITS((n),i) << i))
inline int ISHFTC(int n, int i, int nt)
{
  int mask = (1 << i) - 1;
  int tail = n & mask;
  int r = n >> i;
  r |= tail << (nt - i);
  return r;
}

inline size_t mask(int b1, int b2) { return (1 << (b1-1) | 1 << (b2-1)); }
inline size_t mask(int b) { return (1 << (b-1)); }
inline int mod2(int n) { return (n & 1); }

} // namespace dmtk

#endif // __DMTK_BITS_H__
