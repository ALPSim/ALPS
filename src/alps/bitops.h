/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1999-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Andreas Laeuchli <laeuchli@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

/* $Id$ */

#ifndef ALPS_SRC_ALPS_BITOPS_H
#define ALPS_SRC_ALPS_BITOPS_H
 
#include <alps/config.h>

#ifdef cray
# include <intrinsics.h>
#endif

namespace alps {

//
// Cray intrinsic bit operations : gbit, gbits, maskr, popcnt
//

#ifdef cray
# define gbit   _gbit
# define gbits  _gbits
# define maskr  _maskr
# define popcnt _popcnt
#else

/** extract a bit from a word.
    @param x the word
    @param n position of the bit to be extracted
    @return the n-th bit of the word x */
template <class T, class N>
inline T gbit(T x, N n) { return (x>>n)&1; }

/** extract bits from a word.
    @param x the word
    @param m the number of bits to be extracted
    @param n position of the bit to be extracted
    @return the m bits starting at bit n  */
template <class T, class N>
inline T gbits(T x, N m, long n) { return (x>>n)&((1<<m)-1); }

/** create a right-justified N
    creates a bitpattern with the rightmost i bits set
    @param i the number of bits to be set to 1
    @return a word with the rightmost bits set to 1 */
inline uint32_t maskr(uint16_t i) {return (1u<<i)-1;}
inline static int BX_(long x) { return ((x) - (((x)>>1)&0x77777777)
                             - (((x)>>2)&0x33333333)
                             - (((x)>>3)&0x11111111)); }

/// count the 1-bits in a word
inline long popcnt(uint32_t x) 
{ return (((BX_(x)+(BX_(x)>>4)) & 0x0F0F0F0F) % 255); }

#endif

//
// Fortran-style bit operations : btest, ibits, ibclr, ibset
//

// btest : check p-th bit of i (return value : true for 1 and false for 0)
template <class T, class U>
inline bool btest(T i, U p) { return gbit(i, p) == 1; }

// ibits : extract p-th bit from i
template <class T, class U>
inline T ibits(T i, U p) { return gbit(i, p); }

// ibits : extract n bits [p,p+n-1] from i
template <class T, class U, class V>
inline T ibits(T i, U p, V n) { return gbits(i, p, n); }

// ibclr : clear p-th bit of i to 0
template <class T, class U>
inline T ibclr(T i, U p) { return i & (~(1 << p)); }

// ibset : set p-th bit of i to 1
template <class T, class U>
inline T ibset(T i, U p) { return i | (1 << p); }

// ibset : set p-th bit of i to b (0 or 1)
template <class T, class U, class V>
inline T ibset(T i, U p, V b) { return i & (~(1 << p)) | ((b & 1) << p); }

} // end namespace alps

#endif // ALPS_SRC_ALPS_BITOPS_H
