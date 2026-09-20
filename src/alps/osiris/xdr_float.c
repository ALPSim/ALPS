/*
 * xdr_float.c, Generic XDR routines implementation.
 *
 * Copyright (c) 2010, Oracle America, Inc.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *     * Neither the name of the "Oracle America, Inc." nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *   COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 *   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *   GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 *   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * These are the "floating point" xdr routines used to (de)serialize
 * most common data items.  See xdr.h for more info on the interface to
 * xdr.
 */

#include "xdrcore.h"
#include <float.h>
#include <string.h>

bool_t xdr_int32_t(XDR* stream, int32_t* value) {
    switch (stream->x_op) {
    case XDR_ENCODE: return XDR_PUTINT32(stream, value);
    case XDR_DECODE: return XDR_GETINT32(stream, value);
    case XDR_FREE: return TRUE;
    }
    return FALSE;
}
bool_t xdr_uint32_t(XDR* stream, uint32_t* value) {
    int32_t bits = 0;
    if (stream->x_op == XDR_ENCODE) memcpy(&bits, value, sizeof(bits));
    if (!xdr_int32_t(stream, &bits)) return FALSE;
    if (stream->x_op == XDR_DECODE) memcpy(value, &bits, sizeof(bits));
    return TRUE;
}

/* IEEE-754 bits are serialized through fixed-width XDR integers. No type
   punning, host word-order macros, or obsolete VAX representations. */
bool_t xdr_float(XDR* stream, float* value) {
    uint32_t bits = 0;
    _Static_assert(sizeof(float) == 4 && FLT_RADIX == 2 && FLT_MANT_DIG == 24,
                   "XDR requires IEEE-754 binary32");
    if (stream->x_op == XDR_ENCODE) memcpy(&bits, value, sizeof(bits));
    if (!xdr_uint32_t(stream, &bits)) return FALSE;
    if (stream->x_op == XDR_DECODE) memcpy(value, &bits, sizeof(bits));
    return TRUE;
}
bool_t xdr_double(XDR* stream, double* value) {
    uint64_t bits = 0;
    _Static_assert(sizeof(double) == 8 && FLT_RADIX == 2 && DBL_MANT_DIG == 53,
                   "XDR requires IEEE-754 binary64");
    if (stream->x_op == XDR_ENCODE) memcpy(&bits, value, sizeof(bits));
    uint32_t high = (uint32_t)(bits >> 32), low = (uint32_t)bits;
    if (!xdr_uint32_t(stream, &high) || !xdr_uint32_t(stream, &low)) return FALSE;
    if (stream->x_op == XDR_DECODE) {
        bits = ((uint64_t)high << 32) | low;
        memcpy(value, &bits, sizeof(bits));
    }
    return TRUE;
}
