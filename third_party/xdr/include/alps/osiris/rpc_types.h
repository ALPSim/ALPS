/*
 * Copyright (c) 2010, Oracle America, Inc.
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
 */
/* fixincludes should not add extern "C" to this file */
/*
 * Rpc additions to <sys/types.h>
 */
#ifndef ALPS_OSIRIS_RPC_TYPES_H
#define ALPS_OSIRIS_RPC_TYPES_H
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
typedef int bool_t;
typedef int enum_t;
typedef unsigned char alps_xdr_uchar;
typedef unsigned short alps_xdr_ushort;
typedef unsigned int alps_xdr_uint;
typedef unsigned long alps_xdr_ulong;
typedef int64_t alps_xdr_int64;
typedef uint64_t alps_xdr_uint64;
typedef char* alps_xdr_address;
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif
#define __dontcare__ -1
#define mem_alloc(size) malloc(size)
#define mem_free(ptr, size) free(ptr)
/* A network word is big endian regardless of host byte order. */
static inline uint32_t alps_xdr_network_word(uint32_t word) {
    const uint16_t marker = 1;
    unsigned char first;
    memcpy(&first, &marker, 1);
    if (!first) return word;
    return ((word & UINT32_C(0xff)) << 24) | ((word & UINT32_C(0xff00)) << 8)
        | ((word >> 8) & UINT32_C(0xff00)) | (word >> 24);
}
#endif
