/* 
 * Copyright Matthias Troyer 2005
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
 */

#include <boost/preprocessor/cat.hpp>

#ifndef ALPS_SPRNG_CALL
#error Do not include this file directly
#else

// This header file declares the C API functions for the SPRNG generators

extern "C" {

int ALPS_SPRNG_CALL(get_rn_int) (int *igenptr);
float  ALPS_SPRNG_CALL(get_rn_flt) (int *igenptr);
double  ALPS_SPRNG_CALL(get_rn_dbl) (int *igenptr);
int *  ALPS_SPRNG_CALL(init_rng) (int rng_type,  int gennum, int total_gen,  int seed, int mult);
int  ALPS_SPRNG_CALL(spawn_rng) (int *igenptr, int nspawned, int ***newgens, int checkid);
int  ALPS_SPRNG_CALL(get_seed_rng) (int *genptr);
int  ALPS_SPRNG_CALL(free_rng) (int *genptr);
int  ALPS_SPRNG_CALL(pack_rng) ( int *genptr, char **buffer);
int * ALPS_SPRNG_CALL(unpack_rng) ( char const *packed);
int  ALPS_SPRNG_CALL(print_rng) ( int *igen);

}

#endif 
