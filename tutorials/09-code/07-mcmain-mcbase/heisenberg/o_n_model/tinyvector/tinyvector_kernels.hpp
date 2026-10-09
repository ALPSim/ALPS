#ifndef TINYVECTOR_KERNELS_HPP
#define TINYVECTOR_KERNELS_HPP

#include <immintrin.h>
#include "tinyvector_default.hpp"

struct INTRIN_OPT {};

template<int N>
struct optvec{ typedef tinyvector<double, N, INTRIN_OPT> t; };

/** 
 * determine the vectorization step length, aka how many doubles can we fit into a register? 
 */
#ifdef __AVX__
    const int STEP = 4;
#else
    #ifdef __SSE2__
        const int STEP = 2;
    #else
        const int STEP = 1;
    #endif
#endif

/**
 * recursively call binary operator kernels 
 */
template <class Op, int N>
struct _vectorize {
    template <class vec>
    static inline const void apply(vec & left, const vec & right, const int start = 0) {
        _vectorize<Op, STEP>::apply(left, right, start);
        _vectorize<Op, N - STEP>::apply(left, right, start + STEP);
    }
};

template <class Op>
struct _vectorize<Op, 0>
{
    template <class vec>
    static inline const void apply(vec & left, const vec & right, const int start = 0) { }
};

template <class Op>
struct _vectorize<Op, 1>
{
    template <class vec>
    static inline const void apply(vec & left, const vec & right, const int start = 0) {
        Op::apply(left, right, start);
    }
};

#ifdef __SSE2__

template <class Op>
struct _vectorize<Op, 2>
{
    template <class vec>
    static inline const void apply(vec & left, const vec & right, const int start = 0) {
        double * l = left.data();
        const double * r = right.data();
        __m128d mml, mmr, mms;
        mml = _mm_load_pd(l + start);
        mmr = _mm_load_pd(r + start);
        mms = Op::apply(mml, mmr);
        _mm_store_pd(l + start, mms);
    }
};

#endif

#ifdef __AVX__

template <class Op>
struct _vectorize<Op, 3>
{
    template <class vec>
    static inline const void apply(vec & left, const vec & right, const int start = 0) {
        _vectorize<Op, 2>::apply(left, right, start);
        _vectorize<Op, 1>::apply(left, right, start + 2);
    }
};

template <class Op>
struct _vectorize<Op, 4>
{
    template <class vec>
    static inline const void apply(vec & left, const vec & right, const int start = 0) {
        double * l = left.data();
        const double * r = right.data();
        __m256d mml, mmr, mms;
        mml = _mm256_load_pd(l + start);
        mmr = _mm256_load_pd(r + start);
        mms = Op::apply(mml, mmr);
        _mm256_store_pd(l + start, mms);
    }
};

#endif

struct _plus;
struct _minus;
struct _divide;
struct _multiply;
     
template <int N>
inline const tinyvector<double, N, INTRIN_OPT> & operator+= (tinyvector<double, N, INTRIN_OPT> & left, const tinyvector<double, N, INTRIN_OPT> & right) {
    _vectorize<_plus, N>::apply(left, right, 0);
    return left;
}

template <int N>
inline const tinyvector<double, N, INTRIN_OPT> & operator-= (tinyvector<double, N, INTRIN_OPT> & left, const tinyvector<double, N, INTRIN_OPT> & right) {
    _vectorize<_minus, N>::apply(left, right, 0);
    return left;
}

template <int N>
inline const tinyvector<double, N, INTRIN_OPT> & operator*= (tinyvector<double, N, INTRIN_OPT> & left, const tinyvector<double, N, INTRIN_OPT> & right) {
    _vectorize<_multiply, N>::apply(left, right, 0);
    return left;
}

template <int N>
inline const tinyvector<double, N, INTRIN_OPT> & operator/= (tinyvector<double, N, INTRIN_OPT> & left, const tinyvector<double, N, INTRIN_OPT> & right) {
    _vectorize<_divide, N>::apply(left, right, 0);
    return left;
}
 
struct _plus
{
    template <class T>
    static inline const void apply(T & left, const T & right) {
        left += right;
    }
    template <class vec>
    static inline const void apply(vec & left, const vec & right, const int start) {
        left[start] += right[start];
    }
#ifdef __SSE2__
    static inline const __m128d apply(__m128d & mml, __m128d & mmr) {
        return _mm_add_pd(mml, mmr);
    }
#endif
#ifdef __AVX__
    static inline const __m256d apply(__m256d & mml, __m256d & mmr) {
        return _mm256_add_pd(mml, mmr);
    }
#endif
};

struct _minus
{
    template <class T>
    static inline const void apply(T & left, const T & right) {
        left -= right;
    }
    template <class vec>
    static inline const void apply(vec & left, const vec & right, const int start) {
        left[start] -= right[start];
    }
#ifdef __SSE2__
    static inline const __m128d apply(__m128d & mml, __m128d & mmr) {
        return _mm_sub_pd(mml, mmr);
    }
#endif
#ifdef __AVX__
    static inline const __m256d apply(__m256d & mml, __m256d & mmr) {
        return _mm256_sub_pd(mml, mmr);
    }
#endif
};

struct _multiply
{
    template <class T>
    static inline const void apply(T & left, const T & right) {
        left *= right;
    }
    template <class vec>
    static inline const void apply(vec & left, const vec & right, const int start) {
        left[start] *= right[start];
    }
#ifdef __SSE2__
    static inline const __m128d apply(__m128d & mml, __m128d & mmr) {
        return _mm_mul_pd(mml, mmr);
    }
#endif
#ifdef __AVX__
    static inline const __m256d apply(__m256d & mml, __m256d & mmr) {
        return _mm256_mul_pd(mml, mmr);
    }
#endif
};

struct _divide
{
    template <class T>
    static inline const void apply(T & left, const T & right) {
        left /= right;
    }
    template <class vec>
    static inline const void apply(vec & left, const vec & right, const int start) {
        left[start] /= right[start];
    }
#ifdef __SSE2__
    static inline const __m128d apply(__m128d & mml, __m128d & mmr) {
        return _mm_div_pd(mml, mmr);
    }
#endif
#ifdef __AVX__
    static inline const __m256d apply(__m256d & mml, __m256d & mmr) {
        return _mm256_div_pd(mml, mmr);
    }
#endif
};

#endif
