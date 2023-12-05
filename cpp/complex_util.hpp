//
// Created by haralg on 03.11.23.
//
#ifndef CPP_COMPLEX_UTIL_HPP
#define CPP_COMPLEX_UTIL_HPP

#include <immintrin.h>

typedef double real_t;

struct AVXComplex {
    __m256d real;
    __m256d imag;
};

#define MM_CABS(re,im) _mm256_sqrt_pd(_mm256_fmadd_pd(re, re, _mm256_mul_pd(im, im)))


static inline AVXComplex mm_c_sqrt_unsafe(AVXComplex z) {
    // Vectorized principal complex square root
    // If real part is negative
    // Imaginary part of A is always positive in our case

    __m256d r = MM_CABS(z.real, z.imag);
    __m256d half = _mm256_set1_pd(0.5);

    return AVXComplex{
        _mm256_sqrt_pd(half*(r + z.real)),
        _mm256_sqrt_pd(half*(r - z.real))
    };
}

static inline AVXComplex alpha(__mm256d q, AVXComplex a) {
    // Vectorized alpha function
    // alpha(q, a) = sqrt(a - q^2)
    // Assumes that imag(a) >= 0
    return mm_c_sqrt_unsafe(AVXComplex{
        _mm256_fnmadd_pd(q, q, a.real),
        a.imag
    });
}

#endif //CPP_COMPLEX_UTIL_HPP
