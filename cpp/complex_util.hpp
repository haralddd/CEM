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

AVXComplex mm_c_sqrt(AVXComplex z) {
    // Vectorized complex square root
    // Only works for positive real part, handle other cases separately
    __m256d r = MM_CABS(z.real, z.imag);
    __m256d a = _mm256_sqrt_pd(r);
    __m256d b_re = z.real + r;
    __m256d b_abs = MM_CABS(b_re, z.imag);
    __m256d ret_re = a * b_re / b_abs;
    __m256d ret_im = a * z.imag / b_abs;

    return AVXComplex{ret_re, ret_im};
}

#endif //CPP_COMPLEX_UTIL_HPP
