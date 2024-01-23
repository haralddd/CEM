//
// Created by haralg on 03.11.23.
//
#ifndef CPP_COMPLEX_UTIL_HPP
#define CPP_COMPLEX_UTIL_HPP

#include <immintrin.h>

typedef double real_t;
typedef __m256d m_real_t; // Real type AVX vector
template<unsigned N> struct avx_pow_impl;

struct AVXComplex {
    m_real_t real;
    m_real_t imag;

    AVXComplex squared() {
        return {
            _mm256_fnmadd_pd(imag, imag, real*real),
            _mm256_set1_pd(2.0)*real*imag
        };
    }

    AVXComplex operator*(const AVXComplex &other) {
        return AVXComplex{
            _mm256_fnmadd_pd(imag, other.imag, real*other.real),
            _mm256_fmadd_pd(real, other.imag, imag*other.real)
        };
    }

    static AVXComplex unity() {
        return {
            _mm256_set1_pd(1.0),
            _mm256_set1_pd(0.0)
        };
    }
};

#define MM_CABS_UNSAFE(re,im) _mm256_sqrt_pd(_mm256_fmadd_pd(re, re, _mm256_mul_pd(im,im)))

/*
 * Vectorized principal complex square root
 * UNSAFE:
    * Not overflow safe
    * All results in first quadrant, result(re >= 0, im >= 0)
    * Omitted sgn(z.imag) in result.imag, assumed z.imag > 0
 * Thus: Use this when z has real and imag coeffs close to unity
 */
static inline AVXComplex mm_csqrt_unsafe(AVXComplex z) {

    m_real_t r = MM_CABS_UNSAFE(z.real, z.imag);
    m_real_t half = _mm256_set1_pd(0.5);

    return AVXComplex{
        _mm256_sqrt_pd(half*(r + z.real)),
        _mm256_sqrt_pd(half*(r - z.real))
    };
}

template<unsigned N> struct avx_power_impl {
    static AVXComplex calc(const AVXComplex &x) {
        if (N%2 == 0)
            return power_impl<N/2>::calc(x*x);
        else if (N%3 == 0)
            return power_impl<N/3>::calc(x*x*x);
        return power_impl<N-1>::calc(x)*x;
    }
};

template<> struct avx_power_impl<0> {
    static AVXComplex calc(const AVXComplex &) { return AVXComplex::unity(); }
};

template<> struct avx_power_impl<1> {
    static AVXComplex calc(const AVXComplex &x) { return x; }
};

template<> struct avx_power_impl<2> {
    static AVXComplex calc(const AVXComplex &x) { return x*x; }
};

template<> struct avx_power_impl<3> {
    static AVXComplex calc(const AVXComplex &x) { return x*x*x; }
};

template<> struct avx_power_impl<4> {
    static AVXComplex calc(const AVXComplex &x) {
        AVXComplex x2 = x*x;
        return x2*x2;
    }
};

template<> struct avx_power_impl<5> {
    static AVXComplex calc(const AVXComplex &x) {
        AVXComplex x2 = x*x;
        return x2*x2*x;
    }
};

template<> struct avx_power_impl<6> {
    static AVXComplex calc(const AVXComplex &x) {
        AVXComplex x2 = x*x;
        return x2*x2*x2;
    }
};

template<> struct avx_power_impl<7> {
    static AVXComplex calc(const AVXComplex &x) {
        AVXComplex x2 = x*x;
        AVXComplex x5 = x2*x2*x;
        return x5*x2;
    }
};


template<unsigned N>
inline AVXComplex power(const AVXComplex &x) {
    return avx_power_impl<N>::calc(x);
}


static inline AVXComplex mm_alpha_unsafe(m_real_t q, AVXComplex a) {
    // Vectorized alpha function
    // alpha(q, a) = sqrt(a - q^2)
    // Assumes that imag(a) >= 0
    return mm_csqrt_unsafe(AVXComplex{
        _mm256_fnmadd_pd(q, q, a.real),
        a.imag
    });
}


static inline AVXComplex M_ker(m_real_t p, m_real_t q, AVXComplex a, AVXComplex a0, int n) {
    // Vectorized M_kernel function
    AVXComplex da_neg = mm_sub(a, a0);
    return mm_csqrt_unsafe(AVXComplex{
        _mm256_fnmadd_pd(a.imag, a.imag, q),
        a.imag
    });
}
*/
#endif //CPP_COMPLEX_UTIL_HPP
