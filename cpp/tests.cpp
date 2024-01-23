//
// Created by haralg on 05.12.23.
//

#include "complex_util.hpp"
#include <complex.h>
#include <cstdio>
//#include <cstdlib>
#include <chrono>
#include <Eigen/Dense>
#include <omp.h>

#define TIME_INIT std::chrono::high_resolution_clock::time_point t1, t2
#define TIME_START t1 = std::chrono::high_resolution_clock::now()
#define TIME_END t2 = std::chrono::high_resolution_clock::now(); \
    printf("Time: %f\n", std::chrono::duration<double>(t2 - t1).count())



int test_csqrt(int N);
int test_avx(int N);
int test_eigen(int N);

int main(int argc, char **argv) {
    TIME_INIT;
    int N = atoi(argv[1]);
    char choice = argv[2][0];

    switch (choice) {
        case 'a':
            test_csqrt(N);
            test_avx(N);
            break;
        case 'e':
            test_eigen(N);
            break;
        default:
            printf("Invalid choice\n");
            break;
    }



    return 0;
}

double drand() {
    // Makes random double in (0.01, 1.00)
    return 0.01*(double)(rand() % 100 + 1);
}

int test_csqrt(int N) {
    // C-style complex array sqrt
    TIME_INIT;

    printf("csqrt: alloc\n");
    TIME_START;
    _Complex double input[N];
    _Complex double res[N];
    for (int i = 0; i < N; ++i) {
        input[i] = drand() + drand() * I;
    }
    TIME_END;

    printf("csqrt: test csqrt\n");
    TIME_START;
#pragma omp parallel for simd
    for (int i = 0; i < N; ++i) {
        res[i] = csqrt(input[i]);
    }
    TIME_END;

    return 0;
}

int test_avx(int N) {
    // AVXComplex array square root
    TIME_INIT;

    printf("avx: alloc\n");
    TIME_START;
    AVXComplex input[N / 4];
    AVXComplex res[N / 4];
    for (int i = 0; i < N / 4; ++i) {
        input[i].real = _mm256_set_pd(drand(), drand(), drand(), drand());
        input[i].imag = _mm256_set_pd(drand(), drand(), drand(), drand());
    }
    TIME_END;

    printf("avx: test mm_csqrt_unsafe\n");
    TIME_START;
//#pragma omp parallel for shared(res, input)
    for (int i = 0; i < N / 4; ++i) {
        res[i] = mm_csqrt_unsafe(input[i]);
    }
    TIME_END;

    return 0;

}

int test_eigen(int N) {
    using namespace Eigen;
    // Test Eigen with C-style complex
    Matrix<_Complex double, Dynamic, Dynamic> C;

    for (int i = 0; i < N; ++i) {
        double a = (double) i;
        double b = (double) i;
        C(i, 0) = a + b * I;
    }


    return 0;
}
