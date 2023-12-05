//
// Created by haralg on 05.12.23.
//

#include "complex_util.hpp"
#include <complex.h>
#include <cstdio>
#include <chrono>
#include <Eigen/Dense>

#define TIME_INIT std::chrono::high_resolution_clock::time_point t1, t2;
#define TIME_START t1 = std::chrono::high_resolution_clock::now();
#define TIME_END t2 = std::chrono::high_resolution_clock::now(); \
    printf("Time: %f\n", std::chrono::duration<double>(t2 - t1).count()); \



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

int test_csqrt(int N) {
    // Init C-style complex array
    TIME_INIT;
    printf("csqrt: init array\n");

    _Complex double A[N];
    TIME_START;
    for (int i = 0; i < N; ++i) {
        double a = (double) i;
        double b = (double) i;
        A[i] = a + b * I;
    }
    TIME_END;

    printf("csqrt: test csqrt\n");
    _Complex double A_res[N];
    TIME_START;
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        A_res[i] = csqrt(A[i]);
    }
    TIME_END;

    return 0;
}

int test_avx(int N) {
    // Init AVXComplex array
    TIME_INIT;
    printf("avx: init array\n");

    AVXComplex B[N / 4];
    TIME_START;
    for (int i = 0; i < N; i += 4) {

        double b1 = (double)(i);
        double b2 = (double)(i + 1);
        double b3 = (double)(i + 2);
        double b4 = (double)(i + 3);
        B[i / 4] = AVXComplex
                {
                    _mm256_set_pd(b1, b2, b3, b4),
                    _mm256_set_pd(b1, b2, b3, b4)
                };
    }
    TIME_END;

    printf("avx: test mm_c_sqrt\n");
    AVXComplex B_res[N / 4];
    TIME_START;
#pragma omp parallel for
    for (int i = 0; i < N / 4; ++i) {
        B_res[i] = mm_c_sqrt_unsafe(B[i]);
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
