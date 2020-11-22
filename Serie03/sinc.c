#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <immintrin.h>
#include <math.h>
#include <omp.h>

#include "stopwatch.h"

// Typedef for 'real'. Just another name for double.
typedef double real;

// Compute max. abs. error between array 'x' and 'y' of length 'n'
real max_abs_error(real *x, real *y, uint n)
{
    uint i;
    real error;

    real max_error;

    max_error = fabs(x[0] - y[0]);
    for (i = 1; i < n; i++)
    {
        error = fabs(x[i] - y[i]);
        max_error = error > max_error ? error : max_error;
    }

    return max_error;
}

// Compute y = sinc(x) for arrays 'x' and 'y' of length 'n'.
void libm_sinc(uint n, real *x, real *y)
{
    uint i;
    real xi;

    for (i = 0; i < n; i++)
    {
        xi = x[i];
        y[i] = sin(xi) / xi;
    }
}

// Compute y = sinc(x) for arrays 'x' and 'y' of length 'n' using
// taylor expansion.
void taylor_sinc(uint n, real *x, real *y)
{
    //TODO: Implement the evalutation of the array x with the 'sinc'
    //      function in scalar arithmetics.
}

// Compute y = sinc(x) for arrays 'x' and 'y' of length 'n' using
// taylor expansion.
void taylor_sinc_avx512(uint n, real *x, real *y)
{
    //TODO: Implement the evalutation of the array x with the 'sinc'
    //      function in either AVX or AVX512 arithmetics.
}

// Compute y = sinc(x)+sinc(2x) for arrays 'x' and 'y' of length 'n'.
void libm_sincsum(uint n, real *x, real *y)
{
    uint i;
    real xi, xi2;

    for (i = 0; i < n; i++)
    {
        xi = x[i];
        xi2 = 2.0 * xi;
        y[i] = (sin(xi) / xi) + (sin(xi2) / xi2);
    }
}

// Compute y = sinc(x) for arrays 'x' and 'y' of length 'n' using
// taylor expansion.
void taylor_sincsum(uint n, real *x, real *y)
{
    //TODO: Implement the evalutation of the array x with the 'sinc(x) + sinc(2x)'
    //      function in scalar arithmetics.
}

// Compute y = sinc(x)+sinc(2x) for arrays 'x' and 'y' of length 'n' using
// taylor expansion.
void taylor_sincsum_avx512(uint n, real *x, real *y)
{
    //TODO: Implement the evalutation of the array x with the 'sinc(x) + sinc(2x)'
    //      function in either AVX or AVX512 arithmetics.
}

int main(int argc, char const *argv[])
{

    real *x, *y1, *y2;
    uint n, iter;
    real range;
    uint i, j;
    real scale;
    pstopwatch sw;
    real t, t2;
    real error;

    (void)argc;
    (void)argv;

    n = 1 << 13;
    iter = 1 << 13;
    range = 1.0 * M_PI;
    sw = new_stopwatch();

    x = (real *)_mm_malloc(n * sizeof(real), 64);
    y1 = (real *)_mm_malloc(n * sizeof(real), 64);
    y2 = (real *)_mm_malloc(n * sizeof(real), 64);

    // fill test values

    scale = 1.0 / RAND_MAX;
    for (i = 0; i < n; i++)
    {
        x[i] = 2.0 * range * (((real)rand() * scale) - 0.5);
    }

    printf("Computing function on interval [%+.2f, %+.2f]:\n", -range, range);

    printf("================================\n"
           "             sinc(x)            \n"
           "================================\n");

    printf("Computing sinc(x) by libM:\n");
    /* Cache warmup */
    for (j = 0; j < iter; j++)
    {
        libm_sinc(n, x, y1);
    }
    start_stopwatch(sw);
    for (j = 0; j < iter; j++)
    {
        libm_sinc(n, x, y1);
    }
    t = stop_stopwatch(sw);
    printf("  %.3f ms (%.3f)\n", t * 1.0e3, y1[15]);

    printf("Computing sinc(x) by Taylor:\n");
    /* Cache warmup */
    for (j = 0; j < iter; j++)
    {
        taylor_sinc(n, x, y2);
    }
    start_stopwatch(sw);
    for (j = 0; j < iter; j++)
    {
        taylor_sinc(n, x, y2);
    }
    t2 = stop_stopwatch(sw);
    printf("  %.3f ms (%.3f)\n", t2 * 1.0e3, y2[15]);
    printf("max error:\n");
    error = max_abs_error(y1, y2, n);
    printf("  %.5e\n", error);
    printf("Speedup:\n");
    printf("  %.2f\n", t / t2);

    printf("Computing sinc(x) by Taylor using AVX512:\n");
    /* Cache warmup */
    for (j = 0; j < iter; j++)
    {
        taylor_sinc_avx512(n, x, y2);
    }
    start_stopwatch(sw);
    for (j = 0; j < iter; j++)
    {
        taylor_sinc_avx512(n, x, y2);
    }
    t2 = stop_stopwatch(sw);
    printf("  %.3f ms (%.3f)\n", t2 * 1.0e3, y2[15]);
    printf("max error:\n");
    error = max_abs_error(y1, y2, n);
    printf("  %.5e\n", error);
    printf("Speedup:\n");
    printf("  %.2f\n", t / t2);

    printf("================================\n"
           "        sinc(x) + sinc(2x)      \n"
           "================================\n");

    printf("Computing sinc(x)+sinc(2x) by libM:\n");
    /* Cache warmup */
    for (j = 0; j < iter; j++)
    {
        libm_sincsum(n, x, y1);
    }
    start_stopwatch(sw);
    for (j = 0; j < iter; j++)
    {
        libm_sincsum(n, x, y1);
    }
    t = stop_stopwatch(sw);
    printf("  %.3f ms (%.3f)\n", t * 1.0e3, y1[15]);

    printf("Computing sinc(x)+sinc(2x) by Taylor:\n");
    /* Cache warmup */
    for (j = 0; j < iter; j++)
    {
        taylor_sincsum(n, x, y2);
    }
    start_stopwatch(sw);
    for (j = 0; j < iter; j++)
    {
        taylor_sincsum(n, x, y2);
    }
    t2 = stop_stopwatch(sw);
    printf("  %.3f ms (%.3f)\n", t2 * 1.0e3, y2[15]);
    printf("max error:\n");
    error = max_abs_error(y1, y2, n);
    printf("  %.5e\n", error);
    printf("Speedup:\n");
    printf("  %.2f\n", t / t2);

    printf("Computing sinc(x) by Taylor using AVX512:\n");
    /* Cache warmup */
    for (j = 0; j < iter; j++)
    {
        taylor_sincsum_avx512(n, x, y2);
    }
    start_stopwatch(sw);
    for (j = 0; j < iter; j++)
    {
        taylor_sincsum_avx512(n, x, y2);
    }
    t2 = stop_stopwatch(sw);
    printf("  %.3f ms (%.3f)\n", t2 * 1.0e3, y2[15]);
    printf("max error:\n");
    error = max_abs_error(y1, y2, n);
    printf("  %.5e\n", error);
    printf("Speedup:\n");
    printf("  %.2f\n", t / t2);

    _mm_free(x);
    _mm_free(y1);
    _mm_free(y2);

    return 0;
}
