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
  int k, i;
  int m = 11; // needs to be adjusted to get the needed precission
  real positiveSummand, negativeSummand;
  uint denominator;
  real numerator;
  real x_value;
  // for (i = 0; i < n; i++)
  {
    x_value = x[i];
    numerator = x_value * x_value * x_value; // x^{2*1+1} for k = 1
    denominator = 6; // (2+1)! for k = 1
    positiveSummand = x_value;
    negativeSummand = numerator / denominator;
    printf("=========== \n %.5e\n", x_value);
    // k = 1;
    for (k = 1; k < m-1; k = k + 2)
    {
    printf("Step: %i, PosSum: %.5e, NegSum: %.5e\n", k, positiveSummand, negativeSummand);
    printf("Zähler: %.5e, Nenner: %.5d, x: %.5e\n", numerator, denominator, x_value);
      // x^{2k+1} = x^{2k}*x = x^{2k-2+2}*x = x^{2(k-1) + 2}*x = x^{2k-1}*x^2*x
      numerator *= x_value * x_value * x_value;
      // (2(k+1)+1)! = (2k+2+1)! = (2k+3)! = (2k+1)! * (2k+3) * (2k+2) = (2k+1)! * (4k^2 + 10k + 6) = (2k+1)! * (k(4k + 10) + 6)
      denominator *=  (k * (4*k + 10) + 6);
      negativeSummand += numerator / denominator;
    printf("HALFSTEP Zähler: %.5e, Nenner: %.5d, Summand: %.5e\n", numerator, denominator, negativeSummand);

      numerator *= x_value * x_value * x_value;
      denominator *= (k*(4*k+18)+20); // *= (k^2 + 9k + 20)*= (k+4)*(k+5)
      positiveSummand += numerator / denominator;
      printf("FULLSTEP Zähler: %.5e, Nenner: %.5d, Summand: %.5e\n", numerator, denominator, positiveSummand);
      printf("------\n");
    }
    if (k < m)
    {
      numerator *= x_value * x_value * x_value;
      denominator *= k*4*(k-1);
      negativeSummand += numerator / denominator;
    }
  }

  y[i] = positiveSummand - negativeSummand;
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

    // n = 1 << 13;
    n = 16;
    iter = 1;
    // iter = 1 << 13;
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

    // printf("================================\n"
    //        "        sinc(x) + sinc(2x)      \n"
    //        "================================\n");

    // printf("Computing sinc(x)+sinc(2x) by libM:\n");
    // /* Cache warmup */
    // for (j = 0; j < iter; j++)
    // {
    //     libm_sincsum(n, x, y1);
    // }
    // start_stopwatch(sw);
    // for (j = 0; j < iter; j++)
    // {
    //     libm_sincsum(n, x, y1);
    // }
    // t = stop_stopwatch(sw);
    // printf("  %.3f ms (%.3f)\n", t * 1.0e3, y1[15]);

    // printf("Computing sinc(x)+sinc(2x) by Taylor:\n");
    // /* Cache warmup */
    // for (j = 0; j < iter; j++)
    // {
    //     taylor_sincsum(n, x, y2);
    // }
    // start_stopwatch(sw);
    // for (j = 0; j < iter; j++)
    // {
    //     taylor_sincsum(n, x, y2);
    // }
    // t2 = stop_stopwatch(sw);
    // printf("  %.3f ms (%.3f)\n", t2 * 1.0e3, y2[15]);
    // printf("max error:\n");
    // error = max_abs_error(y1, y2, n);
    // printf("  %.5e\n", error);
    // printf("Speedup:\n");
    // printf("  %.2f\n", t / t2);

    // printf("Computing sinc(x) by Taylor using AVX512:\n");
    // /* Cache warmup */
    // for (j = 0; j < iter; j++)
    // {
    //     taylor_sincsum_avx512(n, x, y2);
    // }
    // start_stopwatch(sw);
    // for (j = 0; j < iter; j++)
    // {
    //     taylor_sincsum_avx512(n, x, y2);
    // }
    // t2 = stop_stopwatch(sw);
    // printf("  %.3f ms (%.3f)\n", t2 * 1.0e3, y2[15]);
    // printf("max error:\n");
    // error = max_abs_error(y1, y2, n);
    // printf("  %.5e\n", error);
    // printf("Speedup:\n");
    // printf("  %.2f\n", t / t2);

    _mm_free(x);
    _mm_free(y1);
    _mm_free(y2);

    return 0;
}
