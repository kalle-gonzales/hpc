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
    //uint i;
    //real xi;
    //real xisq;

    //for (i = 0; i < n; i++)
    //{
    //    xi = x[i];
    //    xisq = xi * xi;
    //    y[i] = (xi * (1.0 - ((xisq) / (6.0)) * (1.0 - ((xisq) / (20.0))) * (1.0 - ((xisq) / (42.0)) * (1.0 - ((xisq) / (72.0)) * (1.0 - ((xisq) / (110.0))))))) / (xi);
    //}

    uint i, j, m = 11;
    real xi, xisq;
    real sin;

    for (i = 0; i < n; i++)
    {
        xi = x[i];
        xisq = xi * xi;
        sin = 1;
        for (j = m; j - 1 > 0; j -= 2)
        {
            sin *= 1 - xisq / (j * (j - 1));
        }

        // normally, we would have to multiply sin with xi, to finish the evaluation of the horner schema. But as we divide sin by xi anyways, we are good as it is.
        y[i] = sin;
    }
}

const __m256d v_one = {1.0, 1.0, 1.0, 1.0};
const __m256d v_two = {2.0, 2.0, 2.0, 2.0};
const __m256d v_six = {6.0, 6.0, 6.0, 6.0};
const __m256d v_twenty = {20.0, 20.0, 20.0, 20.0};
const __m256d v_fortyTwo = {42.0, 42.0, 42.0, 42.0};
const __m256d v_seventyTwo = {72.0, 72.0, 72.0, 72.0};
const __m256d v_oneHundredTen = {110.0, 110.0, 110.0, 110.0};

// Compute y = sinc(x) for arrays 'x' and 'y' of length 'n' using
// taylor expansion.
// we use _m256d to contain 4 doubles. With _m256 we could contain 8 floats, but reals are doubles.
void taylor_sinc_avx512(uint n, real *x, real *y)
{
    __m256d v_y, v_xi, v_xisq, v_sin;

    uint i;

    for (i = 0; i + 3 < n; i += 4)
    {
        v_y = _mm256_loadu_pd(y + i);
        v_xi = _mm256_loadu_pd(x + i);
        v_xisq = _mm256_mul_pd(v_xi, v_xi);
        v_sin = _mm256_mul_pd(v_xi, _mm256_mul_pd(_mm256_sub_pd(v_one, _mm256_div_pd(v_xisq, v_six)), _mm256_mul_pd(_mm256_sub_pd(v_one, _mm256_div_pd(v_xisq, v_twenty)), _mm256_mul_pd(_mm256_sub_pd(v_one, _mm256_div_pd(v_xisq, v_fortyTwo)), _mm256_mul_pd(_mm256_sub_pd(v_one, _mm256_div_pd(v_xisq, v_seventyTwo)), _mm256_sub_pd(v_one, _mm256_div_pd(v_xisq, v_oneHundredTen)))))));
        v_y = _mm256_div_pd(v_sin, v_xi);
        _mm256_storeu_pd(y + i, v_y);
    }

    real xi, xisq;
    for (; i < n; i++)
    {
        xi = x[i];
        xisq = xi * xi;
        y[i] = (xi * (1.0 - ((xisq) / (6.0)) * (1.0 - ((xisq) / (20.0))) * (1.0 - ((xisq) / (42.0)) * (1.0 - ((xisq) / (72.0)) * (1.0 - ((xisq) / (110.0))))))) / (xi);
    }
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
    //uint i;
    //real xi, xi2;
    //real xisq, xi2sq;
    //
    //for (i = 0; i < n; i++)
    //{
    //    xi = x[i];
    //    xisq = xi * xi;
    //    xi2 = 2.0 * xi;
    //    xi2sq = xi2 * xi2;
    //
    //    y[i] = ((xi * (1.0 - ((xisq) / (6.0)) * (1.0 - ((xisq) / (20.0))) * (1.0 - ((xisq) / (42.0)) * (1.0 - ((xisq) / (72.0)) * (1.0 - ((xisq) / (110.0))))))) / (xi)) + ((xi2 * (1.0 - ((xi2sq) / (6.0)) * (1.0 - ((xi2sq) / (20.0))) * (1.0 - ((xi2sq) / (42.0)) * (1.0 - ((xi2sq) / (72.0)) * (1.0 - ((xi2sq) / (110.0))))))) / (xi2));
    //}

    uint i,j,m=11;
    real xi,xi2, xisq, xi2sq;
    real sinxi, sinxi2;

    for(i = 0; i < n; i++) {
        xi     = x[i];
        xi2    = 2.0 * xi;
        xisq   = xi * xi;
        xi2sq  = xi2 * xi2;
        sinxi  = 1;
        sinxi2 = 1;
        for(j = m; j - 1 > 0; j -= 2) {
            sinxi *= 1 - xisq / (j * (j-1));
            sinxi2 *= 1 - xi2sq / (j * (j-1));
        }
        y[i] = sinxi + sinxi2;
    }
}

// Compute y = sinc(x)+sinc(2x) for arrays 'x' and 'y' of length 'n' using
// taylor expansion.
void taylor_sincsum_avx512(uint n, real *x, real *y)
{
    __m256d v_y, v_xi, v_xi2, v_xisq, v_xi2sq, v_sinxi, v_sinxi2, v_sincxi, v_sincxi2;

    uint i;

    for (i = 0; i + 3 < n; i += 4)
    {
        v_y = _mm256_loadu_pd(y + 1);
        v_xi = _mm256_loadu_pd(x + i);
        v_xi2 = _mm256_mul_pd(v_two, v_xi);
        v_xisq = _mm256_mul_pd(v_xi, v_xi);
        v_xi2sq = _mm256_mul_pd(v_xi2, v_xi2);
        v_sinxi = _mm256_mul_pd(v_xi, _mm256_mul_pd(_mm256_sub_pd(v_one, _mm256_div_pd(v_xisq, v_six)), _mm256_mul_pd(_mm256_sub_pd(v_one, _mm256_div_pd(v_xisq, v_twenty)), _mm256_mul_pd(_mm256_sub_pd(v_one, _mm256_div_pd(v_xisq, v_fortyTwo)), _mm256_mul_pd(_mm256_sub_pd(v_one, _mm256_div_pd(v_xisq, v_seventyTwo)), _mm256_sub_pd(v_one, _mm256_div_pd(v_xisq, v_oneHundredTen)))))));
        v_sinxi2 = _mm256_mul_pd(v_xi2, _mm256_mul_pd(_mm256_sub_pd(v_one, _mm256_div_pd(v_xi2sq, v_six)), _mm256_mul_pd(_mm256_sub_pd(v_one, _mm256_div_pd(v_xi2sq, v_twenty)), _mm256_mul_pd(_mm256_sub_pd(v_one, _mm256_div_pd(v_xi2sq, v_fortyTwo)), _mm256_mul_pd(_mm256_sub_pd(v_one, _mm256_div_pd(v_xi2sq, v_seventyTwo)), _mm256_sub_pd(v_one, _mm256_div_pd(v_xi2sq, v_oneHundredTen)))))));
        v_sincxi = _mm256_div_pd(v_sinxi, v_xi);
        v_sincxi2 = _mm256_div_pd(v_sinxi2, v_xi2);
        v_y = _mm256_add_pd(v_sincxi, v_sincxi2);
        _mm256_storeu_pd(y + i, v_y);
    }

    real xi, xi2, xisq, xi2sq;
    for (; i < n; i++)
    {
        xi = x[i];
        xi2 = 2.0 * xi;
        xisq = xi * xi;
        xi2sq = xi2 * xi2;
        y[i] = y[i] = ((xi * (1.0 - ((xisq) / (6.0)) * (1.0 - ((xisq) / (20.0))) * (1.0 - ((xisq) / (42.0)) * (1.0 - ((xisq) / (72.0)) * (1.0 - ((xisq) / (110.0))))))) / (xi)) + ((xi2 * (1.0 - ((xi2sq) / (6.0)) * (1.0 - ((xi2sq) / (20.0))) * (1.0 - ((xi2sq) / (42.0)) * (1.0 - ((xi2sq) / (72.0)) * (1.0 - ((xi2sq) / (110.0))))))) / (xi2));
    }
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
