#include <stdlib.h>
#include <stdio.h>
#include <immintrin.h>
#include "stopwatch.h"

// Base size of the matrices.
static int SIZE = 1024;
// Each matrix dimension has to be a multiple of this blocksize.
static int BLOCKSIZE = 32;

// Number of iterations to measure the performance.
#define N 3

// Fills a matrix x of size n x m with random numbers within range [-1.0, 1.0]
void random_matrix(float *x, int n, int m)
{
	int i, j;

	srand(42);

	for (j = 0; j < m; ++j)
	{
		for (i = 0; i < n; ++i)
		{
			x[i + j * n] = 2.0f * ((float)rand() / (float)RAND_MAX - 0.5f);
		}
	}
}

// Assume A is n x m matrix, B is m x p, therefore C is n x p matrix.
void matmul_naive(float *A, float *B, float *C, int n, int m, int p)
{
	int i, j, k;

	//TODO: Implement the naive algorithm for matrix-matrix-multiplication
}

// Assume A is n x m matrix, B is m x p, therefore C is n x p matrix.
void matmul_cached(float *A, float *B, float *C, int n, int m, int p)
{

	//TODO: Implement a blocked-algorithm for matrix-matrix-multiplication
}

int main()
{
	float *A, *B, *C, time;
	pstopwatch sw;
	int i;
	uint n, m, p;

	sw = new_stopwatch();

	n = SIZE;
	m = 2 * SIZE;
	p = SIZE;

	printf("################################################################################\n");
	printf("#                        matrix multiplication                                 #\n");
	printf("################################################################################\n\n");

	printf("Matrix A: %d x %d --> %.3f MB\n", n, m, n * m * sizeof(float) / 1024.0 / 1024.0);
	printf("Matrix B: %d x %d --> %.3f MB\n", m, p, m * p * sizeof(float) / 1024.0 / 1024.0);
	printf("Matrix C: %d x %d --> %.3f MB\n\n", n, p, n * p * sizeof(float) / 1024.0 / 1024.0);

	// Allocate storage for the matrices in column-major-order.
	A = (float *)_mm_malloc(n * m * sizeof(float), 64); // n x m
	B = (float *)_mm_malloc(m * p * sizeof(float), 64); // m x p
	C = (float *)_mm_malloc(n * p * sizeof(float), 64); // n x p

	// Init matrices with random values.
	random_matrix(A, n, m);
	random_matrix(B, m, p);
	random_matrix(C, n, p);

	matmul_naive(A, B, C, n, m, p);

	// Test naive matrix-matrix-multiplication.
	start_stopwatch(sw);
	for (i = 0; i < N; ++i)
	{
		matmul_naive(A, B, C, n, m, p);
	}
	time = stop_stopwatch(sw);
	printf("time naive: %.3f s\n", time);

	matmul_cached(A, B, C, n, m, p);

	// Test cache-optimized matrix-matrix-multiplication.
	start_stopwatch(sw);
	for (i = 0; i < N; ++i)
	{
		matmul_cached(A, B, C, n, m, p);
	}
	time = stop_stopwatch(sw);
	printf("time cache-optimized: %.3f s\n", time);

	// Cleanup.
	_mm_free(A);
	_mm_free(B);
	_mm_free(C);

	return 0;
}
