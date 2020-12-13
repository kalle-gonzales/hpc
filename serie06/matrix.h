

#ifndef MATRIX_H
#define MATRIX_H

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

typedef struct
{
    double *a; // Coefficients for the matrix, stored column-major-order
    int rows;
    int cols;
    int ld; /* Leading dimension. Is equal to rows for a general matrix,
   * for a submatrix it is equal to 'rows' of the source matrix,
   * because it is needed to index the coefficients of the matrix. */
} matrix;

matrix *new_matrix(int rows, int cols);

matrix *new_sub_matrix(matrix *src, int rows, int roff, int cols, int coff);

void del_matrix(matrix *a);

void copy_matrix(const matrix *a, matrix *b);

void random_matrix(matrix *a);

void clear_matrix(matrix *a);

void random_invertible_matrix(matrix *a);

matrix ***build_block_matrix(matrix *a, int m, int p);

double diff_matrix(const matrix *a, const matrix *b);

void print_matrix(matrix *a);

#endif
