
#include "matrix.h"

matrix *
new_matrix(int rows, int cols)
{
    matrix *a;

    a = (matrix *)_mm_malloc(sizeof(matrix), 64);
    a->rows = rows;
    a->cols = cols;
    a->ld = rows;
    a->a = (double *)_mm_malloc(sizeof(double) * a->ld * cols, 64);

    return a;
}

matrix *
new_sub_matrix(matrix *src, int rows, int roff, int cols, int coff)
{
    matrix *a;

    assert(roff + rows <= src->rows);
    assert(coff + cols <= src->cols);

    a = (matrix *)_mm_malloc(sizeof(matrix), 64);
    a->rows = rows;
    a->cols = cols;
    a->ld = src->ld;
    a->a = src->a + roff + coff * src->ld;

    return a;
}

void del_matrix(matrix *a)
{
    _mm_free(a->a);
    a->a = 0;
    _mm_free(a);
}

void copy_matrix(const matrix *a, matrix *b)
{
    int rows = a->rows;
    int cols = a->cols;
    const double *aa = a->a;
    int lda = a->ld;
    double *ba = b->a;
    int ldb = a->ld;
    int i, j;

    assert(b->rows == a->rows);
    assert(b->cols == a->cols);

    for (j = 0; j < cols; j++)
    {
        for (i = 0; i < rows; i++)
        {
            ba[i + j * ldb] = aa[i + j * lda];
        }
    }
}

void random_matrix(matrix *a)
{
    int n = a->rows;
    double *aa = a->a;
    int lda = a->ld;
    int i, j;

    assert(a->cols == a->rows);

    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            aa[i + j * lda] = 2.0 * rand() / RAND_MAX - 1.0;
        }
    }
}

void clear_matrix(matrix *a)
{
    int n = a->rows;
    double *aa = a->a;
    int lda = a->ld;
    int i, j;

    assert(a->cols == a->rows);

    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            aa[i + j * lda] = 0.0;
        }
    }
}

void random_invertible_matrix(matrix *a)
{
    int n = a->rows;
    double *aa = a->a;
    int lda = a->ld;
    double sum;
    int i, j;

    assert(a->cols == a->rows);

    for (j = 0; j < n; j++)
    {
        sum = 0.0;
        for (i = 0; i < n; i++)
        {
            aa[i + j * lda] = 2.0 * rand() / RAND_MAX - 1.0;
            sum += fabs(aa[i + j * lda]);
        }
        aa[j + j * lda] = sum + 1.0;
    }
}

matrix ***build_block_matrix(matrix *a, int m, int p)
{
    int i, j, roff, coff;

    matrix ***as;

    as = (matrix ***)_mm_malloc(sizeof(matrix **) * p, 64);
    roff = 0;
    for (i = 0; i < p; i++)
    {
        as[i] = (matrix **)_mm_malloc(sizeof(matrix *) * p, 64);
        coff = 0;
        for (j = 0; j < p; j++)
        {
            as[i][j] = new_sub_matrix(a, m, roff, m, coff);
            coff += m;
        }
        assert(coff == a->cols);
        roff += m;
    }
    assert(roff == a->rows);

    return as;
}

double diff_matrix(const matrix *a, const matrix *b)
{
    int i, j;
    double error;

    double maxerror;

    assert(a->rows == b->rows);
    assert(a->cols == b->cols);

    maxerror = 0.0;
    for (j = 0; j < a->cols; j++)
    {
        for (i = 0; i < a->rows; i++)
        {
            error = fabs(a->a[i + j * a->ld] - b->a[i + j * b->ld]);
            if (error > maxerror)
                maxerror = error;
        }
    }

    return maxerror;
}

void print_matrix(matrix *a)
{
    int rows = a->rows;
    int cols = a->cols;
    int lda = a->ld;
    int i, j;

    (void)printf("amatrix(%u,%u,%u)\n", rows, cols, a->ld);
    if (rows == 0 || cols == 0)
        return;

    for (i = 0; i < rows; i++)
    {
        (void)printf("  (%.2e", a->a[i]);
        for (j = 1; j < cols; j++)
            (void)printf(" | %.2e", a->a[i + j * lda]);
        (void)printf(")\n");
    }
}
