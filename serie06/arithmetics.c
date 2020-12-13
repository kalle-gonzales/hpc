
#include "arithmetics.h"

void addmul(const double alpha, const matrix *l, const matrix *r, matrix *a)
{
    int rows = a->rows;
    int cols = a->cols;
    const double *la = l->a;
    int ldl = l->ld;
    const double *ra = r->a;
    int ldr = r->ld;
    double *aa = a->a;
    int lda = a->ld;
    int i, j, k;

    assert(l->rows == a->rows);
    assert(l->cols == r->rows);
    assert(r->cols == a->cols);

    for (k = 0; k < cols; k++)
    {
        for (j = 0; j < l->cols; j++)
        {
            for (i = 0; i < rows; i++)
            {
                aa[i + k * lda] += alpha * la[i + j * ldl] * ra[j + k * ldr];
            }
        }
    }
}

void lrdecomp(matrix *a)
{
    int n = a->rows;
    double *aa = a->a;
    int lda = a->ld;
    double alpha;
    int i, j, k;

    assert(a->cols == a->rows);

    for (k = 0; k < n; k++)
    {
        alpha = 1.0 / aa[k + k * lda];
        for (i = k + 1; i < n; i++)
        {
            aa[i + k * lda] *= alpha;
        }

        for (j = k + 1; j < n; j++)
        {
            for (i = k + 1; i < n; i++)
            {
                aa[i + j * lda] -= aa[i + k * lda] * aa[k + j * lda];
            }
        }
    }
}

void lsolve(const matrix *l, matrix *a)
{
    int n = l->rows;
    const double *la = l->a;
    int ldl = l->ld;
    double *aa = a->a;
    int lda = a->ld;
    int i, j, k;

    assert(l->cols == l->rows);
    assert(a->rows == l->cols);

    for (j = 0; j < a->cols; j++)
    {
        for (k = 0; k < n; k++)
        {
            for (i = k + 1; i < n; i++)
            {
                aa[i + j * lda] -= la[i + k * ldl] * aa[k + j * lda];
            }
        }
    }
}

void rsolve_trans(const matrix *r, matrix *a)
{
    int n = r->rows;
    const double *ra = r->a;
    int ldr = r->ld;
    double *aa = a->a;
    int lda = a->ld;
    double alpha;
    int i, j, k;

    assert(r->cols == r->rows);
    assert(a->cols == r->rows);

    for (k = 0; k < n; k++)
    {
        alpha = 1.0 / ra[k + k * ldr];
        for (i = 0; i < a->rows; i++)
        {
            aa[i + k * lda] *= alpha;
        }

        for (j = k + 1; j < n; j++)
        {
            for (i = 0; i < a->rows; i++)
            {
                aa[i + j * lda] -= aa[i + k * lda] * ra[k + j * ldr];
            }
        }
    }
}

void rsolve(const matrix *r, matrix *a)
{
    int n = r->rows;
    const double *ra = r->a;
    int ldr = r->ld;
    double *aa = a->a;
    int lda = a->ld;
    int i, j, k;

    assert(r->cols == r->rows);
    assert(a->cols == r->rows);

    for (j = 0; j < a->cols; j++)
    {
        for (k = n; k-- > 0;)
        {
            aa[k + j * lda] /= ra[k + k * ldr];
            for (i = 0; i < k; i++)
            {
                aa[i + j * lda] -= ra[i + k * ldr] * aa[k + j * lda];
            }
        }
    }
}

/*
 * Matrices l, r and a are divided into p x p subblocks.
 */
void block_addmul(const double alpha, matrix ***l, matrix ***r, matrix ***a, int p)
{
    int i, j, k;

    for (k = 0; k < p; k++)
    {
        for (j = 0; j < p; j++)
        {
            for (i = 0; i < p; i++)
            {
                addmul(alpha, l[i][j], r[j][k], a[i][k]);
            }
        }
    }
}

/*
 * Matrices l and a are divided into p x p subblocks.
 */
void block_lsolve(matrix ***l, matrix ***a, int p)
{
    int i, j, k;

    for(j = 0; j < p; j++) {
        for(i = 0; i < p; i++) {
            
            lsolve(l[i][i], a[i][j]);

            for(k = i + 1; k < p; k++) {
                addmul(-1.0,l[k][i], a[i][j], a[k][j]);
            }
        }
    }
}

/*
 * Matrices r and a are divided into p x p subblocks.
 */
void block_rsolve(matrix ***r, matrix ***a, int p)
{
    int i, j, k;

    for(j = 0; j < p; j++) {
        for(i = p-1; i >= 0; i--) {
            
            rsolve(r[i][i], a[i][j]);

            for(k = 0; k < i; k++) {
                addmul(-1.0,r[k][i], a[i][j], a[k][j]);
            }
        }
    }
}

/*
 * Matrix a is divided into p x p subblocks.
 */
void block_lrdecomp(matrix ***as, int p)
{
    int k, j, i;

    for (k = 0; k < p; k++)
    {
        lrdecomp(as[k][k]);

        for (j = k + 1; j < p; j++)
        {
            lsolve(as[k][k], as[k][j]);
        }

        for (i = k + 1; i < p; i++)
        {
            rsolve_trans(as[k][k], as[i][k]);
        }

        for (i = k + 1; i < p; i++)
        {
            for (j = k + 1; j < p; j++)
            {
                addmul(-1.0, as[i][k], as[k][j], as[i][j]);
            }
        }
    }
}

/*
 * Matrices l, r and a are divided into p x p subblocks.
 */
void block_parallel_addmul(const double alpha, matrix ***l, matrix ***r, matrix ***a, int p)
{
    int i, j, k;

    #pragma omp parallel 
    {
        #pragma omp single
        for(k = 0; k < p; k++) {
            for(j = 0; j < p; j++) {
                for(i = 0; i < p; i++) {
                #pragma omp task firstprivate(i,j,k) depend(in: l[i][j], r[j][k]) depend(inout: a[i][k])
                addmul(alpha, l[i][j], r[j][k], a[i][k]);
                }
            }
        }
    }
}

/*
 * Matrices l and a are divided into p x p subblocks.
 */
void block_parallel_lsolve(matrix ***l, matrix ***a, int p)
{
    int i, j, k;

    #pragma omp parallel
    {
        #pragma omp single
        for(j = 0; j < p; j++) {
            for(i = 0; i < p; i++) {
                #pragma omp task firstprivate(i,j) depend(in: l[i][i]) depend(inout: a[i][j])
                lsolve(l[i][i], a[i][j]);
                
                for(k = i + 1; k < p; k++) {
                    #pragma omp task firstprivate(i,j,k), depend(in: l[k][i], a[i][j]) depend(inout: a[k][j])
                    addmul(-1.0, l[k][i], a[i][j], a[k][j]);
                }
            }
        }
    }
}

/*
 * Matrices r and a are divided into p x p subblocks.
 */
void block_parallel_rsolve(matrix ***r, matrix ***a, int p)
{
    int i, j, k;

#pragma omp parallel
#pragma omp single
    for(j = 0; j < p; j++) {
        for(i = p - 1; i >= 0; i--) {
#pragma omp task firstprivate(i,j) depend(in: r[i][i]) depend(inout: a[i][j])
            rsolve(r[i][i], a[i][j]);
            
            for(k = 0; k < i; k++) {
#pragma omp task firstprivate(i,j,k), depend(in: r[k][i], a[i][j]) depend(inout: a[k][j])
                addmul(-1.0, r[k][i], a[i][j], a[k][j]);
            }
        }
    }
}
