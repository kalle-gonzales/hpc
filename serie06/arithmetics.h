
#ifndef ARITHMETICS_H
#define ARITHMETICS_H

#include <assert.h>

#include "matrix.h"

void addmul(const double alpha, const matrix *l, const matrix *r, matrix *a);

void lrdecomp(matrix *a);

void lsolve(const matrix *l, matrix *a);

void rsolve(const matrix *r, matrix *a);

void rsolve_trans(const matrix *r, matrix *a);

void block_addmul(const double alpha, matrix ***l, matrix ***r, matrix ***a, int p);

void block_lsolve(matrix ***l, matrix ***a, int p);

void block_rsolve(matrix ***r, matrix ***a, int p);

void block_lrdecomp(matrix ***as, int p);

void block_parallel_addmul(const double alpha, matrix ***l, matrix ***r, matrix ***a, int p);

void block_parallel_lsolve(matrix ***l, matrix ***a, int p);

void block_parallel_rsolve(matrix ***r, matrix ***a, int p);

#endif
