
#include <omp.h>

#include "stopwatch.h"
#include "matrix.h"
#include "arithmetics.h"

int main(int argc, char **argv)
{
    matrix *a, *acopy, *lr, *x, *b;
    matrix ***as, ***bs, ***xs, ***acopys;
    int p, m;
    double maxerror;
    pstopwatch sw;
    double t_lr, t_mul, t_slv, t;

    (void)argc;
    (void)argv;

    p = 32; // Split matrix into p x p many blocks...
    m = 64; // Each block has size m x m.

    sw = new_stopwatch();

    if (argc == 3)
    {
        p = atoi(argv[1]);
        m = atoi(argv[2]);
    }

    printf("Creating matrices of size %d x %d\n", m * p, m * p);
    a = new_matrix(p * m, p * m);
    lr = new_matrix(p * m, p * m);
    acopy = new_matrix(p * m, p * m);
    x = new_matrix(p * m, p * m);
    b = new_matrix(p * m, p * m);

    random_invertible_matrix(a);
    copy_matrix(a, acopy);
    copy_matrix(acopy, lr);
    random_matrix(x);

    printf("Computing reference LR decomposition\n");
    start_stopwatch(sw);
    lrdecomp(lr);
    t_lr = t = stop_stopwatch(sw);
    printf("  %.2f seconds\n", t);

    printf("Splitting into %d x %d submatrices\n", p, p);
    as = build_block_matrix(a, m, p);
    bs = build_block_matrix(b, m, p);
    xs = build_block_matrix(x, m, p);
    acopys = build_block_matrix(acopy, m, p);

    printf("Applying block LR decomposition\n");
    start_stopwatch(sw);
    block_lrdecomp(as, p);
    t = stop_stopwatch(sw);
    printf("  %.2f seconds (x%.2f)\n", t, t_lr / t);

    printf("Comparing results\n");
    maxerror = diff_matrix(a, lr);
    printf("  Maximal error %f\n", maxerror);

    printf("\n===================================\n\n");

    printf("Setup right-hand side (non-blocked)\n");
    start_stopwatch(sw);

    clear_matrix(b);

    addmul(1.0, acopy, x, b);
    t_mul = t = stop_stopwatch(sw);
    printf("  %.2f seconds\n", t);

    printf("Solve linear system (non-blocked)\n");
    start_stopwatch(sw);
    lsolve(a, b);
    rsolve(a, b);
    t_slv = t = stop_stopwatch(sw);
    printf("  %.2f seconds\n", t);

    printf("Comparing solutions\n");
    maxerror = diff_matrix(x, b);
    printf("  Maximal error %f\n", maxerror);

    printf("\n===================================\n\n");

    printf("Setup right-hand side (blocked)\n");
    start_stopwatch(sw);

    clear_matrix(b);

    block_addmul(1.0, acopys, xs, bs, p);
    t = stop_stopwatch(sw);
    printf("  %.2f seconds (x%.2f)\n", t, t_mul / t);

    printf("Solve linear system (blocked)\n");
    start_stopwatch(sw);
    block_lsolve(as, bs, p);
    block_rsolve(as, bs, p);
    t = stop_stopwatch(sw);
    printf("  %.2f seconds (x%.2f)\n", t, t_slv / t);

    printf("Comparing solutions\n");
    maxerror = diff_matrix(x, b);
    printf("  Maximal error %f\n", maxerror);

    printf("\n===================================\n\n");

    printf("Setup right-hand side (in parallel)\n");
    start_stopwatch(sw);

    clear_matrix(b);

    block_parallel_addmul(1.0, acopys, xs, bs, p);
    t = stop_stopwatch(sw);
    printf("  %.2f seconds (x%.2f)\n", t, t_mul / t);

    printf("Solve linear system (in parallel)\n");
    start_stopwatch(sw);
    block_parallel_lsolve(as, bs, p);
    block_parallel_rsolve(as, bs, p);
    t = stop_stopwatch(sw);
    printf("  %.2f seconds (x%.2f)\n", t, t_slv / t);

    printf("Comparing solutions\n");
    maxerror = diff_matrix(x, b);
    printf("  Maximal error %f\n", maxerror);

    del_stopwatch(sw);

    return 0;
}
