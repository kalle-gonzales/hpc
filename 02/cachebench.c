#include <stdio.h>
#include <immintrin.h>
#include "stopwatch.h"

inline void vector_add(const float alpha, const float *x, float *y, uint n) {
  size_t i;

  /* TODO: Impelement the addition of the two vectors x and y, y <-- y + alpha * x */
  
}

int main(int argc, char **argv) {
  size_t n, iter;
  float *x, *y;
  size_t i;
  pstopwatch sw;
  double t, data, flops;
  FILE *file;

  sw = new_stopwatch();

  if (argc == 2) {
    n = atoi(argv[1]);
  } else {
    n = 1 << 20;
  }

  iter = (1ul << 31) / n;

  printf("Testing with vectors of length %lu and %lu iterations\n", n, iter);
  printf("Vectors consume %.2f KB of memory.\n",
      2 * n * sizeof(float) / 1024.0);

  x = (float*) _mm_malloc(n * sizeof(float), 64);
  y = (float*) _mm_malloc(n * sizeof(float), 64);

  srand(42);

  for (i = 0; i < n; ++i) {
    x[i] = 2.0f + (float) rand() / (float) RAND_MAX;
    y[i] = 1.0f + (float) rand() / (float) RAND_MAX;
  }

  /* Cache warm-up */
  vector_add((float) iter / 2, x, y, n);
  vector_add((float) iter / 2, x, y, n);

  /* actual test with filled caches */
  start_stopwatch(sw);
  for (i = 0; i < iter; ++i) {
    vector_add(-1.0f, x, y, n);
  }
  t = stop_stopwatch(sw);
  
  /* Print out the values stdout */
  printf("  %.3f ms\n", t * 1.0e3);
  data = 3 * n * sizeof(float) * iter / (t * 1024.0 * 1024.0 * 1024.0);
  printf("  %.3f GB/s\n", data);
  flops = 2 * n * iter / (t * 1024.0 * 1024.0 * 1024.0);
  printf("  %.3f GFlops\n", flops);

  /* Print out the values to a file 'data.dat' */
  file = fopen("data.dat", "a+");
  fprintf(file, "%lu\t%.3e\n", n, data);
  fclose(file);

  del_stopwatch(sw);
  _mm_free(x);
  _mm_free(y);

  return 0;
}
