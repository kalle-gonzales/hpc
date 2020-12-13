
#include "stopwatch.h"

pstopwatch new_stopwatch() {
  pstopwatch sw;

  sw = (pstopwatch) malloc(sizeof(stopwatch));
#ifndef WIN32
  sw->clk_tck = sysconf(_SC_CLK_TCK);
#endif

  return sw;
}

void del_stopwatch(pstopwatch sw) {
  free(sw);
}

void start_stopwatch(pstopwatch sw) {
#ifdef WIN32
  sw->start = timeGetTime();
#else
  sw->start = omp_get_wtime();
#endif
}

double stop_stopwatch(pstopwatch sw) {
#ifdef WIN32
  sw->current = timeGetTime();
  return (sw->current - sw->start) * 0.001;
#else
  sw->current = omp_get_wtime();
  return (sw->current - sw->start);
#endif
}
