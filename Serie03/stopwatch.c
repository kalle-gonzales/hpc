
#include "stopwatch.h"

pstopwatch new_stopwatch() {
  pstopwatch sw;

  sw = (pstopwatch) malloc(sizeof(stopwatch));
  
  return sw;
}

void del_stopwatch(pstopwatch sw) {
  free(sw);
}

void start_stopwatch(pstopwatch sw) {
  sw->start = omp_get_wtime();
}

double stop_stopwatch(pstopwatch sw) {
  sw->current = omp_get_wtime();
  return (sw->current - sw->start);
}
