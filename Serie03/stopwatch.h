#ifndef STOPWATCH_H_
#define STOPWATCH_H_

#include <stdlib.h>
#include "omp.h"

/* ------------------------------------------------------------
 Timing
 ------------------------------------------------------------ */

typedef struct _stopwatch stopwatch;
typedef stopwatch *pstopwatch;

struct _stopwatch {
  double start;
  double current;
};

pstopwatch new_stopwatch();

void del_stopwatch(pstopwatch sw);

void start_stopwatch(pstopwatch sw);

double stop_stopwatch(pstopwatch sw);

#endif /* STOPWATCH_H_ */
