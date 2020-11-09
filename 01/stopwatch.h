#ifndef STOPWATCH_H_
#define STOPWATCH_H_

#include <stdlib.h>
#include "omp.h"

#ifdef WIN32
#include <Windows.h>
#include <MMSystem.h>
#pragma comment(lib, "winmm")
#else
#include <sys/times.h>
#include <unistd.h>
#endif

/* ------------------------------------------------------------
 Timing
 ------------------------------------------------------------ */

typedef struct _stopwatch stopwatch;
typedef stopwatch *pstopwatch;

struct _stopwatch {
#ifdef WIN32
  DWORD start;
  DWORD current;
#else
  double start;
  double current;
  double clk_tck;
#endif
};

pstopwatch new_stopwatch();

void del_stopwatch(pstopwatch sw);

void start_stopwatch(pstopwatch sw);

double stop_stopwatch(pstopwatch sw);

#endif /* STOPWATCH_H_ */
