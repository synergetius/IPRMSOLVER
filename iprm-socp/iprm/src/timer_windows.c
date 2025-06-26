#include "timer.h"

void start_timer(IPRMTimer* timer)
{
  QueryPerformanceFrequency(&timer->freq);
  QueryPerformanceCounter(&timer->tic);
}

void stop_timer(IPRMTimer* timer) { QueryPerformanceCounter(&timer->toc); }

IPRMFloat get_elapsed_time_sec(IPRMTimer* timer)
{
  return (timer->toc.QuadPart - timer->tic.QuadPart) /
         (IPRMFloat)timer->freq.QuadPart;
}
