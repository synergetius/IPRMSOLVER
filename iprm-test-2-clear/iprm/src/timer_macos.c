#include "timer.h"

void start_timer(IPRMTimer* timer) { timer->tic = mach_absolute_time(); }

void stop_timer(IPRMTimer* timer) { timer->toc = mach_absolute_time(); }

IPRMFloat get_elapsed_time_sec(IPRMTimer* timer)
{
  uint64_t duration;

  duration = timer->toc - timer->tic;

  mach_timebase_info(&(timer->tinfo));
  duration *= timer->tinfo.numer;
  duration /= timer->tinfo.denom;

  return (IPRMFloat)duration / 1e9;
}