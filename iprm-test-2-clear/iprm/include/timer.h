
#ifndef TIMER_H
#define TIMER_H

#include "definitions.h"

#ifdef IS_LINUX
#include <time.h>
typedef struct {
  struct timespec tic;
  struct timespec toc;
} IPRMTimer;
#endif

#ifdef IS_MACOS
#include <mach/mach_time.h>
typedef struct {
  uint64_t tic;
  uint64_t toc;
  mach_timebase_info_data_t tinfo;
} IPRMTimer;
#endif

#ifdef IS_WINDOWS
#define NOGDI
#include <windows.h>
typedef struct {
  LARGE_INTEGER tic;
  LARGE_INTEGER toc;
  LARGE_INTEGER freq;
} IPRMTimer;
#endif

/**
 * @brief Starts timer and sets tic field of struct to the current time.
 *
 * @param timer Pointer to timer struct.
 */
void start_timer(IPRMTimer* timer);

/**
 * @brief Stops timer and sets toc field of struct to the current time.
 *
 * @param timer Pointer to timer struct.
 */
void stop_timer(IPRMTimer* timer);

/**
 * @brief Gets time in seconds recorded by timer. Must be called after
 * start_timer() and stop_timer().
 *
 * @param timer Pointer to timer struct.
 */
IPRMFloat get_elapsed_time_sec(IPRMTimer* timer);

#endif /* #ifndef TIMER_H */