
#ifndef IS_WINDOWS
#define IS_WINDOWS
#endif

#ifndef DEFINITIONS_H
#define DEFINITIONS_H
// 自己添加操作系统的宏定义

// Define QOCOInt and QOCOFloat.
#include <limits.h>
typedef int IPRMInt;
#define IPRMInt_MAX INT_MAX

#ifdef MATLAB
#define printf mexPrintf
#endif

typedef double IPRMFloat;
#ifdef IS_WINDOWS
#define IPRMFloat_MAX 1e308
#else
#define IPRMFloat_MAX __DBL_MAX__
#endif

#define iprm_max(a, b) (((a) > (b)) ? (a) : (b))
#define iprm_min(a, b) (((a) < (b)) ? (a) : (b))
#define iprm_abs(a) (((a) > 0) ? (a) : (-a))
#define safe_div(a, b) (iprm_abs(b) > 1e-15) ? (a / b) : IPRMFloat_MAX
#include <math.h>
#define iprm_sqrt(a) sqrt(a)

#if defined(IPRM_DEBUG) && !defined(IS_WINDOWS)
#include <assert.h>
#include <stdio.h>
#define iprm_assert(a)                                                         \
  do {                                                                         \
    if (!(a)) {                                                                \
      printf("Assertion Failed: %s\n", #a);                                    \
      __asm__ volatile("int $0x03");                                           \
    }                                                                          \
  } while (0)
#else
#define iprm_assert(a)                                                         \
  do {                                                                         \
  } while (0)
#endif

// Need for malloc, calloc, and free.
#include <stdlib.h>
#define iprm_malloc malloc
#define iprm_calloc calloc
#define iprm_free free
#endif /* ifndef DEFINITIONS_H */