#ifndef ENUMS_H
#define ENUMS_H

/**
 * @brief Enum for error codes.
 *
 */
enum iprm_error_code {
  IPRM_NO_ERROR = 0,

  // Error in problem data validation.
  IPRM_DATA_VALIDATION_ERROR,

  // Error in settings validation.
  IPRM_SETTINGS_VALIDATION_ERROR,

  // Error in setup.
  IPRM_SETUP_ERROR,

  // Error in performing amd ordering.
  IPRM_AMD_ERROR,

  // Memory allocation error.
  IPRM_MALLOC_ERROR
};

/**
 * @brief Enum for solver status.
 *
 */
enum iprm_solve_status {
  // Unsolved (Solver needs to be called.)
  IPRM_UNSOLVED = 0,

  // Solved successfully.
  IPRM_SOLVED = 1,

  // Solved Inaccurately.
  IPRM_SOLVED_INACCURATE,

  // Numerical error (occurs when a = 0 and inaccurate stopping criteria not
  // met).
  IPRM_NUMERICAL_ERROR,

  // Maximum number of iterations reached.
  IPRM_MAX_ITER,
};

// clang-format off
static const char *IPRM_ERROR_MESSAGE[] = {
    "", // Error codes start from 1.
    "data validation error",
    "settings validation error",
    "amd error",
    "memory allocation error"
};


static const char *IPRM_SOLVE_STATUS_MESSAGE[] = {
    "unsolved", // Solver not run.
    "solved",
    "solved inaccurately",
    "numerical error",
    "maximum iterations reached",
};
// clang-format on

#endif /* #ifndef ENUMS_H */