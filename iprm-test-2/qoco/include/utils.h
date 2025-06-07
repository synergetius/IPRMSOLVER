
#ifndef UTILS_H
#define UTILS_H

#include "enums.h"
#include "linalg.h"
#include "structs.h"
#include <stdio.h>

/**
 * @brief Prints dimensions, number of nonzero elements, data, column pointers
 * and row indices for a sparse matrix in CSC form.
 *
 * @param M Pointer to QOCOCscMatrix that will be printed.
 */
void print_qoco_csc_matrix(QOCOCscMatrix* M);
/**
 * @brief Prints array of QOCOFloats.
 *
 * @param x Pointer to array.
 * @param n Number of elements in array.
 */
void print_arrayf(QOCOFloat* x, QOCOInt n);
/**
 * @brief Prints array of QOCOInts.
 *
 * @param x Pointer to array.
 * @param n Number of elements in array.
 */
void print_arrayi(QOCOInt* x, QOCOInt n);
void iprm_print_header(IPRMSolver* solver);
void iprm_log_iter(IPRMSolver* solver);
void iprm_print_footer(IPRMSolution* solution, enum qoco_solve_status status);
unsigned char iprm_check_stopping(IPRMSolver* solver);
void iprm_copy_solution(IPRMSolver* solver);
IPRMSettings* iprm_copy_settings(IPRMSettings* settings);

#endif /* #ifndef UTILS_H */