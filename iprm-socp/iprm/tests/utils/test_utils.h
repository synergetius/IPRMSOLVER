#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include "qoco.h"
#include "gtest/gtest.h"

/**
 * @brief Utility function to extend EXPECT_NEAR to vectors.
 * Passes if ||x-y||_\infty < tol.
 *
 * @param x Input vector.
 * @param y Input vector.
 * @param n Length of vectors.
 * @param tol Tolerace.
 *
 */
void expect_eq_vectorf(IPRMFloat* x, IPRMFloat* y, IPRMInt n, IPRMFloat tol);

/**
 * @brief Utility function to extend EXPECT_EQ to IPRMCscMatrix.
 *
 * @param A Input matrix.
 * @param B Input matrix.
 * @param tol Tolerace.
 */
void expect_eq_csc(IPRMCscMatrix* A, IPRMCscMatrix* B, IPRMFloat tol);

/**
 * @brief Utility function to test if
 * |x-y|/|y| <= tol
 *
 * @param x Tested value.
 * @param y True value.
 * @param tol Tolerance.
 */
void expect_rel_error(IPRMFloat x, IPRMFloat y, IPRMFloat tol);

#endif /* #ifndef TEST_UTILS_H */