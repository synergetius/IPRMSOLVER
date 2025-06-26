#include "test_utils.h"

void expect_eq_vectorf(IPRMFloat* x, IPRMFloat* y, IPRMInt n, IPRMFloat tol)
{
  for (IPRMInt i = 0; i < n; ++i) {
    EXPECT_NEAR(x[i], y[i], tol);
  }
}

void expect_eq_csc(IPRMCscMatrix* A, IPRMCscMatrix* B, IPRMFloat tol)
{
  EXPECT_EQ(A->m, B->m);
  EXPECT_EQ(A->n, B->n);
  EXPECT_EQ(A->nnz, B->nnz);

  for (IPRMInt i = 0; i < A->nnz; ++i) {
    EXPECT_EQ(A->i[i], B->i[i]);
  }

  for (IPRMInt i = 0; i < A->n + 1; ++i) {
    EXPECT_EQ(A->p[i], B->p[i]);
  }

  expect_eq_vectorf(A->x, B->x, A->nnz, tol);
}

void expect_rel_error(IPRMFloat x, IPRMFloat y, IPRMFloat tol)
{
  IPRMFloat err = x - y;
  err = iprm_abs(err);
  IPRMFloat yabs = iprm_abs(y);
  err = err / yabs;
  EXPECT_LE(err, tol);
}