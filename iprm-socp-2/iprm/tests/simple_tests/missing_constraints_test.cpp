#include "test_utils.h"
#include "gtest/gtest.h"

#include "qoco.h"

TEST(missing_constraints_test, no_soc_constraints)
{
  QOCOInt p = 2;
  QOCOInt m = 6;
  QOCOInt n = 6;
  QOCOInt l = 6;
  QOCOInt nsoc = 0;

  QOCOFloat Px[] = {1, 2, 3, 4, 5, 6};
  QOCOInt Pnnz = 6;
  QOCOInt Pp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Pi[] = {0, 1, 2, 3, 4, 5};

  QOCOFloat Ax[] = {1, 1, 1, 2};
  QOCOInt Annz = 4;
  QOCOInt Ap[] = {0, 1, 3, 4, 4, 4, 4};
  QOCOInt Ai[] = {0, 0, 1, 1};

  QOCOFloat Gx[] = {-1, -1, -1, -1, -1, -1};
  QOCOInt Gnnz = 6;
  QOCOInt Gp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Gi[] = {0, 1, 2, 3, 4, 5};

  IPRMFloat c[] = {1, 2, 3, 4, 5, 6};
  IPRMFloat b[] = {1, 2};
  IPRMFloat h[] = {0, 0, 0, 0, 0, 0};

  IPRMCscMatrix* P = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* A = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* G = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));

  iprm_set_csc(P, n, n, Pnnz, Px, Pp, Pi);
  iprm_set_csc(A, p, n, Annz, Ax, Ap, Ai);
  iprm_set_csc(G, m, n, Gnnz, Gx, Gp, Gi);

  IPRMFloat xexp[] = {0.2000, 0.8000, 0.6000, -0.0000, 0.0000, 0.0000};
  IPRMFloat sexp[] = {0.2000, 0.8000, 0.6000, 0.0000, 0.0000, 0.0000};
  IPRMFloat yexp[] = {-1.200, -2.400};
  IPRMFloat zexp[] = {0.0000, 0.0000, 0.0000, 4.0000, 5.0000, 6.0000};
  IPRMFloat tol = 1e-4;

  QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
  set_default_settings(settings);
  settings->verbose = 1;

  QOCOSolver* solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

  IPRMInt exit =
      qoco_setup(solver, n, m, p, P, c, A, b, G, h, l, nsoc, nullptr, settings);
  if (exit == IPRM_NO_ERROR) {
    exit = qoco_solve(solver);
  }

  expect_eq_vectorf(solver->sol->x, xexp, n, tol);
  expect_eq_vectorf(solver->sol->s, sexp, m, tol);
  expect_eq_vectorf(solver->sol->y, yexp, p, tol);
  expect_eq_vectorf(solver->sol->z, zexp, n, tol);
  ASSERT_EQ(exit, IPRM_SOLVED);

  qoco_cleanup(solver);
  free(settings);
  free(P);
  free(A);
  free(G);
}

TEST(missing_constraints_test, no_ineq_constraints)
{
  QOCOInt m = 0;
  QOCOInt p = 2;
  QOCOInt n = 6;
  QOCOInt l = 0;
  QOCOInt nsoc = 0;

  QOCOFloat Px[] = {1, 2, 3, 4, 5, 6};
  QOCOInt Pnnz = 6;
  QOCOInt Pp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Pi[] = {0, 1, 2, 3, 4, 5};

  QOCOFloat Ax[] = {1, 1, 1, 2};
  QOCOInt Annz = 4;
  QOCOInt Ap[] = {0, 1, 3, 4, 4, 4, 4};
  QOCOInt Ai[] = {0, 0, 1, 1};

  QOCOFloat c[] = {1, 2, 3, 4, 5, 6};
  QOCOFloat b[] = {1, 2};

  IPRMCscMatrix* P = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* A = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));

  iprm_set_csc(P, n, n, Pnnz, Px, Pp, Pi);
  iprm_set_csc(A, p, n, Annz, Ax, Ap, Ai);

  IPRMFloat xexp[] = {0.2000, 0.8000, 0.6000, -1.0000, -1.0000, -1.0000};
  IPRMFloat yexp[] = {-1.200, -2.400};
  IPRMFloat tol = 1e-4;

  QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
  set_default_settings(settings);
  settings->verbose = 1;

  QOCOSolver* solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

  IPRMInt exit = qoco_setup(solver, n, m, p, P, c, A, b, nullptr, nullptr, l,
                            nsoc, nullptr, settings);
  if (exit == IPRM_NO_ERROR) {
    exit = qoco_solve(solver);
  }

  expect_eq_vectorf(solver->sol->x, xexp, n, tol);
  expect_eq_vectorf(solver->sol->y, yexp, p, tol);
  ASSERT_EQ(exit, IPRM_SOLVED);

  qoco_cleanup(solver);
  free(settings);
  free(P);
  free(A);
}

TEST(missing_constraints_test, no_eq_constraints)
{
  QOCOInt m = 6;
  QOCOInt p = 0;
  QOCOInt n = 6;
  QOCOInt l = 3;
  QOCOInt nsoc = 1;

  QOCOFloat Px[] = {1, 2, 3, 4, 5, 6};
  QOCOInt Pnnz = 6;
  QOCOInt Pp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Pi[] = {0, 1, 2, 3, 4, 5};

  QOCOFloat Gx[] = {-1, -1, -1, -1, -1, -1};
  QOCOInt Gnnz = 6;
  QOCOInt Gp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Gi[] = {0, 1, 2, 3, 4, 5};

  QOCOFloat c[] = {1, 2, 3, 4, 5, 6};
  QOCOFloat h[] = {0, 0, 0, 0, 0, 0};
  QOCOInt q[] = {3};

  IPRMCscMatrix* P = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* G = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));

  iprm_set_csc(P, n, n, Pnnz, Px, Pp, Pi);
  iprm_set_csc(G, m, n, Gnnz, Gx, Gp, Gi);

  IPRMFloat xexp[] = {0.0000, -0.0000, -0.0000, 0.3981, -0.2624, -0.2993};
  IPRMFloat sexp[] = {0.0000, 0.0000, 0.0000, 0.3981, -0.2624, -0.2993};
  IPRMFloat zexp[] = {1.0000, 2.0000, 3.0000, 5.5923, 3.6878, 4.2040};
  IPRMFloat tol = 1e-3;

  QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
  set_default_settings(settings);
  settings->verbose = 1;

  QOCOSolver* solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

  IPRMInt exit = qoco_setup(solver, n, m, p, P, c, nullptr, nullptr, G, h, l,
                            nsoc, q, settings);
  if (exit == IPRM_NO_ERROR) {
    exit = qoco_solve(solver);
  }

  expect_eq_vectorf(solver->sol->x, xexp, n, tol);
  expect_eq_vectorf(solver->sol->s, sexp, m, tol);
  expect_eq_vectorf(solver->sol->z, zexp, n, tol);
  ASSERT_EQ(exit, IPRM_SOLVED);

  qoco_cleanup(solver);
  free(settings);
  free(P);
  free(G);
}

TEST(missing_constraints_test, no_constraints)
{
  QOCOInt m = 0;
  QOCOInt p = 0;
  QOCOInt l = 0;
  QOCOInt n = 6;
  QOCOInt nsoc = 0;

  QOCOFloat Px[] = {1, 2, 3, 4, 5, 6};
  QOCOInt Pnnz = 6;
  QOCOInt Pp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Pi[] = {0, 1, 2, 3, 4, 5};
  QOCOFloat c[] = {1, 2, 3, 4, 5, 6};

  IPRMCscMatrix* P = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  iprm_set_csc(P, n, n, Pnnz, Px, Pp, Pi);

  QOCOFloat xexp[] = {-1, -1, -1, -1, -1, -1};
  QOCOFloat tol = 1e-3;

  QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
  set_default_settings(settings);
  settings->verbose = 1;

  QOCOSolver* solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

  IPRMInt exit = qoco_setup(solver, n, m, p, P, c, nullptr, nullptr, nullptr,
                            nullptr, l, nsoc, nullptr, settings);
  if (exit == IPRM_NO_ERROR) {
    exit = qoco_solve(solver);
  }

  expect_eq_vectorf(solver->sol->x, xexp, n, tol);
  ASSERT_EQ(exit, IPRM_SOLVED);

  qoco_cleanup(solver);
  free(settings);
  free(P);
}

TEST(missing_constraints_test, lp_test_no_P)
{
  //  Solves the following LP
  //  minimize   -x1 - x2
  //  subject to x1 >= 0
  //             x2 >= 0
  //             x1 + x2 <= 1

  QOCOInt p = 0;
  QOCOInt m = 3;
  QOCOInt n = 2;
  QOCOInt l = 3;
  QOCOInt nsoc = 0;

  QOCOFloat Gx[] = {-1, 1, -1, 1};
  QOCOInt Gnnz = 4;
  QOCOInt Gp[] = {0, 2, 4};
  QOCOInt Gi[] = {0, 2, 1, 2};

  QOCOFloat c[] = {-1, -2};
  QOCOFloat h[] = {0, 0, 1};

  IPRMCscMatrix* G = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));

  iprm_set_csc(G, m, n, Gnnz, Gx, Gp, Gi);

  IPRMFloat xexp[] = {0.0, 1.0};
  IPRMFloat sexp[] = {0.0, 1.0, 0.0};
  IPRMFloat zexp[] = {1.0, 0.0, 2.0};
  IPRMFloat tol = 1e-4;

  QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
  set_default_settings(settings);
  settings->verbose = 1;

  QOCOSolver* solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

  IPRMInt exit = qoco_setup(solver, n, m, p, nullptr, c, nullptr, nullptr, G, h,
                            l, nsoc, nullptr, settings);
  if (exit == IPRM_NO_ERROR) {
    exit = qoco_solve(solver);
  }

  expect_eq_vectorf(solver->sol->x, xexp, n, tol);
  expect_eq_vectorf(solver->sol->s, sexp, m, tol);
  expect_eq_vectorf(solver->sol->z, zexp, n, tol);
  ASSERT_EQ(exit, IPRM_SOLVED);

  qoco_cleanup(solver);
  free(settings);
  free(G);
}

TEST(simple_socp_test, p1)
{
  QOCOInt p = 2;
  QOCOInt m = 6;
  QOCOInt n = 6;
  QOCOInt l = 3;
  QOCOInt nsoc = 1;

  QOCOFloat Px[] = {1, 2, 3, 4, 5, 6};
  QOCOInt Pnnz = 6;
  QOCOInt Pp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Pi[] = {0, 1, 2, 3, 4, 5};

  QOCOFloat Ax[] = {1, 1, 1, 2};
  QOCOInt Annz = 4;
  QOCOInt Ap[] = {0, 1, 3, 4, 4, 4, 4};
  QOCOInt Ai[] = {0, 0, 1, 1};

  QOCOFloat Gx[] = {-1, -1, -1, -1, -1, -1};
  QOCOInt Gnnz = 6;
  QOCOInt Gp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Gi[] = {0, 1, 2, 3, 4, 5};

  QOCOFloat c[] = {1, 2, 3, 4, 5, 6};
  IPRMFloat b[] = {1, 2};
  IPRMFloat h[] = {0, 0, 0, 0, 0, 0};
  IPRMInt q[] = {3};

  IPRMCscMatrix* P = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* A = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* G = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));

  iprm_set_csc(P, n, n, Pnnz, Px, Pp, Pi);
  iprm_set_csc(A, p, n, Annz, Ax, Ap, Ai);
  iprm_set_csc(G, m, n, Gnnz, Gx, Gp, Gi);

  IPRMFloat xexp[] = {0.2000, 0.8000, 0.6000, 0.3981, -0.2625, -0.2993};
  IPRMFloat sexp[] = {0.2000, 0.8000, 0.6000, 0.3981, -0.2625, -0.2993};
  IPRMFloat yexp[] = {-1.200, -2.400};
  IPRMFloat zexp[] = {0.0000, 0.0000, 0.0000, 5.5923, 3.6875, 4.2043};
  IPRMFloat tol = 1e-3;

  QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
  set_default_settings(settings);
  settings->verbose = 1;

  QOCOSolver* solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

  IPRMInt exit =
      qoco_setup(solver, n, m, p, P, c, A, b, G, h, l, nsoc, q, settings);
  if (exit == IPRM_NO_ERROR) {
    exit = qoco_solve(solver);
  }

  expect_eq_vectorf(solver->sol->x, xexp, n, tol);
  expect_eq_vectorf(solver->sol->s, sexp, m, tol);
  expect_eq_vectorf(solver->sol->y, yexp, p, tol);
  expect_eq_vectorf(solver->sol->z, zexp, n, tol);
  ASSERT_EQ(exit, IPRM_SOLVED);

  qoco_cleanup(solver);
  free(settings);
  free(P);
  free(A);
  free(G);
}

TEST(simple_socp_test, p2)
{
  QOCOInt p = 2;
  QOCOInt m = 6;
  QOCOInt n = 6;
  QOCOInt l = 1;
  QOCOInt nsoc = 2;

  QOCOFloat Px[] = {1, 2, 3, 4, 5, 6};
  QOCOInt Pnnz = 6;
  QOCOInt Pp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Pi[] = {0, 1, 2, 3, 4, 5};

  QOCOFloat Ax[] = {1, 1, 1, 2};
  QOCOInt Annz = 4;
  QOCOInt Ap[] = {0, 1, 3, 4, 4, 4, 4};
  QOCOInt Ai[] = {0, 0, 1, 1};

  QOCOFloat Gx[] = {-1, -1, -1, -1, -1, -1};
  QOCOInt Gnnz = 6;
  QOCOInt Gp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Gi[] = {0, 1, 2, 3, 4, 5};

  QOCOFloat c[] = {1, 2, 3, 4, 5, 6};
  IPRMFloat b[] = {1, 2};
  IPRMFloat h[] = {0, 0, 0, 0, 0, 0};
  IPRMInt q[] = {2, 3};

  IPRMCscMatrix* P = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* A = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* G = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));

  iprm_set_csc(P, n, n, Pnnz, Px, Pp, Pi);
  iprm_set_csc(A, p, n, Annz, Ax, Ap, Ai);
  iprm_set_csc(G, m, n, Gnnz, Gx, Gp, Gi);

  IPRMFloat xexp[] = {0.2000, 0.8000, 0.6000, 0.3981, -0.2625, -0.2993};
  IPRMFloat sexp[] = {0.2000, 0.8000, 0.6000, 0.3981, -0.2625, -0.2993};
  IPRMFloat yexp[] = {-1.2000, -2.4000};
  IPRMFloat zexp[] = {0.0000, 0.0000, -0.0000, 5.5923, 3.6876, 4.2043};
  IPRMFloat tol = 1e-3;

  QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
  set_default_settings(settings);
  settings->verbose = 1;

  QOCOSolver* solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

  IPRMInt exit =
      qoco_setup(solver, n, m, p, P, c, A, b, G, h, l, nsoc, q, settings);
  if (exit == IPRM_NO_ERROR) {
    exit = qoco_solve(solver);
  }

  expect_eq_vectorf(solver->sol->x, xexp, n, tol);
  expect_eq_vectorf(solver->sol->s, sexp, m, tol);
  expect_eq_vectorf(solver->sol->y, yexp, p, tol);
  expect_eq_vectorf(solver->sol->z, zexp, n, tol);
  ASSERT_EQ(exit, IPRM_SOLVED);

  qoco_cleanup(solver);
  free(settings);
  free(P);
  free(A);
  free(G);
}

TEST(simple_socp_test, p3)
{
  QOCOInt p = 2;
  QOCOInt m = 6;
  QOCOInt n = 6;
  QOCOInt l = 0;
  QOCOInt nsoc = 2;

  QOCOFloat Px[] = {1, 2, 3, 4, 5, 6};
  QOCOInt Pnnz = 6;
  QOCOInt Pp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Pi[] = {0, 1, 2, 3, 4, 5};

  QOCOFloat Ax[] = {1, 1, 1, 2};
  QOCOInt Annz = 4;
  QOCOInt Ap[] = {0, 1, 3, 4, 4, 4, 4};
  QOCOInt Ai[] = {0, 0, 1, 1};

  QOCOFloat Gx[] = {-1, -1, -1, -1, -1, -1};
  QOCOInt Gnnz = 6;
  QOCOInt Gp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Gi[] = {0, 1, 2, 3, 4, 5};

  QOCOFloat c[] = {1, 2, 3, 4, 5, 6};
  IPRMFloat b[] = {1, 2};
  IPRMFloat h[] = {0, 0, 0, 0, 0, 0};
  IPRMInt q[] = {3, 3};

  IPRMCscMatrix* P = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* A = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* G = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));

  iprm_set_csc(P, n, n, Pnnz, Px, Pp, Pi);
  iprm_set_csc(A, p, n, Annz, Ax, Ap, Ai);
  iprm_set_csc(G, m, n, Gnnz, Gx, Gp, Gi);

  IPRMFloat xexp[] = {1.0000, 0.0000, 1.0000, 0.3981, -0.2625, -0.2993};
  IPRMFloat sexp[] = {1.0000, 0.0000, 1.0000, 0.3981, -0.2625, -0.2993};
  IPRMFloat yexp[] = {4.0000, -6.0000};
  IPRMFloat zexp[] = {6.0000, -0.0000, -6.0000, 5.5923, 3.6876, 4.2043};
  IPRMFloat tol = 1e-3;

  QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
  set_default_settings(settings);
  settings->verbose = 1;

  QOCOSolver* solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

  IPRMInt exit =
      qoco_setup(solver, n, m, p, P, c, A, b, G, h, l, nsoc, q, settings);
  if (exit == IPRM_NO_ERROR) {
    exit = qoco_solve(solver);
  }

  expect_eq_vectorf(solver->sol->x, xexp, n, tol);
  expect_eq_vectorf(solver->sol->s, sexp, m, tol);
  expect_eq_vectorf(solver->sol->y, yexp, p, tol);
  expect_eq_vectorf(solver->sol->z, zexp, n, tol);
  ASSERT_EQ(exit, IPRM_SOLVED);

  qoco_cleanup(solver);
  free(settings);
  free(P);
  free(A);
  free(G);
}

TEST(simple_socp_test, TAME)
{
  QOCOInt p = 1;
  QOCOInt m = 2;
  QOCOInt n = 2;
  QOCOInt l = 2;
  QOCOInt nsoc = 0;

  QOCOFloat Px[] = {2, -2, 2};
  QOCOInt Pnnz = 3;
  QOCOInt Pp[] = {0, 1, 3};
  QOCOInt Pi[] = {0, 0, 1};

  QOCOFloat Ax[] = {1, 1};
  QOCOInt Annz = 2;
  QOCOInt Ap[] = {0, 1, 2};
  QOCOInt Ai[] = {0, 0};

  QOCOFloat Gx[] = {-1, -1};
  QOCOInt Gnnz = 2;
  QOCOInt Gp[] = {0, 1, 2};
  QOCOInt Gi[] = {0, 1};

  QOCOFloat c[] = {0, 0};
  IPRMFloat b[] = {1};
  IPRMFloat h[] = {0, 0};

  IPRMCscMatrix* P = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* A = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* G = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));

  iprm_set_csc(P, n, n, Pnnz, Px, Pp, Pi);
  iprm_set_csc(A, p, n, Annz, Ax, Ap, Ai);
  iprm_set_csc(G, m, n, Gnnz, Gx, Gp, Gi);

  QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
  set_default_settings(settings);
  settings->verbose = 1;

  QOCOSolver* solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

  IPRMInt exit =
      qoco_setup(solver, n, m, p, P, c, A, b, G, h, l, nsoc, nullptr, settings);
  if (exit == IPRM_NO_ERROR) {
    exit = qoco_solve(solver);
  }

  QOCOFloat tol = 1e-8;
  IPRMFloat xexp[] = {0.5, 0.5};
  expect_eq_vectorf(solver->sol->x, xexp, n, tol);
  ASSERT_NEAR(solver->sol->obj, 0, tol);
  ASSERT_EQ(exit, IPRM_SOLVED);

  qoco_cleanup(solver);
  free(settings);
  free(P);
  free(A);
  free(G);
}

TEST(simple_socp_test, reduced_tolerance)
{
  QOCOInt p = 2;
  QOCOInt m = 6;
  QOCOInt n = 6;
  QOCOInt l = 3;
  QOCOInt nsoc = 1;

  QOCOFloat Px[] = {1, 2, 3, 4, 5, 6};
  QOCOInt Pnnz = 6;
  QOCOInt Pp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Pi[] = {0, 1, 2, 3, 4, 5};

  QOCOFloat Ax[] = {1, 1, 1, 2};
  QOCOInt Annz = 4;
  QOCOInt Ap[] = {0, 1, 3, 4, 4, 4, 4};
  QOCOInt Ai[] = {0, 0, 1, 1};

  QOCOFloat Gx[] = {-1, -1, -1, -1, -1, -1};
  QOCOInt Gnnz = 6;
  QOCOInt Gp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Gi[] = {0, 1, 2, 3, 4, 5};

  QOCOFloat c[] = {1, 2, 3, 4, 5, 6};
  IPRMFloat b[] = {1, 2};
  IPRMFloat h[] = {0, 0, 0, 0, 0, 0};
  IPRMInt q[] = {3};

  IPRMCscMatrix* P = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* A = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* G = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));

  iprm_set_csc(P, n, n, Pnnz, Px, Pp, Pi);
  iprm_set_csc(A, p, n, Annz, Ax, Ap, Ai);
  iprm_set_csc(G, m, n, Gnnz, Gx, Gp, Gi);

  IPRMFloat xexp[] = {0.2000, 0.8000, 0.6000, 0.3981, -0.2625, -0.2993};
  IPRMFloat sexp[] = {0.2000, 0.8000, 0.6000, 0.3981, -0.2625, -0.2993};
  IPRMFloat yexp[] = {-1.200, -2.400};
  IPRMFloat zexp[] = {0.0000, 0.0000, 0.0000, 5.5923, 3.6875, 4.2043};
  IPRMFloat tol = 1e-4;

  QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
  set_default_settings(settings);
  settings->abstol = 1e-14;
  settings->reltol = 1e-14;
  settings->verbose = 1;

  QOCOSolver* solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

  IPRMInt exit =
      qoco_setup(solver, n, m, p, P, c, A, b, G, h, l, nsoc, q, settings);
  if (exit == IPRM_NO_ERROR) {
    exit = qoco_solve(solver);
  }

  expect_eq_vectorf(solver->sol->x, xexp, n, tol);
  expect_eq_vectorf(solver->sol->s, sexp, m, tol);
  expect_eq_vectorf(solver->sol->y, yexp, p, tol);
  expect_eq_vectorf(solver->sol->z, zexp, n, tol);
  ASSERT_EQ(exit, IPRM_SOLVED_INACCURATE);

  qoco_cleanup(solver);
  free(settings);
  free(P);
  free(A);
  free(G);
}

TEST(simple_socp_test, update_vector_data_test)
{
  QOCOInt p = 2;
  QOCOInt m = 6;
  QOCOInt n = 6;
  QOCOInt l = 3;
  QOCOInt nsoc = 1;

  QOCOFloat Px[] = {1, 2, 3, 4, 5, 6};
  QOCOInt Pnnz = 6;
  QOCOInt Pp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Pi[] = {0, 1, 2, 3, 4, 5};

  QOCOFloat Ax[] = {1, 1, 1, 2};
  QOCOInt Annz = 4;
  QOCOInt Ap[] = {0, 1, 3, 4, 4, 4, 4};
  QOCOInt Ai[] = {0, 0, 1, 1};

  QOCOFloat Gx[] = {-1, -1, -1, -1, -1, -1};
  QOCOInt Gnnz = 6;
  QOCOInt Gp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Gi[] = {0, 1, 2, 3, 4, 5};

  QOCOFloat c[] = {1, 2, 3, 4, 5, 6};
  IPRMFloat b[] = {1, 2};
  IPRMFloat h[] = {0, 0, 0, 0, 0, 0};

  IPRMFloat cnew[] = {0, 0, 0, 0, 0, 0};
  IPRMFloat bnew[] = {4, 5};
  IPRMFloat hnew[] = {1, 1, 1, 1, 1, 1};

  IPRMInt q[] = {3};

  IPRMCscMatrix* P = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* A = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* G = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));

  iprm_set_csc(P, n, n, Pnnz, Px, Pp, Pi);
  iprm_set_csc(A, p, n, Annz, Ax, Ap, Ai);
  iprm_set_csc(G, m, n, Gnnz, Gx, Gp, Gi);

  IPRMFloat xexp[] = {1.9333, 2.0667, 1.4667, 0.2391, -0.1337, -0.1140};
  IPRMFloat sexp[] = {2.9333, 3.0667, 2.4667, 1.2391, 0.8663, 0.8860};
  IPRMFloat yexp[] = {-1.9333, -2.2000};
  IPRMFloat zexp[] = {0.0000, 0.0000, 0.0000, 0.9565, -0.6687, -0.6839};
  IPRMFloat tol = 1e-4;

  QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
  set_default_settings(settings);
  settings->verbose = 1;

  QOCOSolver* solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

  IPRMInt exit =
      qoco_setup(solver, n, m, p, P, c, A, b, G, h, l, nsoc, q, settings);
  if (exit == IPRM_NO_ERROR) {
    exit = qoco_solve(solver);
  }

  update_vector_data(solver, cnew, bnew, hnew);
  exit = qoco_solve(solver);

  expect_eq_vectorf(solver->sol->x, xexp, n, tol);
  expect_eq_vectorf(solver->sol->s, sexp, m, tol);
  expect_eq_vectorf(solver->sol->y, yexp, p, tol);
  expect_eq_vectorf(solver->sol->z, zexp, n, tol);
  ASSERT_EQ(exit, IPRM_SOLVED);

  qoco_cleanup(solver);
  free(settings);
  free(P);
  free(A);
  free(G);
}

TEST(simple_socp_test, update_constraint_data_test)
{
  QOCOInt p = 2;
  QOCOInt m = 6;
  QOCOInt n = 6;
  QOCOInt l = 3;
  QOCOInt nsoc = 1;

  QOCOFloat Px[] = {1, 2, 3, 4, 5, 6};
  QOCOInt Pnnz = 6;
  QOCOInt Pp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Pi[] = {0, 1, 2, 3, 4, 5};

  QOCOFloat Ax[] = {1, 1, 1, 2};
  QOCOFloat Axnew[] = {1, 2, 3, 4};
  QOCOInt Annz = 4;
  QOCOInt Ap[] = {0, 1, 3, 4, 4, 4, 4};
  QOCOInt Ai[] = {0, 0, 1, 1};

  QOCOFloat Gx[] = {-1, -1, -1, -1, -1, -1};
  QOCOInt Gnnz = 6;
  QOCOInt Gp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Gi[] = {0, 1, 2, 3, 4, 5};

  IPRMFloat c[] = {1, 2, 3, 4, 5, 6};
  IPRMFloat b[] = {1, 2};
  IPRMFloat h[] = {0, 0, 0, 0, 0, 0};
  IPRMInt q[] = {3};

  IPRMCscMatrix* P = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* A = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* G = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));

  iprm_set_csc(P, n, n, Pnnz, Px, Pp, Pi);
  iprm_set_csc(A, p, n, Annz, Ax, Ap, Ai);
  iprm_set_csc(G, m, n, Gnnz, Gx, Gp, Gi);

  IPRMFloat xexp[] = {0.0000, 0.5000, 0.1250, 0.3981, -0.2625, -0.2993};
  IPRMFloat sexp[] = {0.0000, 0.5000, 0.1250, 0.3981, -0.2625, -0.2993};
  IPRMFloat yexp[] = {-0.2344, -0.8437};
  IPRMFloat zexp[] = {0.7656, 0.0000, 0.0000, 5.5923, 3.6875, 4.2043};
  IPRMFloat tol = 1e-3;

  QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
  set_default_settings(settings);
  settings->verbose = 1;

  QOCOSolver* solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

  IPRMInt exit =
      qoco_setup(solver, n, m, p, P, c, A, b, G, h, l, nsoc, q, settings);
  if (exit == IPRM_NO_ERROR) {
    exit = qoco_solve(solver);
  }

  update_matrix_data(solver, NULL, Axnew, NULL);

  exit = qoco_solve(solver);

  expect_eq_vectorf(solver->sol->x, xexp, n, tol);
  expect_eq_vectorf(solver->sol->s, sexp, m, tol);
  expect_eq_vectorf(solver->sol->y, yexp, p, tol);
  expect_eq_vectorf(solver->sol->z, zexp, n, tol);
  ASSERT_EQ(exit, IPRM_SOLVED);

  qoco_cleanup(solver);
  free(settings);
  free(P);
  free(A);
  free(G);
}

TEST(simple_socp_test, update_cost_matrix_test)
{
  QOCOInt p = 2;
  QOCOInt m = 6;
  QOCOInt n = 6;
  QOCOInt l = 3;
  QOCOInt nsoc = 1;

  QOCOFloat Px[] = {1, 1, 3, 5, 8};
  QOCOFloat Pxnew[] = {2, 1, 4, 5, 8};

  QOCOInt Pnnz = 5;
  QOCOInt Pp[] = {0, 1, 3, 3, 4, 5, 5};
  QOCOInt Pi[] = {0, 0, 1, 3, 4};

  QOCOFloat Ax[] = {1, 1, 1, 2};
  QOCOFloat Axnew[] = {1, 2, 3, 4};
  QOCOInt Annz = 4;
  QOCOInt Ap[] = {0, 1, 3, 4, 4, 4, 4};
  QOCOInt Ai[] = {0, 0, 1, 1};

  QOCOFloat Gx[] = {-1, -1, -1, -1, -1, -1};
  QOCOInt Gnnz = 6;
  IPRMInt Gp[] = {0, 1, 2, 3, 4, 5, 6};
  IPRMInt Gi[] = {0, 1, 2, 3, 4, 5};

  IPRMFloat c[] = {1, 2, 3, 4, 5, 6};
  IPRMFloat b[] = {1, 2};
  IPRMFloat h[] = {0, 0, 0, 0, 0, 0};
  IPRMInt q[] = {3};

  IPRMCscMatrix* P = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* A = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* G = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));

  iprm_set_csc(P, n, n, Pnnz, Px, Pp, Pi);
  iprm_set_csc(A, p, n, Annz, Ax, Ap, Ai);
  iprm_set_csc(G, m, n, Gnnz, Gx, Gp, Gi);

  IPRMFloat xexp[] = {0.0000, 0.5000, 0.1250, 0.5447, -0.2458, -0.4861};
  IPRMFloat sexp[] = {0.0000, 0.5000, 0.1250, 0.5447, -0.2458, -0.4861};
  IPRMFloat yexp[] = {-0.8750, -0.7500};
  IPRMFloat zexp[] = {0.6250, 0.0000, 0.0000, 6.7234, 3.0338, 6.0000};
  IPRMFloat tol = 1e-4;

  QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
  set_default_settings(settings);
  settings->verbose = 1;

  QOCOSolver* solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

  IPRMInt exit =
      qoco_setup(solver, n, m, p, P, c, A, b, G, h, l, nsoc, q, settings);

  if (exit == IPRM_NO_ERROR) {
    exit = qoco_solve(solver);
  }

  update_matrix_data(solver, Pxnew, Axnew, NULL);

  exit = qoco_solve(solver);

  expect_eq_vectorf(solver->sol->x, xexp, n, tol);
  expect_eq_vectorf(solver->sol->s, sexp, m, tol);
  expect_eq_vectorf(solver->sol->y, yexp, p, tol);
  expect_eq_vectorf(solver->sol->z, zexp, n, tol);
  ASSERT_EQ(exit, IPRM_SOLVED);

  qoco_cleanup(solver);
  free(settings);
  free(P);
  free(A);
  free(G);
}