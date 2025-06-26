#include "iprm.h"

int main()
{
  IPRMInt p = 2;     // Number of affine equality constraints (rows of A).
  IPRMInt m = 6;     // Number of conic constraints (rows of G).
  IPRMInt n = 6;     // Number of optimization variables.
  IPRMInt l = 3;     // Dimension of non-negative orthant.
  IPRMInt nsoc = 1;  // Number of second-order cones.
  IPRMInt q[] = {3}; // Dimension of second-order cones.

  IPRMFloat Px[] = {1, 2, 3, 4, 5, 6};  // Data for upper triangular part of P.
  IPRMInt Pnnz = 6;                     // Number of nonzero elements.
  IPRMInt Pp[] = {0, 1, 2, 3, 4, 5, 6}; // Column pointers.
  IPRMInt Pi[] = {0, 1, 2, 3, 4, 5};    // Row indices.

  IPRMFloat Ax[] = {1, 1, 1, 2};
  IPRMInt Annz = 4;
  IPRMInt Ap[] = {0, 1, 3, 4, 4, 4, 4};
  IPRMInt Ai[] = {0, 0, 1, 1};

  IPRMFloat Gx[] = {-1, -1, -1, -1, -1, -1};
  IPRMInt Gnnz = 6;
  IPRMInt Gp[] = {0, 1, 2, 3, 4, 5, 6};
  IPRMInt Gi[] = {0, 1, 2, 3, 4, 5};

  IPRMFloat c[] = {1, 2, 3, 4, 5, 6};
  IPRMFloat b[] = {1, 2};
  IPRMFloat h[] = {0, 0, 0, 0, 0, 0};

  // Allocate storage for data matrices.
  IPRMCscMatrix* P = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* A = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));
  IPRMCscMatrix* G = (IPRMCscMatrix*)malloc(sizeof(IPRMCscMatrix));

  // Set data matrices.
  iprm_set_csc(P, n, n, Pnnz, Px, Pp, Pi);
  iprm_set_csc(A, p, n, Annz, Ax, Ap, Ai);
  iprm_set_csc(G, m, n, Gnnz, Gx, Gp, Gi);

  IPRMSettings* settings = (IPRMSettings*)malloc(sizeof(IPRMSettings));

  // Set default settings.
  iprm_set_default_settings(settings);
  settings->verbose = 1;

  IPRMSolver* solver = (IPRMSolver*)malloc(sizeof(IPRMSolver));
  printf("1\n");
  // Set up solver.
  IPRMInt exit =
      iprm_setup(solver, n, m, p, P, c, A, b, G, h, l, nsoc, q, settings);
  printf("2\n");
  // Solve problem.
  if (exit == IPRM_NO_ERROR) {
    exit = iprm_solve(solver);
  }
  printf("3\n");
  // Free allocated memory.
  iprm_cleanup(solver);
  free(P);
  free(A);
  free(G);
}