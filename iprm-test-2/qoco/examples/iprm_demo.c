#include "qoco.h"

int main()
{
  /*
  QOCOInt p = 0;
  QOCOInt m = 3;
  QOCOInt n = 2;
  QOCOFloat Px[] = {};
  QOCOInt Pnnz = 0;  
  QOCOInt Pp[] = {0};
  QOCOInt Pi[] = {};

  QOCOFloat Ax[] = {}; // A: 0x2
  QOCOInt Annz = 0;
  QOCOInt Ap[] = {0, 0, 0};
  QOCOInt Ai[] = {};

  QOCOFloat Gx[] = {-1, 1, -1, 1};

  QOCOInt Gnnz = 4;
  QOCOInt Gp[] = {0, 2, 4, 4};
  QOCOInt Gi[] = {0, 2, 1, 2};

  QOCOFloat c[] = {-1, -2}; // dim = n
  QOCOFloat b[] = {}; // dim = p
  QOCOFloat h[] = {0, 0, 1}; // dim = m
  */

  
  QOCOInt p = 2; // Number of affine equality constraints (rows of A)
  QOCOInt m = 6; // Number of affine inequality constraints (rows of G)
  QOCOInt n = 6; // Number of optimization variables
  
  QOCOFloat Px[] = {1, 2, 3, 4, 5, 6};  // Data for upper triangular part of P.
  QOCOInt Pnnz = 6;                     // Number of nonzero elements.
  QOCOInt Pp[] = {0, 1, 2, 3, 4, 5, 6}; // Column pointers.
  QOCOInt Pi[] = {0, 1, 2, 3, 4, 5};    // Row indices.
  // P = diag{1, 2, 3, 4, 5, 6}
  QOCOFloat Ax[] = {1, 1, 1, 2};
  QOCOInt Annz = 4;
  QOCOInt Ap[] = {0, 1, 3, 4, 4, 4, 4};
  QOCOInt Ai[] = {0, 0, 1, 1};
  // A: p x n =
  //[[1 1 0 0 0 0]
  // [0 1 2 0 0 0]]
  QOCOFloat Gx[] = {-1, -1, -1, -1, -1, -1};
  QOCOInt Gnnz = 6;
  QOCOInt Gp[] = {0, 1, 2, 3, 4, 5, 6};
  QOCOInt Gi[] = {0, 1, 2, 3, 4, 5};
  // G: m x n = -I_6
  QOCOFloat c[] = {1, 2, 3, 4, 5, 6}; // dim = n
  //obj: 0.5(x1^2+2x2^2+3x3^2+4x4^2+5x5^2+6x6^2) + x1+2x2+3x3+4x4+5x5+6x6
  QOCOFloat b[] = {1, 2}; // dim = p
  QOCOFloat h[] = {0, 0, 0, 0, 0, 0}; // dim = m
  // h = 0, G = -I_6, 这里采用的是x >= 0 的不等式约束条件
  
  QOCOCscMatrix* P = (QOCOCscMatrix*)malloc(sizeof(QOCOCscMatrix));
  QOCOCscMatrix* A = (QOCOCscMatrix*)malloc(sizeof(QOCOCscMatrix));
  QOCOCscMatrix* G = (QOCOCscMatrix*)malloc(sizeof(QOCOCscMatrix));
  
  qoco_set_csc(P, n, n, Pnnz, Px, Pp, Pi);
  qoco_set_csc(A, p, n, Annz, Ax, Ap, Ai);
  qoco_set_csc(G, m, n, Gnnz, Gx, Gp, Gi);
  
  IPRMSettings* settings = (IPRMSettings*)malloc(sizeof(IPRMSettings));
  iprm_set_default_settings(settings);
  settings->verbose = 1;
  IPRMSolver* solver = (IPRMSolver*)malloc(sizeof(IPRMSolver));
  
  QOCOInt exit =
      iprm_setup(solver, n, m, p, P, c, A, b, G, h, settings);
  

  if (exit == QOCO_NO_ERROR) {
    
    exit = iprm_solve(solver);
  }
  /*
  printf("x: ");
  for (int i = 0;i < n;i++){
    printf("%lf ", solver->sol->x);
  }
  printf("\nxi: ");
  for (int i = 0;i < m;i++){
    printf("%lf ", solver->sol->xi);
  }
  printf("\nt: ");
  for (int i = 0;i < p;i++){
    printf("%lf ", solver->sol->t);
  }
  printf("\ns: ");
  for (int i = 0;i < m;i++){
    printf("%lf ", solver->sol->s);
  }
  printf("\n");
  */
  iprm_cleanup(solver);
  free(P);
  free(A);
  free(G);
  return 0;

  
  
}