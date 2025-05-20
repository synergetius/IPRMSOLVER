
#include "utils.h"

void print_qoco_csc_matrix(QOCOCscMatrix* M)
{
  printf("\nPrinting CSC Matrix:\n");
  printf("m: %d\n", M->m);
  printf("n: %d\n", M->n);
  printf("nnz: %d\n", M->nnz);
  printf("Data: {");
  for (QOCOInt i = 0; i < M->nnz; ++i) {
    printf("%.17g", M->x[i]);
    if (i != M->nnz - 1) {
      printf(",");
    }
  }
  printf("}\n");

  printf("Row Indices: {");
  for (QOCOInt i = 0; i < M->nnz; ++i) {
    printf("%d", M->i[i]);
    if (i != M->nnz - 1) {
      printf(",");
    }
  }
  printf("}\n");

  printf("Column Pointers: {");
  for (QOCOInt i = 0; i < M->n + 1; ++i) {
    printf("%d", M->p[i]);
    if (i != M->n) {
      printf(",");
    }
  }
  printf("}\n");
}

void print_arrayf(QOCOFloat* x, QOCOInt n)
{
  printf("{");
  for (QOCOInt i = 0; i < n; ++i) {
    printf("%.17g", x[i]);
    if (i != n - 1) {
      printf(", ");
    }
  }
  printf("}\n");
}

void print_arrayi(QOCOInt* x, QOCOInt n)
{
  printf("{");
  for (QOCOInt i = 0; i < n; ++i) {
    printf("%d", x[i]);
    if (i != n - 1) {
      printf(", ");
    }
  }
  printf("}\n");
}
void iprm_print_header(IPRMSolver* solver)
{
  IPRMProblemData* data = solver->work->data;
  IPRMSettings* settings = solver->settings;
  printf("\n");
  printf("+-------------------------------------------------------+\n");
  printf("| Problem Data:                                         |\n");
  printf("|     variables:        %-9d                       |\n", data->n);
  printf("|     constraints:      %-9d                       |\n", data->m + data->p);
  printf("|     eq constraints:   %-9d                       |\n", data->p);
  printf("|     ineq constraints: %-9d                       |\n", data->m);
  printf("|     nnz(P):           %-9d                       |\n", data->P->nnz - solver->work->kkt->Pnum_nzadded);
  printf("|     nnz(A):           %-9d                       |\n", data->A->nnz);
  printf("|     nnz(G):           %-9d                       |\n", data->G->nnz);
  printf("| Solver Settings:                                      |\n");
  printf("|     max_iters:        %-9d                       |\n", settings->max_iters);
  printf("|     kkt_static_reg:   %3.2e                      |\n", settings->kkt_static_reg);
  printf("|     kkt_dynamic_reg:  %3.2e                      |\n", settings->kkt_dynamic_reg);
  printf("+-------------------------------------------------------+\n");
  printf("\n");
  printf("+--------+------------+\n");
  printf("|  Iter  |   Phi      |\n");
  printf("+--------+------------+\n");
}
void iprm_log_iter(IPRMSolver* solver)
{
  printf("|   %2d   | %+.2e |\n",
         solver->sol->iters, solver->work->kkt->phi);
  printf("+--------+------------+\n");
}
void iprm_print_footer(IPRMSolution* solution, enum qoco_solve_status status)
{
  printf("\n");
  //printf("status:                %s\n", QOCO_SOLVE_STATUS_MESSAGE[status]);
  printf("number of iterations:  %d\n", solution->iters);
  printf("phi:             %+.3f\n", solution->phi);
  printf("setup time:            %.2e sec\n", solution->setup_time_sec);
  printf("solve time:            %.2e sec\n", solution->solve_time_sec);
  printf("\n");
}
void iprm_copy_solution(IPRMSolver* solver)
{
  copy_arrayf(solver->work->x, solver->sol->x, solver->work->data->n);
  copy_arrayf(solver->work->t, solver->sol->t, solver->work->data->p);
  copy_arrayf(solver->work->s, solver->sol->s, solver->work->data->m);
  copy_arrayf(solver->work->xi, solver->sol->xi, solver->work->data->m);
  solver->sol->phi = solver->work->kkt->phi;
  solver->sol->solve_time_sec =
      get_elapsed_time_sec(&(solver->work->solve_timer));
}
IPRMSettings* iprm_copy_settings(IPRMSettings* settings)
{
  IPRMSettings* new_settings = malloc(sizeof(IPRMSettings));
  new_settings->max_iters = settings->max_iters;
  new_settings->linesearch_iters = settings->linesearch_iters;
  new_settings->iter_ref_iters = settings->iter_ref_iters;
  new_settings->kkt_static_reg = settings->kkt_static_reg;
  new_settings->kkt_dynamic_reg = settings->kkt_dynamic_reg;
  new_settings->verbose = settings->verbose;
  new_settings->mu0 = settings->mu0;
  new_settings->rho0 = settings->rho0;
  new_settings->eta = settings->eta;
  new_settings->gamma0 = settings->gamma0;
  new_settings->delta = settings->delta;
  new_settings->tau = settings->tau;
  new_settings->sigma = settings->sigma;
  new_settings->epsilon = settings->epsilon;

  return new_settings;
}