
#include "utils.h"

void print_iprm_csc_matrix(IPRMCscMatrix* M)
{
  printf("\nPrinting CSC Matrix:\n");
  printf("m: %d\n", M->m);
  printf("n: %d\n", M->n);
  printf("nnz: %d\n", M->nnz);
  printf("Data: {");
  for (IPRMInt i = 0; i < M->nnz; ++i) {
    printf("%.17g", M->x[i]);
    if (i != M->nnz - 1) {
      printf(",");
    }
  }
  printf("}\n");

  printf("Row Indices: {");
  for (IPRMInt i = 0; i < M->nnz; ++i) {
    printf("%d", M->i[i]);
    if (i != M->nnz - 1) {
      printf(",");
    }
  }
  printf("}\n");

  printf("Column Pointers: {");
  for (IPRMInt i = 0; i < M->n + 1; ++i) {
    printf("%d", M->p[i]);
    if (i != M->n) {
      printf(",");
    }
  }
  printf("}\n");
}

void print_arrayf(IPRMFloat* x, IPRMInt n)
{
  printf("{");
  for (IPRMInt i = 0; i < n; ++i) {
    printf("%.17g", x[i]);
    if (i != n - 1) {
      printf(", ");
    }
  }
  printf("}\n");
}

void print_arrayi(IPRMInt* x, IPRMInt n)
{
  printf("{");
  for (IPRMInt i = 0; i < n; ++i) {
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
  printf("+--------+-------------+---------------+-------------+----------+---------+\n");
  printf("|  Iter  |     Obj     |      Phi      |     Pres    |   Dres   |   Gap   |\n");
  printf("+--------+-------------+---------------+-------------+----------+---------+\n");
}
void iprm_log_iter(IPRMSolver* solver)
{
  printf("|   %2d   | %+.5e | %+.5e | %+.5e | %+.5e | %+.5e |\n",
         solver->sol->iters, solver->sol->obj, solver->work->kkt->phi, solver->sol->pres, solver->sol->dres, solver->sol->gap);
  printf("+--------+-------------+---------------+-------------+----------+---------+\n");
}
void iprm_print_footer(IPRMSolution* solution, enum iprm_solve_status status)
{
  printf("\n");
  //printf("status:                %s\n", IPRM_SOLVE_STATUS_MESSAGE[status]);
  printf("number of iterations:  %d\n", solution->iters);
  printf("phi:             %+.3f\n", solution->phi);
  printf("setup time:            %.2e sec\n", solution->setup_time_sec);
  printf("solve time:            %.2e sec\n", solution->solve_time_sec);
  printf("\n");
}
unsigned char iprm_check_stopping(IPRMSolver* solver)
{
	
	IPRMWorkspace* work = solver->work;
	IPRMProblemData* data = solver->work->data;
	
	IPRMFloat eabs = solver->settings->epsilon;
	IPRMFloat erel = solver->settings->epsilon;
	IPRMFloat binf = data->p > 0 ? inf_norm(data->b, data->p) : 0;
	IPRMFloat xiinf = data->m > 0 ? inf_norm(work->xi, data->m) : 0;
	IPRMFloat sinf = data->m > 0 ? inf_norm(work->s, data->m) : 0;
	IPRMFloat cinf = inf_norm(data->c, data->n);
	IPRMFloat hinf = data->m > 0 ? inf_norm(data->h, data->m) : 0;
	
	SpMtv(data->A, work->t, work->xbuff);
	IPRMFloat Attinf = data->p ? inf_norm(work->xbuff, data->n) : 0;
	
	SpMtv(data->G, work->s, work->xbuff);
	IPRMFloat Gtsinf = data->m > 0 ? inf_norm(work->xbuff, data->n) : 0;
	
	USpMv(data->P, work->x, work->xbuff);
	for (IPRMInt i = 0; i < data->n; ++i) {
		work->xbuff[i] -= solver->settings->kkt_static_reg * work->x[i];
	}
	IPRMFloat Pxinf = inf_norm(work->xbuff, data->n);
	IPRMFloat xPx = dot(work->x, work->xbuff, work->data->n);
	
	SpMv(data->A, work->x, work->tbuff);
	IPRMFloat Axinf = data->p ? inf_norm(work->tbuff, data->p) : 0;

	SpMv(data->G, work->x, work->sbuff);
  IPRMFloat Gxinf = data->m ? inf_norm(work->sbuff, data->m) : 0;
  // Compute primal residual.
  IPRMFloat eq_res = inf_norm(&work->kkt->psi[data->n], data->p);
  IPRMFloat ineq_res = inf_norm(&work->kkt->psi[data->n + data->p + data->m], data->m);
  IPRMFloat pres = iprm_max(eq_res, ineq_res);
  solver->sol->pres = pres;
  // Compute dual residual.
  IPRMFloat dres = inf_norm(work->kkt->psi, data->n);
  solver->sol->dres = dres;
  // Compute gap: || z - xi ||_\infty
  IPRMFloat gap = norm_2(&work->kkt->psi[data->n + data->p], data->m); //////// ??
  //IPRMFloat gap = inf_norm(&work->kkt->psi[data->n + data->p], data->m);
  solver->sol->gap = gap;
  // Compute max{Axinf, binf, Gxinf, hinf, xiinf}.
  IPRMFloat pres_rel = iprm_max(Axinf, binf);
  pres_rel = iprm_max(pres_rel, Gxinf);
  pres_rel = iprm_max(pres_rel, hinf);
  pres_rel = iprm_max(pres_rel, xiinf);
  // Compute max{Pxinf, Attinf, Gtsinf, cinf}.
  IPRMFloat dres_rel = iprm_max(Pxinf, Attinf);
  dres_rel = iprm_max(dres_rel, Gtsinf);
  dres_rel = iprm_max(dres_rel, cinf);
  // Compute max{1, abs(pobj), abs(dobj)}.
  IPRMFloat ctx = dot(work->data->c, work->x, work->data->n);
  IPRMFloat btt = dot(work->data->b, work->t, work->data->p);
  IPRMFloat hts = dot(work->data->h, work->s, work->data->m);
  IPRMFloat pobj = 0.5 * xPx + ctx;
  IPRMFloat dobj = -0.5 * xPx - btt - hts;
  pobj = iprm_abs(pobj);
  dobj = iprm_abs(dobj);
  IPRMFloat gap_rel = iprm_max(1, pobj);
  gap_rel = iprm_max(gap_rel, dobj);
  
  
  //printf("pres: %.12e, eabs + erel * pres_rel:%.12e\n", pres, eabs + erel * pres_rel);
  //printf("dres: %.12e, eabs + erel * dres_rel:%.12e\n", dres, eabs + erel * dres_rel);
  //printf("solver->sol->gap:%.12e, eabs + erel * gap_rel:%.12e\n", solver->sol->gap, eabs + erel * gap_rel);
  
  //printf("target pres:%+.5e dres:%+.5e gap:%+.5e\n", eabs + erel * pres_rel, eabs + erel * dres_rel, eabs + erel * gap_rel);
  
  // IPRM论文中的终止条件
  /*
  if (solver->work->mu < solver->settings->epsilon && work->kkt->phi < solver->settings->epsilon){
    solver->sol->status = IPRM_SOLVED;
    return 1;
  }
  */
 
  // IPRSOCP论文中的终止条件
  if (isnan(work->kkt->phi) || isinf(work->kkt->phi)){
    solver->sol->status = IPRM_NUMERICAL_ERROR;
    return 1;
  }
  if (solver->work->mu < solver->settings->epsilon && work->kkt->square_psi < solver->settings->epsilon){
    solver->sol->status = IPRM_SOLVED;
    return 1;
  }
 
  // QOCO论文中的终止条件
  // if (pres < eabs + erel * pres_rel && dres < eabs + erel * dres_rel &&
      // solver->sol->gap < eabs + erel * gap_rel) {
    // solver->sol->status = IPRM_SOLVED;
    // printf("%+.5e %+.5e\n", solver->sol->gap, eabs + erel * gap_rel);
    // printf("check solved\n");
    // return 1;
  // }
  return 0;
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