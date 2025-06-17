#include "qoco_api.h"
#include "amd.h"
#include "linalg.h"
#include <stdio.h>
QOCOInt iprm_setup(IPRMSolver* solver, QOCOInt n, QOCOInt m, QOCOInt p,
  QOCOCscMatrix* P, QOCOFloat* c, QOCOCscMatrix* A,
  QOCOFloat* b, QOCOCscMatrix* G, QOCOFloat* h,
  IPRMSettings* settings)
{

  // Start setup timer.
  QOCOTimer setup_timer;
  start_timer(&setup_timer);
  // Validate problem data.
  if (iprm_validate_data(P, c, A, b, G, h)) {
    return qoco_error(QOCO_DATA_VALIDATION_ERROR);
  }
  // Validate settings.
  if (iprm_validate_settings(settings)) {
    return qoco_error(QOCO_SETTINGS_VALIDATION_ERROR);
  }
  solver->settings = iprm_copy_settings(settings);
  // Allocate workspace.
  //printf("#a.3\n");
  solver->work = qoco_malloc(sizeof(IPRMWorkspace));
  // Malloc error.
  if (!(solver->work)) {
    return QOCO_MALLOC_ERROR;
  }
 // printf("#a.4\n");
  solver->work->data = qoco_malloc(sizeof(IPRMProblemData));
  // Malloc error
  if (!(solver->work->data)) {
    return QOCO_MALLOC_ERROR;
  }
  // Copy problem data.
  solver->work->data->m = m;
  solver->work->data->n = n;
  solver->work->data->p = p;
  solver->work->data->A = new_qoco_csc_matrix(A);
  solver->work->data->G = new_qoco_csc_matrix(G);
  solver->work->data->c = qoco_malloc(n * sizeof(QOCOFloat));
  solver->work->data->b = qoco_malloc(p * sizeof(QOCOFloat));
  solver->work->data->h = qoco_malloc(m * sizeof(QOCOFloat));
  copy_arrayf(c, solver->work->data->c, n);
  copy_arrayf(b, solver->work->data->b, p);
  copy_arrayf(h, solver->work->data->h, m);
  // Copy P.
  if (P) {
    solver->work->data->P = new_qoco_csc_matrix(P);
  }
  else {
    solver->work->data->P = NULL;
  }
  solver->work->data->At = create_transposed_matrix(solver->work->data->A);
  solver->work->data->Gt = create_transposed_matrix(solver->work->data->G);
  solver->work->kkt = qoco_malloc(sizeof(IPRMKKT));
  // 对P的静态正则化
  solver->work->kkt->Pnzadded_idx = qoco_calloc(n, sizeof(QOCOInt));
  if (P) {
    solver->work->kkt->Pnum_nzadded =
        regularize(solver->work->data->P, solver->settings->kkt_static_reg,
                   solver->work->kkt->Pnzadded_idx);
  }
  else {
    solver->work->data->P =
        construct_identity(n, solver->settings->kkt_static_reg);
    solver->work->kkt->Pnum_nzadded = n;
  }
  iprm_allocate_kkt(solver->work);
  solver->work->kkt->PregtoKKT =
      qoco_calloc(solver->work->data->P->nnz, sizeof(QOCOInt));
  solver->work->kkt->AtoKKT =
      qoco_calloc(solver->work->data->A->nnz, sizeof(QOCOInt));
  solver->work->kkt->GtoKKT =
      qoco_calloc(solver->work->data->G->nnz, sizeof(QOCOInt));
  solver->work->kkt->ZtoKKT = qoco_calloc(m, sizeof(QOCOInt));
  solver->work->kkt->psi = qoco_malloc((n + p + m + m) * sizeof(QOCOFloat));
  solver->work->kkt->psi_buff = qoco_malloc((n + p + m + m) * sizeof(QOCOFloat));
  solver->work->kkt->psi_buff2 = qoco_malloc((n + p + m + m) * sizeof(QOCOFloat));
  solver->work->kkt->rhs = qoco_malloc((n + p + m) * sizeof(QOCOFloat));
  solver->work->kkt->npm_buff = qoco_malloc((n + m + p) * sizeof(QOCOFloat));
  solver->work->kkt->npm_buff2 = qoco_malloc((n + m + p) * sizeof(QOCOFloat));
  solver->work->kkt->xts = qoco_malloc((n + m + p) * sizeof(QOCOFloat));

  iprm_construct_kkt(solver);
  solver->work->x = qoco_malloc(n * sizeof(QOCOFloat));
  solver->work->t = qoco_malloc(p * sizeof(QOCOFloat));
  solver->work->s = qoco_malloc(m * sizeof(QOCOFloat)); /// !!!
  solver->work->xi = qoco_malloc(m * sizeof(QOCOFloat));
  solver->work->z = qoco_malloc(m * sizeof(QOCOFloat));
  solver->work->y = qoco_malloc(m * sizeof(QOCOFloat));
  solver->work->mu = 0;
  solver->work->rho = 0;
  solver->work->xbuff = qoco_malloc(n * sizeof(QOCOFloat));
  solver->work->tbuff = qoco_malloc(p * sizeof(QOCOFloat));
  solver->work->ubuff1 = qoco_malloc(m * sizeof(QOCOFloat));
  solver->work->Dxi = qoco_malloc(m * sizeof(QOCOFloat));
  // Number of columns of KKT matrix.
  QOCOInt Kn = solver->work->kkt->K->n;
  // Allocate memory for QDLDL.
  solver->work->kkt->etree = qoco_malloc(sizeof(QOCOInt) * Kn);
  solver->work->kkt->Lnz = qoco_malloc(sizeof(QOCOInt) * Kn);
  solver->work->kkt->Lp = qoco_malloc(sizeof(QOCOInt) * (Kn + 1));
  solver->work->kkt->D = qoco_malloc(sizeof(QOCOFloat) * Kn);
  solver->work->kkt->Dinv = qoco_malloc(sizeof(QOCOFloat) * Kn);
  solver->work->kkt->iwork = qoco_malloc(sizeof(QOCOInt) * 3 * Kn);
  solver->work->kkt->bwork = qoco_malloc(sizeof(unsigned char) * Kn);
  solver->work->kkt->fwork = qoco_malloc(sizeof(QOCOFloat) * Kn);

  QOCOCscMatrix* K = solver->work->kkt->K;
  solver->work->kkt->p = qoco_malloc(K->n * sizeof(QOCOInt));
  solver->work->kkt->pinv = qoco_malloc(K->n * sizeof(QOCOInt));
  QOCOInt amd_status = amd_order(K->n, K->p, K->i, solver->work->kkt->p,
                                 (double*)NULL, (double*)NULL);
  if (amd_status < 0) {
    return qoco_error(QOCO_AMD_ERROR);
  }
  invert_permutation(solver->work->kkt->p, solver->work->kkt->pinv, K->n);
  // Permute KKT matrix.
  QOCOInt* KtoPKPt = qoco_malloc(K->nnz * sizeof(QOCOInt));
  QOCOCscMatrix* PKPt = csc_symperm(K, solver->work->kkt->pinv, KtoPKPt);
  // Update mappings
  for (QOCOInt i = 0; i < solver->work->data->P->nnz; ++i) {
    solver->work->kkt->PregtoKKT[i] = KtoPKPt[solver->work->kkt->PregtoKKT[i]];
  }
  for (QOCOInt i = 0; i < solver->work->data->A->nnz; ++i) {
    solver->work->kkt->AtoKKT[i] = KtoPKPt[solver->work->kkt->AtoKKT[i]];
  }
  for (QOCOInt i = 0; i < solver->work->data->G->nnz; ++i) {
    solver->work->kkt->GtoKKT[i] = KtoPKPt[solver->work->kkt->GtoKKT[i]];
  }
  for (QOCOInt i = 0; i < m; ++i) {
    solver->work->kkt->ZtoKKT[i] =
        KtoPKPt[solver->work->kkt->ZtoKKT[i]];
  }
  free_qoco_csc_matrix(solver->work->kkt->K);
  qoco_free(KtoPKPt);
  solver->work->kkt->K = PKPt;
  // Compute elimination tree.
  QOCOInt sumLnz =
      QDLDL_etree(Kn, solver->work->kkt->K->p, solver->work->kkt->K->i,
                  solver->work->kkt->iwork, solver->work->kkt->Lnz,
                  solver->work->kkt->etree);
  if (sumLnz < 0) {
    return QOCO_SETUP_ERROR;
  }
  solver->work->kkt->Li = qoco_malloc(sizeof(QOCOInt) * sumLnz);
  solver->work->kkt->Lx = qoco_malloc(sizeof(QOCOFloat) * sumLnz);
  // Allocate solution struct.
  solver->sol = qoco_malloc(sizeof(IPRMSolution));
  solver->sol->x = qoco_malloc(n * sizeof(QOCOFloat));
  solver->sol->t = qoco_malloc(p * sizeof(QOCOFloat));
  solver->sol->s = qoco_malloc(m * sizeof(QOCOFloat));
  solver->sol->xi = qoco_malloc(m * sizeof(QOCOFloat));
  solver->sol->iters = 0;
  solver->sol->status = QOCO_UNSOLVED;
  stop_timer(&setup_timer);
  solver->sol->setup_time_sec = get_elapsed_time_sec(&setup_timer);

  return QOCO_NO_ERROR;
}


void qoco_set_csc(QOCOCscMatrix* A, QOCOInt m, QOCOInt n, QOCOInt Annz,
                  QOCOFloat* Ax, QOCOInt* Ap, QOCOInt* Ai)
{
  A->m = m;
  A->n = n;
  A->nnz = Annz;
  A->x = Ax;
  A->p = Ap;
  A->i = Ai;
}

void iprm_set_default_settings(IPRMSettings* settings)
{
  settings->max_iters = 200;
  settings->linesearch_iters = 20;
  settings->kkt_static_reg = 1e-8; 
  settings->kkt_dynamic_reg = 1e-8;
  settings->verbose = 0;
  settings->mu0 = 0.1;
  settings->rho0 = 1;
  settings->eta = 10;
  settings->gamma0 = 0.1; ///
  settings->delta = 0.5;
  settings->tau = 0.01;
  settings->sigma = 0.01;
  settings->epsilon = 1e-7;
  settings->iter_ref_iters = 1;
}

QOCOInt iprm_solve(IPRMSolver* solver)
{
  start_timer(&(solver->work->solve_timer));

  IPRMWorkspace* work = solver->work;
  if (iprm_validate_settings(solver->settings)) {
    return qoco_error(QOCO_SETTINGS_VALIDATION_ERROR);
  }
  if (solver->settings->verbose) {
    iprm_print_header(solver);
  }

  iprm_initialize(solver);
  iprm_eval_adjust(solver);
  iprm_check_stopping(solver);
  if (solver->settings->verbose) {
    iprm_log_iter(solver);
  }
  for (QOCOInt i = 1; i <= solver->settings->max_iters; ++i) {
    
    work->Dmu = -work->mu + work->gamma * work->kkt->phi;
    iprm_update_blocks(solver);
    iprm_construct_kkt_rhs(work);
	//printf("||rhs_s||_inf =  %.12e\n", inf_norm(&work->kkt->rhs[work->data->n + work->data->p], work->data->m));
    iprm_kkt_solve(solver);
    // Dx, Dt, Ds 存储在work->kkt->xts中
    SpMv(work->data->G, work->kkt->xts, work->Dxi);
    for (QOCOInt j = 0;j < work->data->m;j++){
      work->Dxi[j] = - work->Dxi[j] - 
        work->kkt->psi[work->data->n + work->data->p + work->data->m + j];
    }
    // printf("||Dxi||_inf = %.12e\n", inf_norm(work->Dxi, work->data->m));
    // printf("||Dx||_inf = %.12e\n", inf_norm(work->kkt->xts, work->data->n));
    // printf("||Dt||_inf = %.12e\n", 
      // inf_norm(&work->kkt->xts[work->data->n], work->data->p));
    // printf("||Ds||_inf = %.12e\n", 
      // inf_norm(&work->kkt->xts[work->data->n + work->data->p], work->data->m));
    // printf("mu = %.12e\n", work->mu);
    // printf("rho = %.12e\n", work->rho);
    // printf("Dmu = %.12e\n", work->Dmu);
    /*
    for (int j = 0;j < work->data->m;j++){
      printf("%.12e ", qoco_sqrt(work->mu) / work->z[j]);
    }
    printf("\n");
    */
    iprm_linesearch_update_fast(solver); // 线搜索并更新步长
    work->rho = qoco_max(
      work->rho, 
      solver->settings->sigma * 
      inf_norm(work->s, work->data->m) /              // 是否要用safe_div？
      qoco_max(1, norm_2(work->x, work->data->n)));
    iprm_eval_adjust(solver);
    solver->sol->iters = i;
	if (iprm_check_stopping(solver)){
      // QOCO 的checkstopping中还有更多状态，包括数值错误的异常处理
      stop_timer(&(work->solve_timer));
      iprm_copy_solution(solver);
      if (solver->settings->verbose) {
        iprm_print_footer(solver->sol, solver->sol->status);
      }
      return solver->sol->status;
    }
    if (solver->settings->verbose) {
      iprm_log_iter(solver);
    }
  }
  stop_timer(&(solver->work->solve_timer));
  iprm_copy_solution(solver);
  solver->sol->status = QOCO_MAX_ITER;
  
  if (solver->settings->verbose) {
    iprm_print_footer(solver->sol, solver->sol->status);
  }
  return QOCO_MAX_ITER;
}
QOCOInt iprm_cleanup(IPRMSolver* solver)
{
  free_qoco_csc_matrix(solver->work->data->P);
  free_qoco_csc_matrix(solver->work->data->A);
  free_qoco_csc_matrix(solver->work->data->G);
  free_qoco_csc_matrix(solver->work->data->At);
  free_qoco_csc_matrix(solver->work->data->Gt);
  qoco_free(solver->work->data->b);
  qoco_free(solver->work->data->c);
  qoco_free(solver->work->data->h);
  qoco_free(solver->work->data);
  qoco_free(solver->work->kkt->rhs);
  qoco_free(solver->work->kkt->psi);
  qoco_free(solver->work->kkt->psi_buff);
  qoco_free(solver->work->kkt->psi_buff2);
  qoco_free(solver->work->kkt->xts);
  qoco_free(solver->work->kkt->npm_buff);
  qoco_free(solver->work->kkt->npm_buff2);
  qoco_free(solver->work->x);
  qoco_free(solver->work->xi);
  qoco_free(solver->work->t);
  qoco_free(solver->work->s);
  qoco_free(solver->work->z);
  qoco_free(solver->work->y);

  qoco_free(solver->work->xbuff);
  qoco_free(solver->work->tbuff);
  qoco_free(solver->work->ubuff1);
  qoco_free(solver->work->Dxi);

  free_qoco_csc_matrix(solver->work->kkt->K);
  qoco_free(solver->work->kkt->p);
  qoco_free(solver->work->kkt->pinv);

  qoco_free(solver->work->kkt->ZtoKKT);
  qoco_free(solver->work->kkt->Pnzadded_idx);
  qoco_free(solver->work->kkt->PregtoKKT);
  qoco_free(solver->work->kkt->AtoKKT);
  qoco_free(solver->work->kkt->GtoKKT);
  qoco_free(solver->work->kkt->etree);
  qoco_free(solver->work->kkt->Lnz);
  qoco_free(solver->work->kkt->Lp);
  qoco_free(solver->work->kkt->D);
  qoco_free(solver->work->kkt->Dinv);
  qoco_free(solver->work->kkt->iwork);
  qoco_free(solver->work->kkt->bwork);
  qoco_free(solver->work->kkt->fwork);
  qoco_free(solver->work->kkt->Li);
  qoco_free(solver->work->kkt->Lx);
  qoco_free(solver->work->kkt);

  qoco_free(solver->sol->x);
  qoco_free(solver->sol->xi);
  qoco_free(solver->sol->t);
  qoco_free(solver->sol->s);
  qoco_free(solver->sol);

  qoco_free(solver->work);
  qoco_free(solver->settings);
  qoco_free(solver);
  return 1;
}
