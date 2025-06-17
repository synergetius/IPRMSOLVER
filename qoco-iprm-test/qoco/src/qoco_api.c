/**
 * @file qoco_api.c
 * @author Govind M. Chari <govindchari1@gmail.com>
 *
 * @section LICENSE
 *
 * Copyright (c) 2024, Govind M. Chari
 * This source code is licensed under the BSD 3-Clause License
 */

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
  //printf("#a.1\n");
  if (iprm_validate_data(P, c, A, b, G, h)) {
    return qoco_error(QOCO_DATA_VALIDATION_ERROR);
  }
  // Validate settings.
  //printf("#a.2\n");
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
  //printf("#a.5\n");

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

  //printf("#a.6\n");
  // Copy P.
  if (P) {
    solver->work->data->P = new_qoco_csc_matrix(P);
  }
  else {
    solver->work->data->P = NULL;
  }
  //printf("#a.7\n");
  // 暂时跳过ruiz equilibration
  solver->work->data->At = create_transposed_matrix(solver->work->data->A);
  //printf("#a.7.1\n");
  solver->work->data->Gt = create_transposed_matrix(solver->work->data->G);
  //printf("#a.7.2\n");
  solver->work->kkt = qoco_malloc(sizeof(IPRMKKT));
  // 对P的静态正则化
  //printf("#a.8\n");
  solver->work->kkt->Pnzadded_idx = qoco_calloc(n, sizeof(QOCOInt));
  //printf("#a.9\n");
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
  //printf("#a.10\n");
  // KKT 的allocate要服务于KKT的构造
  // KKT的构造又要区分哪些是静态的部分、哪些是每一轮中动态更新的部分
  // 上三角中，P和-A^T是静态的，不用随时更新
  // 残差需要动态更新；K的四角（对应上三角部分的三个角块）需要动态更新
  // 左上角块含有静态与动态部分的加法，可以将P保留在内存中，而每一轮动态加到K矩阵中

  iprm_allocate_kkt(solver->work);
  solver->work->kkt->PregtoKKT =
      qoco_calloc(solver->work->data->P->nnz, sizeof(QOCOInt));
  solver->work->kkt->AtoKKT =
      qoco_calloc(solver->work->data->A->nnz, sizeof(QOCOInt));
  solver->work->kkt->GtoKKT =
      qoco_calloc(solver->work->data->G->nnz, sizeof(QOCOInt));
  solver->work->kkt->ZtoKKT = qoco_calloc(m, sizeof(QOCOInt));
  solver->work->kkt->psi = qoco_malloc((n + p + m + m) * sizeof(QOCOFloat));
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
  solver->sol = qoco_malloc(sizeof(QOCOSolution));
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
QOCOInt qoco_setup(QOCOSolver* solver, QOCOInt n, QOCOInt m, QOCOInt p,
                   QOCOCscMatrix* P, QOCOFloat* c, QOCOCscMatrix* A,
                   QOCOFloat* b, QOCOCscMatrix* G, QOCOFloat* h, QOCOInt l,
                   QOCOInt nsoc, QOCOInt* q, QOCOSettings* settings)
{
  // Start setup timer.
  QOCOTimer setup_timer;
  start_timer(&setup_timer);

  // Validate problem data.
  if (qoco_validate_data(P, c, A, b, G, h, l, nsoc, q)) {
    return qoco_error(QOCO_DATA_VALIDATION_ERROR);
  }

  // Validate settings.
  if (qoco_validate_settings(settings)) {
    return qoco_error(QOCO_SETTINGS_VALIDATION_ERROR);
  }

  solver->settings = copy_settings(settings);

  // Allocate workspace.
  solver->work = qoco_malloc(sizeof(QOCOWorkspace));

  // Malloc error.
  if (!(solver->work)) {
    return QOCO_MALLOC_ERROR;
  }

  solver->work->data = qoco_malloc(sizeof(QOCOProblemData));
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
  solver->work->data->q = qoco_malloc(nsoc * sizeof(QOCOInt));

  copy_arrayf(c, solver->work->data->c, n);
  copy_arrayf(b, solver->work->data->b, p);
  copy_arrayf(h, solver->work->data->h, m);
  copy_arrayi(q, solver->work->data->q, nsoc);

  solver->work->data->l = l;
  solver->work->data->nsoc = nsoc;

  // Copy P.
  if (P) {
    solver->work->data->P = new_qoco_csc_matrix(P);
  }
  else {
    solver->work->data->P = NULL;
  }

  // Equilibrate data.
  solver->work->kkt = qoco_malloc(sizeof(QOCOKKT));
  solver->work->kkt->delta = qoco_malloc((n + p + m) * sizeof(QOCOFloat));
  solver->work->kkt->Druiz = qoco_malloc(n * sizeof(QOCOFloat));
  solver->work->kkt->Eruiz = qoco_malloc(p * sizeof(QOCOFloat));
  solver->work->kkt->Fruiz = qoco_malloc(m * sizeof(QOCOFloat));
  solver->work->kkt->Dinvruiz = qoco_malloc(n * sizeof(QOCOFloat));
  solver->work->kkt->Einvruiz = qoco_malloc(p * sizeof(QOCOFloat));
  solver->work->kkt->Finvruiz = qoco_malloc(m * sizeof(QOCOFloat));
  solver->work->data->At = create_transposed_matrix(solver->work->data->A);
  solver->work->data->Gt = create_transposed_matrix(solver->work->data->G);
  ruiz_equilibration(solver);

  // Regularize P.
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

  // Allocate KKT struct.
  allocate_kkt(solver->work);
  solver->work->kkt->nt2kkt = qoco_calloc(solver->work->Wnnz, sizeof(QOCOInt));
  solver->work->kkt->ntdiag2kkt = qoco_calloc(m, sizeof(QOCOInt));
  solver->work->kkt->PregtoKKT =
      qoco_calloc(solver->work->data->P->nnz, sizeof(QOCOInt));
  solver->work->kkt->AtoKKT =
      qoco_calloc(solver->work->data->A->nnz, sizeof(QOCOInt));
  solver->work->kkt->GtoKKT =
      qoco_calloc(solver->work->data->G->nnz, sizeof(QOCOInt));
  solver->work->kkt->rhs = qoco_malloc((n + m + p) * sizeof(QOCOFloat));
  solver->work->kkt->kktres = qoco_malloc((n + m + p) * sizeof(QOCOFloat));
  solver->work->kkt->xyz = qoco_malloc((n + m + p) * sizeof(QOCOFloat));
  solver->work->kkt->xyzbuff1 = qoco_malloc((n + m + p) * sizeof(QOCOFloat));
  solver->work->kkt->xyzbuff2 = qoco_malloc((n + m + p) * sizeof(QOCOFloat));
  construct_kkt(solver);
  //print_K(solver->work->kkt->K);
  // Allocate primal and dual variables.
  solver->work->x = qoco_malloc(n * sizeof(QOCOFloat));
  solver->work->s = qoco_malloc(m * sizeof(QOCOFloat));
  solver->work->y = qoco_malloc(p * sizeof(QOCOFloat));
  solver->work->z = qoco_malloc(m * sizeof(QOCOFloat));
  solver->work->mu = 0.0;

  // Allocate Nesterov-Todd scalings and scaled variables.
  QOCOInt Wnnzfull = solver->work->data->l;
  for (QOCOInt i = 0; i < solver->work->data->nsoc; ++i) {
    Wnnzfull += solver->work->data->q[i] * solver->work->data->q[i];
  }

  solver->work->W = qoco_malloc(solver->work->Wnnz * sizeof(QOCOFloat));
  solver->work->Wfull = qoco_malloc(Wnnzfull * sizeof(QOCOFloat));
  for (int i = 0; i < Wnnzfull; ++i) {
    solver->work->Wfull[i] = 0.0;
  }
  solver->work->Wnnzfull = Wnnzfull;
  solver->work->Winv = qoco_malloc(solver->work->Wnnz * sizeof(QOCOFloat));
  solver->work->Winvfull = qoco_malloc(Wnnzfull * sizeof(QOCOFloat));
  solver->work->WtW = qoco_malloc(solver->work->Wnnz * sizeof(QOCOFloat));
  solver->work->lambda = qoco_malloc(m * sizeof(QOCOFloat));
  QOCOInt qmax = 0;
  if (solver->work->data->nsoc) {
    qmax = max_arrayi(solver->work->data->q, solver->work->data->nsoc);
  }
  solver->work->sbar = qoco_malloc(qmax * sizeof(QOCOFloat));
  solver->work->zbar = qoco_malloc(qmax * sizeof(QOCOFloat));
  solver->work->xbuff = qoco_malloc(n * sizeof(QOCOFloat));
  solver->work->ybuff = qoco_malloc(p * sizeof(QOCOFloat));
  solver->work->ubuff1 = qoco_malloc(m * sizeof(QOCOFloat));
  solver->work->ubuff2 = qoco_malloc(m * sizeof(QOCOFloat));
  solver->work->ubuff3 = qoco_malloc(m * sizeof(QOCOFloat));
  solver->work->Ds = qoco_malloc(m * sizeof(QOCOFloat));

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

  // Compute AMD ordering.
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

  // Update mappings from NT matrix to permuted matrix.
  for (QOCOInt i = 0; i < solver->work->Wnnz; ++i) {
    solver->work->kkt->nt2kkt[i] = KtoPKPt[solver->work->kkt->nt2kkt[i]];
  }
  for (QOCOInt i = 0; i < m; ++i) {
    solver->work->kkt->ntdiag2kkt[i] =
        KtoPKPt[solver->work->kkt->ntdiag2kkt[i]];
  }

  for (QOCOInt i = 0; i < solver->work->data->P->nnz; ++i) {
    solver->work->kkt->PregtoKKT[i] = KtoPKPt[solver->work->kkt->PregtoKKT[i]];
  }

  for (QOCOInt i = 0; i < solver->work->data->A->nnz; ++i) {
    solver->work->kkt->AtoKKT[i] = KtoPKPt[solver->work->kkt->AtoKKT[i]];
  }

  for (QOCOInt i = 0; i < solver->work->data->G->nnz; ++i) {
    solver->work->kkt->GtoKKT[i] = KtoPKPt[solver->work->kkt->GtoKKT[i]];
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
  solver->sol = qoco_malloc(sizeof(QOCOSolution));
  solver->sol->x = qoco_malloc(n * sizeof(QOCOFloat));
  solver->sol->s = qoco_malloc(m * sizeof(QOCOFloat));
  solver->sol->y = qoco_malloc(p * sizeof(QOCOFloat));
  solver->sol->z = qoco_malloc(m * sizeof(QOCOFloat));
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

void set_default_settings(QOCOSettings* settings)
{
  settings->max_iters = 200;
  settings->bisect_iters = 5;
  settings->ruiz_iters = 0;
  settings->iter_ref_iters = 1;
  settings->kkt_static_reg = 1e-8;
  settings->kkt_dynamic_reg = 1e-8;
  settings->abstol = 1e-7;
  settings->reltol = 1e-7;
  settings->abstol_inacc = 1e-5;
  settings->reltol_inacc = 1e-5;
  settings->verbose = 0;
}
void iprm_set_default_settings(IPRMSettings* settings)
{
  settings->max_iters = 10;
  settings->linesearch_iters = 20;
  settings->kkt_static_reg = 1e-8; //
  settings->kkt_dynamic_reg = 1e-8;
  settings->verbose = 0;
  settings->mu0 = 0.1;
  settings->rho0 = 1;
  settings->eta = 10;
  settings->gamma0 = 0.001;
  settings->delta = 0.5;
  settings->tau = 0.01;
  settings->sigma = 0.01;
  settings->epsilon = 1e-5;
  settings->iter_ref_iters = 1;
}

QOCOInt qoco_update_settings(QOCOSolver* solver,
                             const QOCOSettings* new_settings)
{
  if (qoco_validate_settings(new_settings)) {
    return qoco_error(QOCO_SETTINGS_VALIDATION_ERROR);
  }

  solver->settings->max_iters = new_settings->max_iters;
  solver->settings->bisect_iters = new_settings->bisect_iters;
  solver->settings->ruiz_iters = new_settings->ruiz_iters;
  solver->settings->iter_ref_iters = new_settings->iter_ref_iters;
  solver->settings->kkt_static_reg = new_settings->kkt_static_reg;
  solver->settings->kkt_dynamic_reg = new_settings->kkt_dynamic_reg;
  solver->settings->abstol = new_settings->abstol;
  solver->settings->reltol = new_settings->reltol;
  solver->settings->abstol_inacc = new_settings->abstol_inacc;
  solver->settings->abstol_inacc = new_settings->abstol_inacc;
  solver->settings->verbose = new_settings->verbose;

  return 0;
}

void update_vector_data(QOCOSolver* solver, QOCOFloat* cnew, QOCOFloat* bnew,
                        QOCOFloat* hnew)
{
  solver->sol->status = QOCO_UNSOLVED;
  QOCOProblemData* data = solver->work->data;

  // Update cost vector.
  if (cnew) {
    for (QOCOInt i = 0; i < data->n; ++i) {
      data->c[i] = solver->work->kkt->k * solver->work->kkt->Druiz[i] * cnew[i];
    }
  }

  // Update equality constraint vector.
  if (bnew) {
    for (QOCOInt i = 0; i < data->p; ++i) {
      data->b[i] = solver->work->kkt->Eruiz[i] * bnew[i];
    }
  }

  // Update conic constraint vector.
  if (hnew) {
    for (QOCOInt i = 0; i < data->m; ++i) {
      data->h[i] = solver->work->kkt->Fruiz[i] * hnew[i];
    }
  }
}

void update_matrix_data(QOCOSolver* solver, QOCOFloat* Pxnew, QOCOFloat* Axnew,
                        QOCOFloat* Gxnew)
{
  solver->sol->status = QOCO_UNSOLVED;
  QOCOProblemData* data = solver->work->data;
  QOCOKKT* kkt = solver->work->kkt;

  // Undo regularization.
  unregularize(data->P, solver->settings->kkt_static_reg);

  // Unequilibrate P.
  scale_arrayf(data->P->x, data->P->x, kkt->kinv, data->P->nnz);
  row_col_scale(data->P, kkt->Dinvruiz, kkt->Dinvruiz);

  // Unequilibrate c.
  scale_arrayf(data->c, data->c, kkt->kinv, data->n);
  ew_product(data->c, kkt->Dinvruiz, data->c, data->n);

  // Unequilibrate A.
  row_col_scale(data->A, kkt->Einvruiz, kkt->Dinvruiz);

  // Unequilibrate G.
  row_col_scale(data->G, kkt->Finvruiz, kkt->Dinvruiz);

  // Unequilibrate b.
  ew_product(data->b, kkt->Einvruiz, data->b, data->p);

  // Unequilibrate h.
  ew_product(data->h, kkt->Finvruiz, data->h, data->m);

  // Update P and avoid nonzeros that were added for regularization.
  if (Pxnew) {
    QOCOInt avoid =
        kkt->Pnum_nzadded > 0 ? kkt->Pnzadded_idx[0] : data->P->nnz + 1;
    QOCOInt offset = 0;
    for (QOCOInt i = 0; i < data->P->nnz - kkt->Pnum_nzadded; ++i) {
      if (i == avoid) {
        offset++;
        avoid = offset > kkt->Pnum_nzadded ? kkt->Pnzadded_idx[offset]
                                           : data->P->nnz + 1;
      }
      data->P->x[i + offset] = Pxnew[i];
    }
  }

  // Update A.
  if (Axnew) {
    for (QOCOInt i = 0; i < data->A->nnz; ++i) {
      data->A->x[i] = Axnew[i];
    }
  }

  // Update G.
  if (Gxnew) {
    for (QOCOInt i = 0; i < data->G->nnz; ++i) {
      data->G->x[i] = Gxnew[i];
    }
  }

  // Equilibrate new matrix data.
  ruiz_equilibration(solver);

  // Regularize P.
  unregularize(data->P, -solver->settings->kkt_static_reg);

  // Update P in KKT matrix.
  for (QOCOInt i = 0; i < data->P->nnz; ++i) {
    solver->work->kkt->K->x[solver->work->kkt->PregtoKKT[i]] = data->P->x[i];
  }

  // Update A in KKT matrix.
  for (QOCOInt i = 0; i < data->A->nnz; ++i) {
    solver->work->kkt->K->x[solver->work->kkt->AtoKKT[i]] = data->A->x[i];
  }

  // Update G in KKT matrix.
  for (QOCOInt i = 0; i < data->G->nnz; ++i) {
    solver->work->kkt->K->x[solver->work->kkt->GtoKKT[i]] = data->G->x[i];
  }
}
void print_vars(IPRMWorkspace* work){
  /////
  printf("x: ");
  for (int i = 0;i < work->data->n;i++){
    printf("%lf ", work->x[i]);
  }
  printf("\nxi: ");
  for (int i = 0;i < work->data->m;i++){
    printf("%lf ", work->xi[i]);
  }
  printf("\nt: ");
  for (int i = 0;i < work->data->p;i++){
    printf("%lf ", work->t[i]);
  }
  printf("\ns: ");
  for (int i = 0;i < work->data->m;i++){
    printf("%lf ", work->s[i]);
  }
  printf("\n");
  /////
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
  for (QOCOInt i = 1; i <= solver->settings->max_iters; ++i) {
    if (work->mu <= solver->settings->epsilon && 
      work->kkt->phi <= solver->settings->epsilon
    ){
      solver->sol->status = QOCO_SOLVED; // 暂时放在这里。
      // QOCO 的checkstopping中还有更多状态，包括数值错误的异常处理
      stop_timer(&(work->solve_timer));
      iprm_copy_solution(solver);
      if (solver->settings->verbose) {
        iprm_print_footer(solver->sol, solver->sol->status);
      }
      return solver->sol->status;
    }
    work->Dmu = -work->mu + work->gamma * work->kkt->phi;
    iprm_update_blocks(solver);
    iprm_construct_kkt_rhs(work);
    iprm_kkt_solve(solver);
    // Dx, Dt, \hat{Ds} 存储在work->kkt->xts中
    for (QOCOInt j = 0;j < work->data->m;j++){
      work->kkt->xts[work->data->n + work->data->p + j] *= work->mu;
      //把 xts 中的\hat{Ds} 改为 Ds 
    }
    SpMv(work->data->G, work->kkt->xts, work->Dxi);
    for (QOCOInt j = 0;j < work->data->m;j++){
      work->Dxi[j] = - work->Dxi[j] - 
        work->kkt->psi[work->data->n + work->data->p + work->data->m + j];
    }
    iprm_linesearch_update(solver); // 线搜索并更新步长
    work->rho = qoco_max(
      work->rho, 
      solver->settings->sigma * 
      inf_norm(work->s, work->data->m) / 
      qoco_max(1, norm_2(work->x, work->data->n)));
    ////////


    iprm_eval_adjust(solver);

    ////////

    solver->sol->iters = i;
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

QOCOInt qoco_solve(QOCOSolver* solver)
{
  start_timer(&(solver->work->solve_timer));

  // Validate settings.
  if (qoco_validate_settings(solver->settings)) {
    return qoco_error(QOCO_SETTINGS_VALIDATION_ERROR);
  }

  if (solver->settings->verbose) {
    print_header(solver);
  }

  // Get initializations for primal and dual variables.
  initialize_ipm(solver);
  for (QOCOInt i = 1; i <= solver->settings->max_iters; ++i) {

    // Compute kkt residual.
    compute_kkt_residual(solver);

    // Compute mu.
    compute_mu(solver->work);

    // Check stopping criteria.
    if (check_stopping(solver)) {
      stop_timer(&(solver->work->solve_timer));
      unscale_variables(solver->work);
      copy_solution(solver);
      if (solver->settings->verbose) {
        print_footer(solver->sol, solver->sol->status);
      }
      return solver->sol->status;
    }

    // Compute Nesterov-Todd scalings.
    compute_nt_scaling(solver->work);

    // Update Nestrov-Todd block of KKT matrix.
    update_nt_block(solver);

    // Perform predictor-corrector.
    predictor_corrector(solver);

    // Update iteration count.
    solver->sol->iters = i;

    // Log solver progress to console if we are solving in verbose mode.
    if (solver->settings->verbose) {
      log_iter(solver);
    }
  }

  stop_timer(&(solver->work->solve_timer));
  unscale_variables(solver->work);
  copy_solution(solver);
  solver->sol->status = QOCO_MAX_ITER;
  if (solver->settings->verbose) {
    print_footer(solver->sol, solver->sol->status);
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
QOCOInt qoco_cleanup(QOCOSolver* solver)
{

  // Free problem data.
  free_qoco_csc_matrix(solver->work->data->P);
  free_qoco_csc_matrix(solver->work->data->A);
  free_qoco_csc_matrix(solver->work->data->G);
  free_qoco_csc_matrix(solver->work->data->At);
  free_qoco_csc_matrix(solver->work->data->Gt);

  qoco_free(solver->work->data->b);
  qoco_free(solver->work->data->c);
  qoco_free(solver->work->data->h);
  qoco_free(solver->work->data->q);
  qoco_free(solver->work->data);

  // Free primal and dual variables.
  qoco_free(solver->work->kkt->rhs);
  qoco_free(solver->work->kkt->kktres);
  qoco_free(solver->work->kkt->xyz);
  qoco_free(solver->work->kkt->xyzbuff1);
  qoco_free(solver->work->kkt->xyzbuff2);
  qoco_free(solver->work->x);
  qoco_free(solver->work->s);
  qoco_free(solver->work->y);
  qoco_free(solver->work->z);

  // Free Nesterov-Todd scalings and scaled variables.
  qoco_free(solver->work->W);
  qoco_free(solver->work->Wfull);
  qoco_free(solver->work->Winv);
  qoco_free(solver->work->Winvfull);
  qoco_free(solver->work->WtW);
  qoco_free(solver->work->lambda);
  qoco_free(solver->work->sbar);
  qoco_free(solver->work->zbar);
  qoco_free(solver->work->xbuff);
  qoco_free(solver->work->ybuff);
  qoco_free(solver->work->ubuff1);
  qoco_free(solver->work->ubuff2);
  qoco_free(solver->work->ubuff3);
  qoco_free(solver->work->Ds);

  // Free KKT struct.
  free_qoco_csc_matrix(solver->work->kkt->K);
  qoco_free(solver->work->kkt->p);
  qoco_free(solver->work->kkt->pinv);
  qoco_free(solver->work->kkt->delta);
  qoco_free(solver->work->kkt->Druiz);
  qoco_free(solver->work->kkt->Eruiz);
  qoco_free(solver->work->kkt->Fruiz);
  qoco_free(solver->work->kkt->Dinvruiz);
  qoco_free(solver->work->kkt->Einvruiz);
  qoco_free(solver->work->kkt->Finvruiz);
  qoco_free(solver->work->kkt->nt2kkt);
  qoco_free(solver->work->kkt->ntdiag2kkt);
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

  // Free solution struct.
  qoco_free(solver->sol->x);
  qoco_free(solver->sol->s);
  qoco_free(solver->sol->y);
  qoco_free(solver->sol->z);
  qoco_free(solver->sol);

  qoco_free(solver->work);
  qoco_free(solver->settings);
  qoco_free(solver);

  return 1;
}