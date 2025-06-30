#include "iprm_api.h"
#include "amd.h"
#include "linalg.h"
#include <stdio.h>
IPRMInt iprm_setup(IPRMSolver* solver, IPRMInt n, IPRMInt m, IPRMInt p,
  IPRMCscMatrix* P, IPRMFloat* c, IPRMCscMatrix* A,
  IPRMFloat* b, IPRMCscMatrix* G, IPRMFloat* h, IPRMInt l,
  IPRMInt nsoc, IPRMInt* q, IPRMSettings* settings)
{

  // Start setup timer.
  IPRMTimer setup_timer;
  start_timer(&setup_timer);
  // Validate problem data.
  if (iprm_validate_data(P, c, A, b, G, h, l, nsoc, q)) {
    return iprm_error(IPRM_DATA_VALIDATION_ERROR);
  }
  // Validate settings.
  if (iprm_validate_settings(settings)) {
    return iprm_error(IPRM_SETTINGS_VALIDATION_ERROR);
  }
  solver->settings = iprm_copy_settings(settings);
  // Allocate workspace.
  //printf("#a.3\n");
  solver->work = iprm_malloc(sizeof(IPRMWorkspace));
  // Malloc error.
  if (!(solver->work)) {
    return IPRM_MALLOC_ERROR;
  }
 // printf("#a.4\n");
  solver->work->data = iprm_malloc(sizeof(IPRMProblemData));
  // Malloc error
  if (!(solver->work->data)) {
    return IPRM_MALLOC_ERROR;
  }
  // Copy problem data.
  solver->work->data->m = m;
  solver->work->data->n = n;
  solver->work->data->p = p;
  solver->work->data->A = new_iprm_csc_matrix(A);
  solver->work->data->G = new_iprm_csc_matrix(G);
  solver->work->data->c = iprm_malloc(n * sizeof(IPRMFloat));
  solver->work->data->b = iprm_malloc(p * sizeof(IPRMFloat));
  solver->work->data->h = iprm_malloc(m * sizeof(IPRMFloat));
  solver->work->data->q = iprm_malloc(nsoc * sizeof(IPRMInt));
  copy_arrayf(c, solver->work->data->c, n);
  copy_arrayf(b, solver->work->data->b, p);
  copy_arrayf(h, solver->work->data->h, m);
  copy_arrayi(q, solver->work->data->q, nsoc);
  solver->work->data->l = l;
  solver->work->data->nsoc = nsoc;
  // Copy P.
  if (P) {
    solver->work->data->P = new_iprm_csc_matrix(P);
  }
  else {
    solver->work->data->P = NULL;
  }
  solver->work->data->At = create_transposed_matrix(solver->work->data->A);
  solver->work->data->Gt = create_transposed_matrix(solver->work->data->G);
  solver->work->kkt = iprm_malloc(sizeof(IPRMKKT));
  // 对P的静态正则化
  solver->work->kkt->Pnzadded_idx = iprm_calloc(n, sizeof(IPRMInt));
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
  //printf("1.1\n");
  iprm_allocate_kkt(solver->work);
  ///printf("1.2\n");
  solver->work->kkt->PregtoKKT =
      iprm_calloc(solver->work->data->P->nnz, sizeof(IPRMInt));
  solver->work->kkt->AtoKKT =
      iprm_calloc(solver->work->data->A->nnz, sizeof(IPRMInt));
  solver->work->kkt->GtoKKT =
      iprm_calloc(solver->work->data->G->nnz, sizeof(IPRMInt));
  solver->work->kkt->B55toKKT = iprm_calloc(l, sizeof(IPRMInt));
  solver->work->kkt->B26toKKT = iprm_calloc(m, sizeof(IPRMInt));
  solver->work->kkt->B36toKKT = iprm_calloc(nsoc, sizeof(IPRMInt));
  solver->work->kkt->B66toKKT = iprm_calloc(m, sizeof(IPRMInt));
  //solver->work->kkt->ZtoKKT = iprm_calloc(m, sizeof(IPRMInt));
  solver->work->kkt->psi = iprm_malloc((n + p + m + m) * sizeof(IPRMFloat));
  solver->work->kkt->psi_buff = iprm_malloc((n + p + m + m) * sizeof(IPRMFloat));
  solver->work->kkt->psi_buff2 = iprm_malloc((n + p + m + m) * sizeof(IPRMFloat));
  solver->work->kkt->rhs = iprm_malloc(solver->work->kkt->K->n * sizeof(IPRMFloat)); // 
  solver->work->kkt->kktbuff1 = iprm_malloc(solver->work->kkt->K->n * sizeof(IPRMFloat));
  solver->work->kkt->kktbuff2 = iprm_malloc(solver->work->kkt->K->n * sizeof(IPRMFloat));
  solver->work->kkt->kktbuff3 = iprm_malloc(solver->work->kkt->K->n * sizeof(IPRMFloat));
  solver->work->kkt->w = iprm_malloc(solver->work->kkt->K->n * sizeof(IPRMFloat));
  //solver->work->kkt->npm_buff = iprm_malloc((n + m + p) * sizeof(IPRMFloat));
  //solver->work->kkt->npm_buff2 = iprm_malloc((n + m + p) * sizeof(IPRMFloat));
  //solver->work->kkt->xts = iprm_malloc((n + m + p) * sizeof(IPRMFloat));
  //printf("1.3\n");
  iprm_construct_kkt(solver);
  //printf("1.4\n");
  solver->work->x = iprm_malloc(n * sizeof(IPRMFloat));
  solver->work->t = iprm_malloc(p * sizeof(IPRMFloat));
  solver->work->s = iprm_malloc(m * sizeof(IPRMFloat)); /// !!!
  solver->work->xi = iprm_malloc(m * sizeof(IPRMFloat));
  solver->work->z = iprm_malloc(m * sizeof(IPRMFloat));
  solver->work->y = iprm_malloc(m * sizeof(IPRMFloat));
  solver->work->u = iprm_malloc(m * sizeof(IPRMFloat));
  solver->work->v = iprm_malloc(m * sizeof(IPRMFloat));
  solver->work->mu = 0;
  solver->work->rho = 0;
  solver->work->xbuff = iprm_malloc(n * sizeof(IPRMFloat));
  solver->work->tbuff = iprm_malloc(p * sizeof(IPRMFloat));
  solver->work->sbuff = iprm_malloc(m * sizeof(IPRMFloat));
  solver->work->Dxi = iprm_malloc(m * sizeof(IPRMFloat));
  // Number of columns of KKT matrix.
  //printf("1.5\n");
  IPRMInt Kn = solver->work->kkt->K->n;
  // Allocate memory for QDLDL.
  solver->work->kkt->etree = iprm_malloc(sizeof(IPRMInt) * Kn);
  solver->work->kkt->Lnz = iprm_malloc(sizeof(IPRMInt) * Kn);
  solver->work->kkt->Lp = iprm_malloc(sizeof(IPRMInt) * (Kn + 1));
  solver->work->kkt->D = iprm_malloc(sizeof(IPRMFloat) * Kn);
  solver->work->kkt->Dinv = iprm_malloc(sizeof(IPRMFloat) * Kn);
  solver->work->kkt->iwork = iprm_malloc(sizeof(IPRMInt) * 3 * Kn);
  solver->work->kkt->bwork = iprm_malloc(sizeof(unsigned char) * Kn);
  solver->work->kkt->fwork = iprm_malloc(sizeof(IPRMFloat) * Kn);

  IPRMCscMatrix* K = solver->work->kkt->K;
  ///////////////////////////////////////////
  
  solver->work->kkt->p = iprm_malloc(K->n * sizeof(IPRMInt));
  solver->work->kkt->pinv = iprm_malloc(K->n * sizeof(IPRMInt));
  
  IPRMInt amd_status = amd_order(K->n, K->p, K->i, solver->work->kkt->p,
                                 (double*)NULL, (double*)NULL);
  // p没有问题
  
 // for (int i = 0;i < K->n;i++){
  //  printf("%d ", solver->work->kkt->p[i]);
  //}
  //printf("\n");
  
  ///////////////// 发现问题：排序之后的K丢掉了有根号2数值的元素
  if (amd_status < 0) {
    return iprm_error(IPRM_AMD_ERROR);
  }
  //printf("1.6\n");
  invert_permutation(solver->work->kkt->p, solver->work->kkt->pinv, K->n);
  // Permute KKT matrix.
  IPRMInt* KtoPKPt = iprm_malloc(K->nnz * sizeof(IPRMInt));
  IPRMCscMatrix* PKPt = csc_symperm(K, solver->work->kkt->pinv, KtoPKPt);
  //printf("KtoPKPt:\n");
  //for (int i = 0;i < K->nnz;i++){
  //  printf("%d ", KtoPKPt[i]);
  //}
  //printf("\n");
  // Update mappings
  //printf("1.7\n");
  for (IPRMInt i = 0; i < solver->work->data->P->nnz; ++i) {
    solver->work->kkt->PregtoKKT[i] = KtoPKPt[solver->work->kkt->PregtoKKT[i]];
  }
  for (IPRMInt i = 0; i < solver->work->data->A->nnz; ++i) {
    solver->work->kkt->AtoKKT[i] = KtoPKPt[solver->work->kkt->AtoKKT[i]];
  }
  for (IPRMInt i = 0; i < solver->work->data->G->nnz; ++i) {
    solver->work->kkt->GtoKKT[i] = KtoPKPt[solver->work->kkt->GtoKKT[i]];
  }
  for (IPRMInt i = 0; i < l; ++i) {
    solver->work->kkt->B55toKKT[i] =
        KtoPKPt[solver->work->kkt->B55toKKT[i]];
  }
  //printf("1.8\n");
  for (IPRMInt i = l; i < m; ++i) {
    solver->work->kkt->B26toKKT[i] =
        KtoPKPt[solver->work->kkt->B26toKKT[i]];
  }
  for (IPRMInt i = l; i < m; ++i) {
    solver->work->kkt->B66toKKT[i] =
        KtoPKPt[solver->work->kkt->B66toKKT[i]];
  }
  for (IPRMInt i = 0; i < nsoc; ++i) {
    solver->work->kkt->B36toKKT[i] =
        KtoPKPt[solver->work->kkt->B36toKKT[i]];
  }
 // printf("1.9\n");
  //print_K(solver->work->kkt->K);
  //printf("Knnz = %d\n", solver->work->kkt->K->nnz);
  free_iprm_csc_matrix(solver->work->kkt->K);
  iprm_free(KtoPKPt);
  solver->work->kkt->K = PKPt;
  
  //////////////////////////////////////////
  //print_K(solver->work->kkt->K);
  //print_iprm_csc_matrix(solver->work->kkt->K);
  
  // Compute elimination tree.
  IPRMInt sumLnz =
      QDLDL_etree(Kn, solver->work->kkt->K->p, solver->work->kkt->K->i,
                  solver->work->kkt->iwork, solver->work->kkt->Lnz,
                  solver->work->kkt->etree);
  if (sumLnz < 0) {
    return IPRM_SETUP_ERROR;
  }
  solver->work->kkt->Li = iprm_malloc(sizeof(IPRMInt) * sumLnz);
  solver->work->kkt->Lx = iprm_malloc(sizeof(IPRMFloat) * sumLnz);
  // Allocate solution struct.
  solver->sol = iprm_malloc(sizeof(IPRMSolution));
  solver->sol->x = iprm_malloc(n * sizeof(IPRMFloat));
  solver->sol->t = iprm_malloc(p * sizeof(IPRMFloat));
  solver->sol->s = iprm_malloc(m * sizeof(IPRMFloat));
  solver->sol->xi = iprm_malloc(m * sizeof(IPRMFloat));
  solver->sol->iters = 0;
  solver->sol->status = IPRM_UNSOLVED;
  stop_timer(&setup_timer);
  solver->sol->setup_time_sec = get_elapsed_time_sec(&setup_timer);

  return IPRM_NO_ERROR;
}


void iprm_set_csc(IPRMCscMatrix* A, IPRMInt m, IPRMInt n, IPRMInt Annz,
                  IPRMFloat* Ax, IPRMInt* Ap, IPRMInt* Ai)
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
  settings->max_iters = 20;
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
  settings->iter_ref_iters = 1;//10;
}

IPRMInt iprm_solve(IPRMSolver* solver)
{
  start_timer(&(solver->work->solve_timer));

  IPRMWorkspace* work = solver->work;
  if (iprm_validate_settings(solver->settings)) {
    return iprm_error(IPRM_SETTINGS_VALIDATION_ERROR);
  }
  if (solver->settings->verbose) {
    iprm_print_header(solver);
  }
  //printf("2.1\n");
  
  iprm_initialize(solver);
  //printf("finished init\n");
  //print_K(work->kkt->K);
  //printf("2.2\n");
  
  iprm_eval_adjust(solver);
  //printf("2.3\n");
  iprm_check_stopping(solver);
  if (solver->settings->verbose) {
    iprm_log_iter(solver);
  }
  for (IPRMInt i = 1; i <= solver->settings->max_iters; ++i) {
    
    work->Dmu = -work->mu + work->gamma * work->kkt->phi;
    iprm_update_blocks(solver);
    
    iprm_construct_kkt_rhs(work);
    

	  // printf("||rhs_s||_inf =  %.12e\n", inf_norm(&work->kkt->rhs[work->data->n + work->data->nsoc + work->data->nsoc + work->data->p], work->data->m));
    iprm_kkt_solve(solver);
    
    // Dx, Dt, Ds 存储在work->kkt->w中
    SpMv(work->data->G, work->kkt->w, work->Dxi);
    for (IPRMInt j = 0;j < work->data->m;j++){
      work->Dxi[j] = - work->Dxi[j] - 
        work->kkt->psi[work->data->n + work->data->p + work->data->m + j];
    }
    
    // printf("||x||_inf = %.12e\n", inf_norm(work->x, work->data->n));
    // printf("||t||_inf = %.12e\n", inf_norm(work->t, work->data->p));
    // printf("||s||_inf = %.12e\n", inf_norm(work->s, work->data->m));
    // printf("||xi||_inf = %.12e\n", inf_norm(work->xi, work->data->m));
    // printf("||z||_inf = %.12e\n", inf_norm(work->z, work->data->m));
    // printf("||Dxi||_inf = %.12e\n", inf_norm(work->Dxi, work->data->m));
    // printf("||Dx||_inf = %.12e\n", inf_norm(work->kkt->w, work->data->n));
    // printf("||Dt||_inf = %.12e\n", 
      // inf_norm(&work->kkt->w[work->data->n + work->data->nsoc + work->data->nsoc], work->data->p));
    // printf("||Ds||_inf = %.12e\n", 
      // inf_norm(&work->kkt->w[work->data->n + work->data->nsoc + work->data->nsoc + work->data->p], work->data->m));
    // printf("mu = %.12e\n", work->mu);
    // printf("rho = %.12e\n", work->rho);
    // printf("Dmu = %.12e\n", work->Dmu);
    /*
    for (int j = 0;j < work->data->m;j++){
      printf("%.12e ", iprm_sqrt(work->mu) / work->z[j]);
    }
    printf("\n");
    */
    iprm_linesearch_update(solver); // 线搜索并更新步长
    
    work->rho = iprm_max(
      work->rho, 
      solver->settings->sigma * 
      inf_norm(work->s, work->data->m) /              // 是否要用safe_div？
      iprm_max(1, norm_2(work->x, work->data->n)));
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
  solver->sol->status = IPRM_MAX_ITER;
  
  if (solver->settings->verbose) {
    iprm_print_footer(solver->sol, solver->sol->status);
  }
  return IPRM_MAX_ITER;
}
IPRMInt iprm_cleanup(IPRMSolver* solver)
{
  free_iprm_csc_matrix(solver->work->data->P);
  free_iprm_csc_matrix(solver->work->data->A);
  free_iprm_csc_matrix(solver->work->data->G);
  free_iprm_csc_matrix(solver->work->data->At);
  free_iprm_csc_matrix(solver->work->data->Gt);
  iprm_free(solver->work->data->b);
  iprm_free(solver->work->data->c);
  iprm_free(solver->work->data->h);
  iprm_free(solver->work->data);
  iprm_free(solver->work->kkt->rhs);
  iprm_free(solver->work->kkt->psi);
  iprm_free(solver->work->kkt->psi_buff);
  iprm_free(solver->work->kkt->psi_buff2);
  //iprm_free(solver->work->kkt->xts);
  iprm_free(solver->work->kkt->w);
  iprm_free(solver->work->kkt->kktbuff1);
  iprm_free(solver->work->kkt->kktbuff2);
  iprm_free(solver->work->kkt->kktbuff3);
  //iprm_free(solver->work->kkt->npm_buff);
  //iprm_free(solver->work->kkt->npm_buff2);
  iprm_free(solver->work->x);
  iprm_free(solver->work->xi);
  iprm_free(solver->work->t);
  iprm_free(solver->work->s);
  iprm_free(solver->work->z);
  iprm_free(solver->work->y);
  iprm_free(solver->work->u);
  iprm_free(solver->work->v);

  iprm_free(solver->work->xbuff);
  iprm_free(solver->work->tbuff);
  iprm_free(solver->work->sbuff);
  iprm_free(solver->work->Dxi);

  free_iprm_csc_matrix(solver->work->kkt->K);
  iprm_free(solver->work->kkt->p);
  iprm_free(solver->work->kkt->pinv);

  //iprm_free(solver->work->kkt->ZtoKKT);
  iprm_free(solver->work->kkt->Pnzadded_idx);
  iprm_free(solver->work->kkt->PregtoKKT);
  iprm_free(solver->work->kkt->AtoKKT);
  iprm_free(solver->work->kkt->GtoKKT);
  iprm_free(solver->work->kkt->B55toKKT);
  iprm_free(solver->work->kkt->B26toKKT);
  iprm_free(solver->work->kkt->B36toKKT);
  iprm_free(solver->work->kkt->B66toKKT);
  iprm_free(solver->work->kkt->etree);
  iprm_free(solver->work->kkt->Lnz);
  iprm_free(solver->work->kkt->Lp);
  iprm_free(solver->work->kkt->D);
  iprm_free(solver->work->kkt->Dinv);
  iprm_free(solver->work->kkt->iwork);
  iprm_free(solver->work->kkt->bwork);
  iprm_free(solver->work->kkt->fwork);
  iprm_free(solver->work->kkt->Li);
  iprm_free(solver->work->kkt->Lx);
  iprm_free(solver->work->kkt);

  iprm_free(solver->sol->x);
  iprm_free(solver->sol->xi);
  iprm_free(solver->sol->t);
  iprm_free(solver->sol->s);
  iprm_free(solver->sol);

  iprm_free(solver->work);
  iprm_free(solver->settings);
  iprm_free(solver);
  return 1;
}
