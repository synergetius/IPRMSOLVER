

#include "kkt.h"
#include "utils.h"

void print_K(IPRMCscMatrix* K){
  // csc_symperm所得到的同一列的非零元素的行可以不按顺序！所以需要完整遍历才能打印全
  for (int j = 0;j < K->n;j++){
    int t = K->p[j];
    for (int i = 0;i < K->m;i++){
      int fl = 0;
      for (int t = K->p[j];t < K->p[j + 1];t++){
        if(K->i[t] == i){
          fl = 1;
          printf("%.2lf\t", K->x[t]);
          break;
        }
      }
      if (!fl){
        printf("0\t");
      }
      /*
      if (t < K->p[j + 1] && i == K->i[t]){
        printf("%.2lf\t", K->x[t]);
        t += 1;
      }
      else{
        printf("0\t");
      }
      */
    }
    printf("\n");
  }
}
void iprm_allocate_kkt(IPRMWorkspace* work)
{
  work->kkt->K = iprm_malloc(sizeof(IPRMCscMatrix));
  work->kkt->K->m = work->data->n + work->data->nsoc + work->data->nsoc + work->data->p + work->data->m;
  work->kkt->K->n = work->data->n + work->data->nsoc + work->data->nsoc + work->data->p + work->data->m;
  work->kkt->K->nnz = 
      work->data->P->nnz + // (1, 1): Preg
      work->data->nsoc + // (2, 2): I_{nsoc}
      work->data->nsoc + // (3, 3): -I_{nsoc}
      work->data->A->nnz + // (1, 4): A^T
      work->data->p + // (4, 4)正则化块
      work->data->G->nnz + // (1, 5), (1, 6)
      (work->data->m - work->data->l) + // (2, 6)
      work->data->nsoc + // (3, 6) e_i只有首个元素非零
      work->data->m; // (5, 5)和(6, 6)块
  work->kkt->K->x = iprm_calloc(work->kkt->K->nnz, sizeof(IPRMFloat));
  work->kkt->K->i = iprm_calloc(work->kkt->K->nnz, sizeof(IPRMInt));
  work->kkt->K->p = iprm_calloc((work->kkt->K->n + 1), sizeof(IPRMInt));
}

void iprm_construct_kkt(IPRMSolver* solver)
{
  IPRMWorkspace* work = solver->work;
  IPRMInt nz = 0;
  IPRMInt col = 0;
  // Add P block
  for (IPRMInt k = 0; k < work->data->P->nnz; ++k) {
    work->kkt->PregtoKKT[k] = nz;
    work->kkt->K->x[nz] = work->data->P->x[k];
    work->kkt->K->i[nz] = work->data->P->i[k];
    nz += 1;
  }
  for (IPRMInt k = 0; k < work->data->P->n + 1; ++k) {
    work->kkt->K->p[col] = work->data->P->p[k];
    col += 1;
  }
  // Add (2, 2) = I_{nsoc} block
  for (IPRMInt Icol = 0; Icol < work->data->nsoc; ++Icol){
    work->kkt->K->x[nz] = 1.0;
    work->kkt->K->i[nz] = work->data->n + Icol;
    nz += 1;
    work->kkt->K->p[col] = work->kkt->K->p[col - 1] + 1;
    col += 1;
  }
  // Add (4, 4) = -I_{nsoc} block
  for (IPRMInt nIcol = 0; nIcol < work->data->nsoc; ++nIcol){
    work->kkt->K->x[nz] = -1.0;
    work->kkt->K->i[nz] = work->data->n + work->data->nsoc + nIcol;
    nz += 1;
    work->kkt->K->p[col] = work->kkt->K->p[col - 1] + 1;
    col += 1;
  }
  // Add A^T block
  for (IPRMInt Atcol = 0; Atcol < work->data->At->n; ++Atcol) {
    IPRMInt nzadded = 0;
    for (IPRMInt k = work->data->At->p[Atcol]; k < work->data->At->p[Atcol + 1];
         ++k) {
      // If the nonzero is in row i of A then add
      work->kkt->AtoKKT[k] = nz;
      work->kkt->K->x[nz] = work->data->At->x[k];
      work->kkt->K->i[nz] = work->data->At->i[k];
      nz += 1;
      nzadded += 1;
    }

    // Add -e * Id regularization.
    work->kkt->K->x[nz] = -solver->settings->kkt_static_reg;
    work->kkt->K->i[nz] = work->data->n + work->data->nsoc + work->data->nsoc + Atcol;
    nz += 1;
    nzadded += 1;
    work->kkt->K->p[col] = work->kkt->K->p[col - 1] + nzadded;
    col += 1;
  }
  // Add non-negative orthant part of sqrt(mu) * G^T and Block (5, 5)
  IPRMInt diag = 0;
  for (IPRMInt Gtcol = 0; Gtcol < work->data->l; ++Gtcol) {
    // Counter for number of nonzeros from G added to this column of KKT matrix
    IPRMInt nzadded = 0;
    for (IPRMInt k = work->data->Gt->p[Gtcol]; k < work->data->Gt->p[Gtcol + 1];
         ++k) {
      work->kkt->GtoKKT[k] = nz;
      work->kkt->K->x[nz] = work->data->Gt->x[k];
      work->kkt->K->i[nz] = work->data->Gt->i[k];
      nz += 1;
      nzadded += 1;
    }
    // -Z_0^2 暂时以-I占位
    work->kkt->K->x[nz] = -1.0;
    work->kkt->K->i[nz] = work->data->n + work->data->nsoc + work->data->nsoc + work->data->p + Gtcol;
    work->kkt->K->p[col] = work->kkt->K->p[col - 1] + nzadded + 1;

    // Mapping from NT matrix entries to KKT matrix entries.
    work->kkt->B55toKKT[diag] = nz;
    diag++;
    nz += 1;
    col += 1;
  }
  // Add second-order cone parts of G^T and other blocks of column 6
  IPRMInt idx = work->data->l;
  for (IPRMInt c = 0; c < work->data->nsoc; ++c) {
    for (IPRMInt Gtcol = idx; Gtcol < idx + work->data->q[c]; ++Gtcol) {
      // Loop over columns of G

      // nzadded 统计这一列中添加的元素数量
      IPRMInt nzadded = 0;
      for (IPRMInt k = work->data->Gt->p[Gtcol];
           k < work->data->Gt->p[Gtcol + 1]; ++k) {
        work->kkt->GtoKKT[k] = nz;
        work->kkt->K->x[nz] = work->data->Gt->x[k];
        work->kkt->K->i[nz] = work->data->Gt->i[k];
        nz += 1;
        nzadded += 1;
      }
      // 考虑mapping的格式：
      // P, A^T, G^T 都直接使用非零元素序数(nz)的映射
      // -Z_0^2使用Z序数到非零元素序数的映射
      // (2, 6), (6, 6)使用Z序数到非零元素序数的映射
      // (3, 6) 使用SOC序数到非零元素序数的映射

      // Add block (2, 6)
      work->kkt->B26toKKT[Gtcol] = nz;
      if (Gtcol == idx){ // 暂时用soc的单位e_i占位
        work->kkt->K->x[nz] = -iprm_sqrt(2);
      }
      else{
        work->kkt->K->x[nz] = 0;//5.05;
      }
      work->kkt->K->i[nz] = work->data->n + c;
      nz += 1;
      nzadded += 1;
      // Add block (3, 6)
      if (Gtcol == idx){ // 是该二阶锥对应的首个元素
        work->kkt->B36toKKT[c] = nz;
        work->kkt->K->x[nz] = iprm_sqrt(2);
        work->kkt->K->i[nz] = work->data->n + work->data->nsoc + c;
        nz += 1;
        nzadded += 1;
      }
      // Add block (6, 6)
      work->kkt->B66toKKT[Gtcol] = nz;
      work->kkt->K->x[nz] = -1;
      work->kkt->K->i[nz] = work->data->n + work->data->p + work->data->nsoc + work->data->nsoc + Gtcol;
      nz += 1;
      nzadded += 1;

      work->kkt->K->p[col] = work->kkt->K->p[col - 1] + nzadded;
      col += 1;
    }
    idx += work->data->q[c];
  } 
  
}



void iprm_initialize(IPRMSolver* solver){
  IPRMWorkspace* work = solver->work;
  // 应该把z设成单位，这样对应于初始化所需的3x3分块形式
  for (IPRMInt i = 0;i < work->data->l;++i){
    work->z[i] = 1.0;
  }
  IPRMInt idx = work->data->l;
  for (IPRMInt i = 0;i < work->data->nsoc;++i){
    for (IPRMInt j = 1;j < work->data->q[i];j++){
      work->z[idx + j] = 0;
    }
    work->z[idx] = 1.0;
    idx += work->data->q[i];
  }
  work->mu = 1.0;
  
  iprm_update_blocks(solver);
  
  copy_and_negate_arrayf(work->data->c, work->kkt->rhs, work->data->n);
  for (IPRMInt i = work->data->n; i < work->data->n + work->data->nsoc + work->data->nsoc; ++i) {
    work->kkt->rhs[i] = 0;
  }
  copy_arrayf(work->data->b, &work->kkt->rhs[work->data->n + work->data->nsoc + work->data->nsoc], work->data->p);
  copy_arrayf(work->data->h, &work->kkt->rhs[work->data->n + work->data->nsoc + work->data->nsoc + work->data->p], work->data->m);
  /*
  for (idx = 0; idx < work->data->n; ++idx) {
    work->kkt->rhs[idx] = -work->data->c[idx];
  }
  
  for (IPRMInt i = 0; i < work->data->nsoc + work->data->nsoc; ++i) {
    work->kkt->rhs[idx] = 0;
    idx += 1;
  }
  for (IPRMInt i = 0; i < work->data->p; ++i) {
    work->kkt->rhs[idx] = work->data->b[i];
    idx += 1;
  }
  for (IPRMInt i = 0; i < work->data->m; ++i) {
    work->kkt->rhs[idx] = work->data->h[i];
    idx += 1;
  }
    */
  iprm_kkt_solve(solver);
  //printf("b\n");
  copy_arrayf(work->kkt->w, work->x, work->data->n);
  copy_arrayf(&work->kkt->w[work->data->n + work->data->nsoc + work->data->nsoc], work->t,
              work->data->p);
  copy_arrayf(
      &work->kkt->w[work->data->n + work->data->nsoc + work->data->nsoc + work->data->p],
      work->s, work->data->m);
  //printf("c\n");
  work->mu = solver->settings->mu0;
  work->rho = solver->settings->rho0;
  //printf("d\n");
  SpMv(work->data->G, work->x, work->xi);
  axpy(work->xi, work->data->h, work->xi, -1, work->data->m);
  
  /*
  printf("x: %lf\n", norm_2(work->x, work->data->n));
  printf("t: %lf\n", norm_2(work->t, work->data->p));
  printf("s: %lf\n", norm_2(work->s, work->data->m));
  printf("xi: %lf\n", norm_2(work->xi, work->data->m));
  */
  solver->sol->iters = 0;
}






void iprm_update_blocks(IPRMSolver* solver){
  
  IPRMWorkspace* work = solver->work;
  // Block (1, 5), (1, 6)
  for (IPRMInt i = 0; i < work->data->Gt->nnz; ++i) {
    work->kkt->K->x[work->kkt->GtoKKT[i]] = iprm_sqrt(work->mu) * work->data->Gt->x[i];
  }
  //printf("done(1, 5)(1, 6)\n");
  // Block (2, 6)
  
  IPRMInt idx = work->data->l;
  for (IPRMInt c = 0; c < work->data->nsoc;++c){
    for (int i = idx;i < idx + work->data->q[c];i++){
      work->kkt->K->x[work->kkt->B26toKKT[i]] = -iprm_sqrt(2.0) * work->z[i];
      //printf("#: %d %.10e\n", work->kkt->B26toKKT[i], work->kkt->K->x[work->kkt->B26toKKT[i]]);
    }
    idx += work->data->q[c];
  }
  //printf("done(2, 6)\n");
  // block (3, 6)
  idx = work->data->l;
  for (IPRMInt c = 0; c < work->data->nsoc;++c){
    // 计算det(z_c)
    IPRMFloat det = soc_det(&work->z[idx], work->data->q[c]);
    //printf("det = %.12e\n", det);
    work->kkt->K->x[work->kkt->B36toKKT[c]] = iprm_sqrt(2 * det);
    idx += work->data->q[c];
  }
  //printf("done(3, 6)\n");
  // block (5, 5)
  for (IPRMInt i = 0;i < work->data->l;i++){
    //work->kkt->K->x[work->kkt->B55toKKT[i]] = -work->z[i] * work->z[i];
    work->kkt->K->x[work->kkt->B55toKKT[i]] = -work->z[i] * work->z[i] - solver->settings->kkt_static_reg;
  }
  //printf("done(5, 5)\n");
  // block (6, 6)
  idx = work->data->l;
  for (IPRMInt c = 0; c < work->data->nsoc;++c){
    // 计算det(z_c)
    IPRMFloat det = soc_det(&work->z[idx], work->data->q[c]);
    //printf("det = %.12e\n", det);
    for (int i = idx;i < idx + work->data->q[c];i++){
      //work->kkt->K->x[work->kkt->B66toKKT[i]] = -det;
      work->kkt->K->x[work->kkt->B66toKKT[i]] = -det - solver->settings->kkt_static_reg; // det加正则化
    }
    idx += work->data->q[c];
  }
  //printf("done(6, 6)\n");
}
void iprm_construct_kkt_rhs(IPRMWorkspace* work)
{
  
  copy_and_negate_arrayf(work->kkt->psi, work->kkt->rhs, work->data->n);
  for (IPRMInt i = work->data->n;i < work->data->n + work->data->nsoc + work->data->nsoc;++i){
    work->kkt->rhs[i] = 0;
  }
  copy_and_negate_arrayf(&work->kkt->psi[work->data->n], &work->kkt->rhs[work->data->n + work->data->nsoc + work->data->nsoc], work->data->p);
  // 
  for (IPRMInt j = 0;j < work->data->l;j++){
    work->kkt->rhs[work->data->n + work->data->nsoc + work->data->nsoc + work->data->p + j] = -(
      work->rho * work->z[j] * work->z[j] * work->kkt->psi[work->data->n + work->data->p + j] / iprm_sqrt(work->mu)
      + iprm_sqrt(work->mu) * work->kkt->psi[work->data->n + work->data->p + j]
      + work->Dmu * work->z[j] / iprm_sqrt(work->mu)
      + iprm_sqrt(work->mu) * work->kkt->psi[work->data->n + work->data->p + work->data->m + j]
    );
  }
  //
  IPRMInt idx = work->data->l;
  for (IPRMInt c = 0;c < work->data->nsoc;c++){
    IPRMFloat det = soc_det(&work->z[idx], work->data->q[c]);
    IPRMFloat zr = dot(&work->z[idx], &work->kkt->psi[work->data->n + work->data->p + idx], work->data->q[c]);
    for (int j = idx;j < idx + work->data->q[c];j++){
      work->kkt->rhs[work->data->n + work->data->nsoc + work->data->nsoc + work->data->p + j] = -(
        (work->rho * det / iprm_sqrt(work->mu) + iprm_sqrt(work->mu)) * work->kkt->psi[work->data->n + work->data->p + j]
        + (2 * work->rho * zr + work->Dmu) * work->z[j] / iprm_sqrt(work->mu) 
        + iprm_sqrt(work->mu) * work->kkt->psi[work->data->n + work->data->p + work->data->m + j]
      );
    }
    // 补上e_i项的系数（括号内减相当于去括号后的加）
    work->kkt->rhs[work->data->n + work->data->nsoc + work->data->nsoc + work->data->p + idx] +=
      2 * work->rho * det * work->kkt->psi[work->data->n + work->data->p + idx] / iprm_sqrt(work->mu);
    idx += work->data->q[c];
  }
}
void iprm_kkt_multiply(IPRMSolver* solver, IPRMFloat* src, IPRMFloat* tar){
  // 检查了kkt_solve没发现问题，而rhs的重构误差仍然很大，怀疑是此函数中有问题
  // 计算K * src （K要使用无正则化的版本）
  IPRMWorkspace* work = solver->work;
  IPRMProblemData* data = solver->work->data;
  IPRMInt idx;
  // Row 1
  USpMv(work->data->P, src, tar);
  axpy(src, tar, tar, -solver->settings->kkt_static_reg, data->n);
  if (data->p > 0){
    SpMtv(data->A, &src[data->n + data->nsoc + data->nsoc], work->xbuff);
    axpy(work->xbuff, tar, tar, 1.0, data->n);
  }
  if (data->m > 0){
    SpMtv(data->G, &src[data->n + data->nsoc + data->nsoc + data->p], work->xbuff);
    axpy(work->xbuff, tar, tar, iprm_sqrt(work->mu), data->n);
  }
  // Row 2
  copy_arrayf(&src[data->n], &tar[data->n], data->nsoc);
  idx = data->l;
  for (int c = 0;c < data->nsoc;c++){
    IPRMFloat d = dot(&work->z[idx], &src[data->n + data->nsoc + data->nsoc + data->p + idx], data->q[c]);
    tar[data->n + c] -= iprm_sqrt(2.0) * d;
    idx += data->q[c];
  }
  // Row 3
  copy_and_negate_arrayf(&src[data->n + data->nsoc], &tar[data->n + data->nsoc], data->nsoc);
  idx = data->l;
  for (int c = 0;c < data->nsoc;c++){
    IPRMFloat det = soc_det(&work->z[idx], data->q[c]);
    tar[data->n + data->nsoc + c] += iprm_sqrt(2.0 * det) * src[data->n + data->nsoc + data->nsoc + data->p + idx];
    idx += data->q[c];
  }
  // Row 4
  if (data->p > 0){
    SpMv(data->A, src, &tar[data->n + data->nsoc + data->nsoc]);
  }
  // Row 5, 6
  if (data->m > 0){
    SpMv(data->G, src, work->sbuff);
    for (int i = 0;i < data->m;i++){
      tar[data->n + data->nsoc + data->nsoc + data->p + i] = iprm_sqrt(work->mu) * work->sbuff[i];
    }
  }
  // Row 5
  for (int i = 0; i < data->l;i++){
    tar[data->n + data->nsoc + data->nsoc + data->p + i] -= work->z[i] * work->z[i] * src[data->n + data->nsoc + data->nsoc + data->p + i];
  }
  // Row 6
  idx = data->l;
  for (int c = 0;c < data->nsoc;c++){
    // item 2
    axpy(&work->z[idx], 
      &tar[data->n + data->nsoc + data->nsoc + data->p + idx], 
      &tar[data->n + data->nsoc + data->nsoc + data->p + idx], 
      -iprm_sqrt(2.0) * src[data->n + c], data->q[c]);
    IPRMFloat det = soc_det(&work->z[idx], data->q[c]);
    // item 3
    tar[data->n + data->nsoc + data->nsoc + data->p + idx] += iprm_sqrt(2 * det) * src[data->n + data->nsoc + c];
    // item 4
    axpy(&src[data->n + data->nsoc + data->nsoc + data->p + idx], 
      &tar[data->n + data->nsoc + data->nsoc + data->p + idx],
      &tar[data->n + data->nsoc + data->nsoc + data->p + idx],
      -det, data->q[c]
    );
    idx += data->q[c];
  }

}
void iprm_kkt_solve(IPRMSolver* solver){
  IPRMWorkspace* work = solver->work;
  IPRMProblemData* data = work->data;
  // Factor KKT matrix.
  QDLDL_factor(work->kkt->K->n, work->kkt->K->p, work->kkt->K->i,
               work->kkt->K->x, work->kkt->Lp, work->kkt->Li, work->kkt->Lx,
               work->kkt->D, work->kkt->Dinv, work->kkt->Lnz, work->kkt->etree,
               work->kkt->bwork, work->kkt->iwork, work->kkt->fwork,
               solver->work->kkt->p, solver->work->data->n + solver->work->data->nsoc,
               solver->settings->kkt_dynamic_reg);
    // 这里用到work->kkt->p
  //printf("factored\n");
	// printf("||L||_inf = %.12e\n", inf_norm(work->kkt->Lx, work->kkt->Lp[work->kkt->K->n]));
  IPRMKKT* kkt = solver->work->kkt;

  // Permute and store rhs
  for (IPRMInt k = 0; k < kkt->K->n; ++k) {
    kkt->kktbuff1[k] = kkt->rhs[kkt->p[k]];
  }
  copy_arrayf(kkt->kktbuff1, kkt->rhs, kkt->K->n); // 这里rhs就是排过序的了
  //printf("copied\n");
  QDLDL_solve(kkt->K->n, kkt->Lp, kkt->Li, kkt->Lx, kkt->Dinv, kkt->kktbuff1);
  // Iterative refinement.
  //printf("solved\n");
  for (IPRMInt i = 0;i < solver->settings->iter_ref_iters;++i){
    //计算 r = rhs - K w
    for (IPRMInt k = 0; k < kkt->K->n; ++k) {
      kkt->w[kkt->p[k]] = kkt->kktbuff1[k];
    }
    // 现在w是未排列的，要根据K分块计算
    // 结果存入kktbuff2
    iprm_kkt_multiply(solver, kkt->w, kkt->kktbuff2);
    
    for (IPRMInt k = 0; k < kkt->K->n; ++k) {
      kkt->w[k] = kkt->kktbuff2[kkt->p[k]]; // 存储为排序后顺序
    }
    for (IPRMInt k = 0; k < kkt->K->n; ++k) {
      kkt->w[k] = kkt->rhs[k] - kkt->w[k]; // rhs - Kw
    }
    QDLDL_solve(kkt->K->n, kkt->Lp, kkt->Li, kkt->Lx, kkt->Dinv, kkt->w);
    axpy(kkt->kktbuff1, kkt->w, kkt->kktbuff1, 1.0, kkt->K->n);
    // 得到新的已排序的解
  }

  for (IPRMInt i = 0; i < kkt->K->n; ++i) {
    kkt->w[kkt->p[i]] = kkt->kktbuff1[i]; // 得到最终的解（未排序）
  }
  // 检验所得解是否是方程组的解
  for (IPRMInt k = 0; k < kkt->K->n; ++k) {
    kkt->kktbuff1[kkt->p[k]] = kkt->rhs[k]; // 未排序的rhs
  }
  iprm_kkt_multiply(solver, kkt->w, kkt->kktbuff2);
  printf("reconstructed rhs: %.10e\n", inf_norm(kkt->kktbuff2, kkt->K->n));
  axpy(kkt->kktbuff1, kkt->kktbuff2, kkt->kktbuff2, -1, kkt->K->n);
  printf("rhs error: %.10e\n", inf_norm(kkt->kktbuff2, kkt->K->n));
}