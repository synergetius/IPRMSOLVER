/**
 * @file kkt.c
 * @author Govind M. Chari <govindchari1@gmail.com>
 *
 * @section LICENSE
 *
 * Copyright (c) 2024, Govind M. Chari
 * This source code is licensed under the BSD 3-Clause License
 */

#include "kkt.h"
#include "utils.h"


void iprm_allocate_kkt(IPRMWorkspace* work)
{
  work->kkt->K = qoco_malloc(sizeof(QOCOCscMatrix));
  work->kkt->K->m = work->data->n + work->data->m + work->data->p;
  work->kkt->K->n = work->data->n + work->data->m + work->data->p;
  work->kkt->K->nnz = 
      work->data->P->nnz + // 静态正则化后的P
      work->data->A->nnz + // 
      work->data->G->nnz + // 
      work->data->p + // (2, 2)块静态正则化
      work->data->m; // (3, 3)块
  work->kkt->K->x = qoco_calloc(work->kkt->K->nnz, sizeof(QOCOFloat));
  work->kkt->K->i = qoco_calloc(work->kkt->K->nnz, sizeof(QOCOInt));
  work->kkt->K->p = qoco_calloc((work->kkt->K->n + 1), sizeof(QOCOInt));
}

void iprm_construct_kkt(IPRMSolver* solver)
{
  IPRMWorkspace* work = solver->work;
  QOCOInt nz = 0;
  QOCOInt col = 0;
  // Add P block
  for (QOCOInt k = 0; k < work->data->P->nnz; ++k) {
    work->kkt->PregtoKKT[k] = nz;
    work->kkt->K->x[nz] = work->data->P->x[k];
    work->kkt->K->i[nz] = work->data->P->i[k];
    nz += 1;
  }
  for (QOCOInt k = 0; k < work->data->P->n + 1; ++k) {
    work->kkt->K->p[col] = work->data->P->p[k];
    col += 1;
  }
  // Add A^T block
  for (QOCOInt Atcol = 0; Atcol < work->data->At->n; ++Atcol) {
    QOCOInt nzadded = 0;
    for (QOCOInt k = work->data->At->p[Atcol]; k < work->data->At->p[Atcol + 1];
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
    work->kkt->K->i[nz] = work->data->n + Atcol;
    nz += 1;
    nzadded += 1;
    work->kkt->K->p[col] = work->kkt->K->p[col - 1] + nzadded;
    col += 1;
  }
  // Add G^T block and (3, 3) diagonal block
  // 此处赋值仅作参考（先保留，因为不确定是否和AMD排序等预处理有关），实际不会用到这里的值，因为μG^T和-μZ^2都是动态变化的
  QOCOInt diag = 0;
  for (QOCOInt Gtcol = 0; Gtcol < work->data->m; ++Gtcol) {
    QOCOInt nzadded = 0;
    for (QOCOInt k = work->data->Gt->p[Gtcol]; k < work->data->Gt->p[Gtcol + 1];
         ++k) {
      work->kkt->GtoKKT[k] = nz;
      work->kkt->K->x[nz] = work->data->Gt->x[k]; //
      work->kkt->K->i[nz] = work->data->Gt->i[k];
      nz += 1;
      nzadded += 1;
    }
    work->kkt->K->x[nz] = -1; //
    work->kkt->K->i[nz] = work->data->n + work->data->p + Gtcol;
    work->kkt->K->p[col] = work->kkt->K->p[col - 1] + nzadded + 1;

    work->kkt->ZtoKKT[diag] = nz;
    diag++;
    nz += 1;
    col += 1;
  }
}
/*
void print_K(QOCOCscMatrix* K){
  for (int j = 0;j < K->n;j++){
    int t = K->p[j];
    for (int i = 0;i < K->m;i++){
      if (t < K->p[j + 1] && i == K->i[t]){
        printf("%.2lf\t", K->x[t]);
        t += 1;
      }
      else{
        printf("0\t");
      }
    }
    printf("\n");
  }
}
*/
void iprm_initialize(IPRMSolver* solver){
  IPRMWorkspace* work = solver->work;

  for (QOCOInt i = 0; i < work->data->m; ++i) {
    work->kkt->K->x[work->kkt->ZtoKKT[i]] = -1.0;
  }
  for (QOCOInt i = 0; i < work->data->Gt->nnz; ++i) {
    work->kkt->K->x[work->kkt->GtoKKT[i]] = work->data->Gt->x[i];
    
  }
  work->mu = 1.0; // 要与K矩阵一致
  for (QOCOInt i = 0;i < work->data->m;++i){
    work->z[i] = 1.0;
  }

  

  QOCOInt idx = 0;
  for (idx = 0; idx < work->data->n; ++idx) {
    work->kkt->rhs[idx] = -work->data->c[idx];
  }
  for (QOCOInt i = 0; i < work->data->p; ++i) {
    work->kkt->rhs[idx] = work->data->b[i];
    idx += 1;
  }
  for (QOCOInt i = 0; i < work->data->m; ++i) {
    work->kkt->rhs[idx] = work->data->h[i];
    idx += 1;
  }

  iprm_kkt_solve(solver);

  copy_arrayf(work->kkt->xts, work->x, work->data->n);
  copy_arrayf(&work->kkt->xts[work->data->n], work->t,
              work->data->p);
  copy_arrayf(
      &work->kkt->xts[work->data->n + work->data->p],
      work->s, work->data->m);
  copy_and_negate_arrayf(
      &work->kkt->xts[work->data->n + work->data->p],
      work->xi, work->data->m); // 使用QOCO的初始化
  /*
  printf("x: %lf\n", norm_2(work->x, work->data->n));
  printf("t: %lf\n", norm_2(work->t, work->data->p));
  printf("s: %lf\n", norm_2(work->s, work->data->m));
  printf("xi: %lf\n", norm_2(work->xi, work->data->m));
  */
  work->mu = solver->settings->mu0;
  work->rho = solver->settings->rho0;
}






void iprm_update_blocks(IPRMSolver* solver){
  // Block (1, 3)
  IPRMWorkspace* work = solver->work;
  for (QOCOInt i = 0; i < work->data->Gt->nnz; ++i) {
    work->kkt->K->x[work->kkt->GtoKKT[i]] = work->mu * work->data->Gt->x[i];
  }
  for (QOCOInt j = 0; j < work->data->m;j++){
    work->kkt->K->x[work->kkt->ZtoKKT[j]] = -work->mu * work->z[j] * work->z[j] - solver->settings->kkt_static_reg;
  }
}
void iprm_construct_kkt_rhs(IPRMWorkspace* work)
{
  copy_and_negate_arrayf(work->kkt->psi, work->kkt->rhs, work->data->n + work->data->p);
  for (QOCOInt i = 0;i < work->data->m;++i){
    work->kkt->rhs[work->data->n + work->data->p + i] = -(
        work->rho * (work->z[i] + work->y[i]) * work->z[i] 
            * work->kkt->psi[work->data->n + work->data->p + i]
        + work->Dmu * work->z[i]
        + work->mu * work->kkt->psi[work->data->n + work->data->p + work->data->m + i]
    );
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
               solver->work->kkt->p, solver->work->data->n,
               solver->settings->kkt_dynamic_reg);
  IPRMKKT* kkt = solver->work->kkt;

  // Permute and store rhs
  for (QOCOInt k = 0; k < kkt->K->n; ++k) {
    kkt->npm_buff[k] = kkt->rhs[kkt->p[k]];
  }
  copy_arrayf(kkt->npm_buff, kkt->rhs, kkt->K->n); // 这里rhs就是排过序的了
  QDLDL_solve(kkt->K->n, kkt->Lp, kkt->Li, kkt->Lx, kkt->Dinv, kkt->npm_buff);
  // Iterative refinement.
  for (QOCOInt i = 0;i < solver->settings->iter_ref_iters;++i){
    //计算 r = rhs - K Dxts
    for (QOCOInt k = 0; k < kkt->K->n; ++k) {
      kkt->xts[kkt->p[k]] = kkt->npm_buff[k];
    }
    // 现在xts是未排列的，要根据K分块计算
    // 结果存入npm_buff2
    USpMv(data->P,kkt->xts, kkt->npm_buff2);
    // 抵消静态正则化的影响
    axpy(kkt->xts, kkt->npm_buff2, kkt->npm_buff2, -solver->settings->kkt_static_reg, data->n);
    if (data->p > 0){
      SpMtv(data->A, &kkt->xts[data->n],work->xbuff);
      axpy(work->xbuff, kkt->npm_buff2, kkt->npm_buff2, 1.0, data->n);
      SpMv(data->A, kkt->xts, &kkt->npm_buff2[data->n]);
    }
    if (data->m > 0){
      SpMtv(data->G, &kkt->xts[data->n + data->p], work->xbuff);
      axpy(work->xbuff, kkt->npm_buff2, kkt->npm_buff2, work->mu, data->n);
      SpMv(data->G, kkt->xts, &kkt->npm_buff2[data->n + data->p]);
      
      for (int j = 0;j < data->m;j++){
        kkt->npm_buff2[data->n + data->p + j] = work->mu * (
          kkt->npm_buff2[data->n + data->p + j] - 
          work->z[j] * work->z[j] * kkt->xts[data->n + data->p + j]
        ); 
      }
    }
    for (QOCOInt k = 0; k < kkt->K->n; ++k) {
      kkt->xts[k] = kkt->npm_buff2[kkt->p[k]]; // 存储为排序后顺序
    }
    for (QOCOInt k = 0; k < kkt->K->n; ++k) {
      kkt->xts[k] = kkt->rhs[k] - kkt->xts[k]; // r - Kd
    }
    QDLDL_solve(kkt->K->n, kkt->Lp, kkt->Li, kkt->Lx, kkt->Dinv, kkt->xts);
    axpy(kkt->npm_buff, kkt->xts, kkt->npm_buff, 1.0, kkt->K->n);
    // 得到新的已排序的解
  }

  for (QOCOInt i = 0; i < kkt->K->n; ++i) {
    kkt->xts[kkt->p[i]] = kkt->npm_buff[i]; // 得到最终的解（未排序）
  }
}