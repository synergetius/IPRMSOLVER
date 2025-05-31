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
void iprm_allocate_kkt(IPRMWorkspace* work){
  work->kkt->K = qoco_malloc(sizeof(QOCOCscMatrix));
  work->kkt->K->m = work->data->n + work->data->p;
  work->kkt->K->n = work->data->n + work->data->p;
  work->kkt->K->nnz = work->apkbpipe->ApkB->nnz + work->data->A->nnz + work->data->p;
  work->kkt->K->x = qoco_calloc(work->kkt->K->nnz, sizeof(QOCOFloat));
  work->kkt->K->i = qoco_calloc(work->kkt->K->nnz, sizeof(QOCOInt));
  work->kkt->K->p = qoco_calloc((work->kkt->K->n + 1), sizeof(QOCOInt));
}
void iprm_construct_kkt(IPRMSolver* solver){
  IPRMWorkspace* work = solver->work;
  QOCOCscMatrix* P = work->data->P;
  QOCOCscMatrix* PpkMtM = work->apkbpipe->ApkB;
  
  // 构建K的左上角块的稀疏结构并建立Preg+MtM到K的映射
  goMtMPipe(work->mtmpipe);
  goApkBPipe(work->apkbpipe, 1.0);
  for (int j = 0;j < PpkMtM->n;j++){
    work->kkt->K->p[j] = PpkMtM->p[j];
    for (int e = PpkMtM->p[j];e < PpkMtM->p[j + 1];e++){
      work->kkt->K->i[e] = PpkMtM->i[e];
      work->kkt->PpkMtMtoKKT[e] = e;
      work->kkt->K->x[e] = PpkMtM->x[e];
    }
  }
  work->kkt->K->p[PpkMtM->n] = PpkMtM->nnz;
  // 加入At和正则化块
  int nz = PpkMtM->nnz; 
  int col = work->data->n + 1; // 从此开始
  for (QOCOInt Atcol = 0; Atcol < work->data->At->n; ++Atcol) {
    QOCOInt nzadded = 0;
    for (
      QOCOInt k = work->data->At->p[Atcol]; 
      k < work->data->At->p[Atcol + 1];++k) {
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
    //下面这一句会遍历到col = work->data->n + work->data->At->n
    //从而把所有的K->p都记好
    work->kkt->K->p[col] = work->kkt->K->p[col - 1] + nzadded;
    col += 1;
  }
}




void print_K2(QOCOCscMatrix* K){
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
void iprm_initialize(IPRMSolver* solver){
  ///////// 打印发现initialize结果不符合原始残差和对偶残差接近零的条件，有BUG
  //printf("start init\n");
  IPRMWorkspace* work = solver->work;
  // 设置系数矩阵
  work->mu = 1.0;
  work->rho = 1.0;
  for (QOCOInt i = 0;i < work->data->m;++i){
    work->rhoy[i] = 1.0;
    work->y[i] = 1.0;
  }
  iprm_update_blocks(solver);
  //printf("done update blocks\n");
  // 设置等号右边的向量
  SpMtv(work->data->G, work->data->h, work->kkt->rhs);
  for (int i = 0; i < work->data->n; i++) {
    work->kkt->rhs[i] -= work->data->c[i];
  }
  for (int i = 0; i < work->data->p; i++) {
    work->kkt->rhs[work->data->n + i] = work->data->b[i];
  }
  //printf("done set rhs\n");
  // 求解x, t
  iprm_kkt_solve(solver);
  //printf("done solve\n");
  copy_arrayf(work->kkt->xt, work->x, work->data->n);
  copy_arrayf(&work->kkt->xt[work->data->n], work->t,
              work->data->p);
  //printf("#\n");
  // 计算s, xi
  SpMv(work->data->G, work->x, work->s);
  axpy(work->data->h, work->s, work->s, -1.0, work->data->m);

  SpMv(work->data->G, work->x, work->xi);
  axpy(work->xi, work->data->h, work->xi, -1, work->data->m);
  // 设置mu, rho的初始值
  work->mu = solver->settings->mu0;
  work->rho = solver->settings->rho0;
  //
  iprm_compute_zy(work);
  iprm_compute_psi_phi(solver);
  /*
  printf("c: ");
  for (int i = 0;i < work->data->n;i++){
    printf("%lf ", work->data->c[i]);
  }
  printf("\n");
  printf("b: ");
  for (int i = 0;i < work->data->p;i++){
    printf("%lf ", work->data->b[i]);
  }
  printf("\n");
  printf("primal res: ");
  for (int i = 0;i < work->data->n;i++){
    printf("%lf ", work->kkt->psi[i]);
  }
  printf("\n");
  printf("Ax - b: ");
  for (int i = 0;i < work->data->p;i++){
    printf("%lf ", work->kkt->psi[work->data->n + i]);
  }
  printf("\n");
  */
  solver->sol->iters = 0;
}

void iprm_update_blocks(IPRMSolver* solver){
  IPRMWorkspace* work = solver->work;
  // 先计算MtM，然后将其除以mu的结果和P分别加到K矩阵的左上角
  QOCOCscMatrix* G = work->data->G;
  // 计算M = rho * Y * G
  for (int j = 0;j < G->n;j++){
    for (int e = G->p[j];e < G->p[j + 1];e++){
      int i = G->i[e];
      work->mtmpipe->M->x[e] = work->rhoy[i] * G->x[e];
    }
  }
  // 计算MtM
  goMtMPipe(work->mtmpipe);
  // 计算Preg + 1/mu * MtM
  
  goApkBPipe(work->apkbpipe, 1 / work->mu);
  // 复制到K矩阵中
  for (int e = 0;e < work->apkbpipe->ApkB->nnz;e++){
    work->kkt->K->x[work->kkt->PpkMtMtoKKT[e]] 
    = work->apkbpipe->ApkB->x[e];
  }
}
void iprm_construct_kkt_rhs(IPRMWorkspace* work)
{
  IPRMProblemData* data = work->data;
  // 先将一、二行分别设为-r_k^d和-r_k^h
  copy_and_negate_arrayf(work->kkt->psi, work->kkt->rhs, data->n + data->p);
  // 计算括号内的式子
  for (int i = 0;i < data->m;i++){
    work->sbuff[i] = 
      work->rhoy[i] * work->rhoy[i] * 
        (work->kkt->psi[data->n + data->p + i]
        +work->kkt->psi[data->n + data->p + data->m + i]) / work->mu
      + work->rho * work->kkt->psi[data->n + data->p + i]
      + work->Dmu * work->rhoy[i] / work->mu;
  }
  // 计算G^T乘以括号内的式子
  SpMtv(data->G, work->sbuff, work->xbuff);
  // 将-G^T(...)加到第一行上
  axpy(work->xbuff, work->kkt->rhs, work->kkt->rhs, -1.0, data->n);

}

void iprm_kkt_solve(IPRMSolver* solver){
  IPRMWorkspace* work = solver->work;
  IPRMProblemData* data = work->data;
  // Factor KKT matrix.
  //printf("start iprm kkt solve\n");
  QDLDL_factor(work->kkt->K->n, work->kkt->K->p, work->kkt->K->i,
               work->kkt->K->x, work->kkt->Lp, work->kkt->Li, work->kkt->Lx,
               work->kkt->D, work->kkt->Dinv, work->kkt->Lnz, work->kkt->etree,
               work->kkt->bwork, work->kkt->iwork, work->kkt->fwork,
               solver->work->kkt->p, solver->work->data->n,
               solver->settings->kkt_dynamic_reg);
  //printf("done factor\n");
  IPRMKKT* kkt = solver->work->kkt;

  // 将rhs排序，结果分别在xtbuff1和rhs中各存一份
  for (QOCOInt k = 0; k < kkt->K->n; ++k) {
    kkt->xtbuff1[k] = kkt->rhs[kkt->p[k]];
  }
  copy_arrayf(kkt->xtbuff1, kkt->rhs, kkt->K->n);
  // 解2x2方程组，得到排序后的结果在xtbuff1中
  QDLDL_solve(kkt->K->n, kkt->Lp, kkt->Li, kkt->Lx, kkt->Dinv, kkt->xtbuff1);
  // Iterative refinement.
  for (QOCOInt i = 0;i < solver->settings->iter_ref_iters;++i){
    //将上次迭代所求解得到的排序后结果（xtbuff1）转换成未排序的，存入xt
    for (QOCOInt k = 0; k < kkt->K->n; ++k) {
      kkt->xt[kkt->p[k]] = kkt->xtbuff1[k];
    }
    
    /* 计算K * xt（未正则化的K），结果存入xtbuff2*/
    // 计算P * x，并抵消静态正则化的影响
    USpMv(data->P, kkt->xt, kkt->xtbuff2);
    axpy(kkt->xt, kkt->xtbuff2, kkt->xtbuff2, -solver->settings->kkt_static_reg, data->n);
    // 计算MtM * x并将其1/mu倍加到xtbuff2第一行
    USpMv(work->mtmpipe->MtM, kkt->xt, work->xbuff);
    axpy(work->xbuff, kkt->xtbuff2, kkt->xtbuff2, 1 / work->mu, data->n);
    if (data->p > 0){
      // 计算 A^T * t并加到xtbuff2第一行
      SpMtv(data->A, &kkt->xt[data->n], work->xbuff);
      axpy(work->xbuff, kkt->xtbuff2, kkt->xtbuff2, 1.0, data->n);
      // 计算 A * x 并存到xtbuff2第二行
      SpMv(data->A, kkt->xt, &kkt->xtbuff2[data->n]);
    }
    // 将xtbuff2中未排序结果 排序并转移到xt
    for (QOCOInt k = 0; k < kkt->K->n; ++k) {
      kkt->xt[k] = kkt->xtbuff2[kkt->p[k]];
    }
    ///////
    // rhs - K * （上一代的）xt 计算后存储在xt 中
    for (QOCOInt k = 0; k < kkt->K->n; ++k) {
      kkt->xt[k] = kkt->rhs[k] - kkt->xt[k]; // r - Kd
    }
    // 求解迭代方程组，结果（矫正项）存入kkt->xt
    QDLDL_solve(kkt->K->n, kkt->Lp, kkt->Li, kkt->Lx, kkt->Dinv, kkt->xt);
    // 将矫正项加到上一代的解中（已排序），得到新一代解（已排序）
    axpy(kkt->xtbuff1, kkt->xt, kkt->xtbuff1, 1.0, kkt->K->n);
  }
  // 得到最终未排序的解（存储在xt中）
  for (QOCOInt i = 0; i < kkt->K->n; ++i) {
    kkt->xt[kkt->p[i]] = kkt->xtbuff1[i];
  }
}