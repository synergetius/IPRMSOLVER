

#include "cone.h"
#include "utils.h"
#include <math.h>


void iprm_compute_zy(IPRMWorkspace* work){
  // 计算z, y
  // 不使用任何buff
  for (QOCOInt j = 0;j < work->data->m;j++){
    QOCOFloat d = work->s[j] - work->rho * work->xi[j];
    QOCOFloat e = qoco_sqrt(d * d + 4 * work->rho * work->mu);
    work->z[j] = (e - d) / (2 * work->rho);
    work->rhoy[j] = (e + d) / 2; ///// 计算一个有用的中间值
    work->y[j] = (e + d) / (2 * work->rho);   /////////// 考虑这里是否有必要用safe_div
  }

}
void iprm_compute_psi_phi(IPRMSolver* solver)
{
  // 计算psi和phi（KKT方程残差）
  // 仅使用work->xbuff
  IPRMWorkspace* work = solver->work;
  IPRMProblemData* data = solver->work->data;

  USpMv(data->P, work->x, work->kkt->psi);
  for (int j = 0; j < data->n; ++j) {
    work->kkt->psi[j] += data->c[j] - solver->settings->kkt_static_reg * work->x[j];
    // 加c到第一行，并抵消P静态正则化的影响
  }

  if (data->p > 0){
    SpMtv(data->A, work->t, work->xbuff); // A^t t
    axpy(work->kkt->psi, work->xbuff, work->kkt->psi, 1.0, data->n); // 加到第一行
    SpMv(data->A, work->x, &work->kkt->psi[data->n]); // Ax 加到第二行 
    axpy(data->b, &work->kkt->psi[data->n], &work->kkt->psi[data->n], -1.0, data->p); //加-b到第二行
  }
 
  if (data->m > 0){
    SpMtv(data->G, work->s, work->xbuff); //G^t s
    axpy(work->kkt->psi, work->xbuff, work->kkt->psi, 1.0, data->n); // 加到第一行
    SpMv(data->G, work->x, &work->kkt->psi[data->n + data->p + data->m]);  // 加到第四行
    axpy(work->xi, work->z, &work->kkt->psi[data->n + data->p], -1.0, data->m); //加z-ξ到第三行  ////////////// 注意这里用到的z需要在compute_zy中预先计算
    axpy(work->xi, &work->kkt->psi[data->n + data->p + data->m], 
      &work->kkt->psi[data->n + data->p + data->m], 1.0, data->m); // 加ξ到第四行
    axpy(data->h, &work->kkt->psi[data->n + data->p + data->m], 
      &work->kkt->psi[data->n + data->p + data->m], -1.0, data->m); // 加-h到第四行
  }

  work->kkt->phi = 0;
  for (int j = 0;j < data->n + data->p + data->m + data->m;j++){
    work->kkt->phi += work->kkt->psi[j] * work->kkt->psi[j]; 
  }
  work->kkt->phi /= 2;
}
void iprm_compute_obj(IPRMSolver* solver){
  IPRMWorkspace* work = solver->work;
  // Compute objective.
  QOCOFloat obj = dot(work->x, work->data->c, work->data->n);
  USpMv(work->data->P, work->x, work->xbuff);

  // Correct for regularization in P.
  QOCOFloat regularization_correction = 0.0;
  for (QOCOInt i = 0; i < work->data->n; ++i) {
    regularization_correction +=
        solver->settings->kkt_static_reg * work->x[i] * work->x[i];
  }
  obj += 0.5 *
         (dot(work->xbuff, work->x, work->data->n) - regularization_correction);

  solver->sol->obj = obj;
}
void iprm_eval_adjust(IPRMSolver* solver){

  IPRMWorkspace* work = solver->work;
  iprm_compute_zy(work);
  iprm_compute_psi_phi(solver);
  while (work->mu > solver->settings->eta * work->kkt->phi && work->mu > solver->settings->epsilon){
    work->mu /= solver->settings->eta;
    iprm_compute_zy(work);
    iprm_compute_psi_phi(solver);
  }
  work->gamma = qoco_min(solver->settings->gamma0, work->mu / work->kkt->phi);
  iprm_compute_obj(solver);
}


void iprm_linesearch_update(IPRMSolver* solver){
  // 线搜索的同时更新x, mu, t, s, xi
  // 由于嵌套了函数调用，避免使用work->xbuff
  IPRMWorkspace* work = solver->work;
  IPRMProblemData* data = work->data;
  IPRMKKT* kkt = work->kkt;
  IPRMSettings* settings = solver->settings;

  QOCOFloat* Dx = kkt->xt;
  QOCOFloat* Dt = &kkt->xt[data->n];
  QOCOFloat* Ds = work->Ds;
  QOCOFloat* Dxi = work->Dxi;
  
  QOCOFloat* x0 = kkt->xtbuff1;
  QOCOFloat* t0 = &kkt->xtbuff1[data->n];

  QOCOFloat* s0 = work->sbuff;
  QOCOFloat* xi0 = work->xibuff; 
  QOCOFloat mu0 = work->mu;
  for (int i = 0;i < data->n;i++){
    x0[i] = work->x[i];
  }
  for (int i = 0;i < data->p;i++){
    t0[i] = work->t[i];
  }
  for (int i = 0;i < data->m;i++){
    s0[i] = work->s[i];
  }
  for (int i = 0;i < data->m;i++){
    xi0[i] = work->xi[i];
  }
  QOCOFloat a = 1.0;
  QOCOFloat phi0 = work->kkt->phi; // 更新前的phi值
  for(QOCOInt i = 0;i < settings->linesearch_iters;i++){
    work->mu = mu0 + a * work->Dmu;
    for (int j = 0;j < work->data->n;j++){
      work->x[j] = x0[j] + a * Dx[j];
    }
    for (int j = 0;j < work->data->p;j++){
      work->t[j] = t0[j] + a * Dt[j];
    }
    for (int j = 0;j < work->data->m;j++){
      work->s[j] = s0[j] + a * Ds[j];
    }
    for (int j = 0;j < work->data->m;j++){
      work->xi[j] = xi0[j] + a * Dxi[j];
    }
    iprm_compute_zy(work);
    iprm_compute_psi_phi(solver);
    iprm_compute_obj(solver);
    if (work->kkt->phi <= (1 - 2 * settings->tau * a) * phi0){
	    printf("alpha = %lf\n", a);
      return; // 达到更新标准
    }

    a *= settings->delta;

  }
  
  // 没有找到，重启步长
  // 对比一下结果
  
  a = 1;
  printf("alpha = %lf\n", a);
  //work->mu = mu0 + a * work->Dmu; // mu不更新，只更新其他的
  for (int j = 0;j < work->data->n;j++){
    work->x[j] = x0[j] + a * Dx[j];
  }
  for (int j = 0;j < work->data->p;j++){
    work->t[j] = t0[j] + a * Dt[j];
  }
  for (int j = 0;j < work->data->m;j++){
    work->s[j] = s0[j] + a * Ds[j];
  }
  for (int j = 0;j < work->data->m;j++){
    work->xi[j] = xi0[j] + a * Dxi[j];
  }
  iprm_compute_zy(work);
  iprm_compute_psi_phi(solver);
  iprm_compute_obj(solver);
  
}