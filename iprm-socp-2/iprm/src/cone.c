

#include "cone.h"
#include "utils.h"
#include <math.h>

IPRMFloat soc_det(const IPRMFloat* u, IPRMInt n){
  IPRMFloat det = u[0] * u[0];
  for (int i = 1;i < n;i++){
    det -= u[i] * u[i];
  }
  return det;
}
void soc_product(const IPRMFloat* u, const IPRMFloat* v, IPRMFloat* p,
                 IPRMInt n)
{
  // u和v不可与p相同
  p[0] = dot(u, v, n);
  for (IPRMInt i = 1; i < n; ++i) {
    p[i] = u[0] * v[i] + v[0] * u[i];
  }
}
void cone_product(const IPRMFloat* u, const IPRMFloat* v, IPRMFloat* p,
                  IPRMInt l, IPRMInt nsoc, const IPRMInt* q)
{
  IPRMInt idx;
  // Compute LP cone product.
  for (idx = 0; idx < l; ++idx) {
    p[idx] = u[idx] * v[idx];
  }

  // Compute second-order cone product.
  for (IPRMInt i = 0; i < nsoc; ++i) {
    soc_product(&u[idx], &v[idx], &p[idx], q[i]);
    idx += q[i];
  }
}
void soc_sqrt(const IPRMFloat* u, IPRMFloat* p, IPRMInt n){
  // u可与p相同
  IPRMFloat norm = norm_2(&u[1], n - 1);
  IPRMFloat a = iprm_sqrt(u[0] + norm);
  IPRMFloat b = iprm_sqrt(u[0] - norm);
  p[0] = (a + b) / 2;
  if (norm < 1e-15){
    // 防止除零
    //printf("too small.\n");
    if (n > 1){
      p[1] = (a - b) / 2;
      for (int i = 2;i < n;i++){
        p[i] = 0;
      }
    }
  }
  else{
    for (int i = 1;i < n;i++){
      p[i] = u[i] * (a - b) / (2 * norm);
    }
  }
}
void cone_sqrt(const IPRMFloat* u, IPRMFloat* p,
                  IPRMInt l, IPRMInt nsoc, const IPRMInt* q)
{
  IPRMInt idx;
  for (idx = 0; idx < l; ++idx) {
    p[idx] = iprm_sqrt(u[idx]);
  }
  for (IPRMInt i = 0; i < nsoc; ++i) {
    soc_sqrt(&u[idx], &p[idx], q[i]);
    idx += q[i];
  }
}
void iprm_compute_zy(IPRMWorkspace* work){
  axpy(work->xi, work->s, work->v, -work->rho, work->data->m);
  cone_product(work->v, work->v, work->u, work->data->l, work->data->nsoc, work->data->q);
  IPRMInt idx;
  for (idx = 0;idx < work->data->l;idx++){
    work->u[idx] += 4 * work->rho * work->mu;
  }
  for (IPRMInt c = 0;c < work->data->nsoc;c++){
    work->u[idx] += 4 * work->rho * work->mu;
    idx += work->data->q[c];
  }
  cone_sqrt(work->u, work->u, work->data->l, work->data->nsoc, work->data->q);
  for (IPRMInt i = 0; i < work->data->m; ++i) {
    work->z[i] = (work->u[i] - work->v[i]) / (2 * work->rho);
    work->y[i] = (work->u[i] + work->v[i]) / (2 * work->rho);
  }
}
void iprm_compute_psi_phi(IPRMSolver* solver)
{
  // 计算psi和phi（原始的KKT方程残差）
  IPRMWorkspace* work = solver->work;
  IPRMProblemData* data = solver->work->data;
  //printf("#1\n");
  USpMv(data->P, work->x, work->kkt->psi);
  for (int j = 0; j < data->n; ++j) {
    work->kkt->psi[j] += data->c[j] - solver->settings->kkt_static_reg * work->x[j];
    // 加c到第一行，并抵消P静态正则化的影响
  }
  //printf("#2\n");
  if (data->p > 0){
    SpMtv(data->A, work->t, work->xbuff); // A^t t
    axpy(work->kkt->psi, work->xbuff, work->kkt->psi, 1.0, data->n); // 加到第一行
    SpMv(data->A, work->x, &work->kkt->psi[data->n]); // Ax 加到第二行 
    axpy(data->b, &work->kkt->psi[data->n], &work->kkt->psi[data->n], -1.0, data->p); //加-b到第二行
  }
  //printf("#3\n");
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
  //printf("#4\n");
  
  // 采用IPRSOCP中的实现
  work->kkt->phi = 0;
  for (int j = data->n + data->p;j < data->n + data->p + data->m;j++){
    work->kkt->phi += work->kkt->psi[j] * work->kkt->psi[j]; 
  }
  work->kkt->phi /= 2;
  
  work->kkt->square_psi = 0;
  for (int j = 0;j < data->n + data->p + data->m + data->m;j++){
    work->kkt->square_psi += work->kkt->psi[j] * work->kkt->psi[j]; 
  }
  work->kkt->square_psi /= 2;
  
  //printf("finished psi phi\n");
}

void iprm_compute_obj(IPRMSolver* solver){
  IPRMWorkspace* work = solver->work;
  // Compute objective.
  IPRMFloat obj = dot(work->x, work->data->c, work->data->n);
  USpMv(work->data->P, work->x, work->xbuff);

  // Correct for regularization in P.
  IPRMFloat regularization_correction = 0.0;
  for (IPRMInt i = 0; i < work->data->n; ++i) {
    regularization_correction +=
        solver->settings->kkt_static_reg * work->x[i] * work->x[i];
  }
  obj += 0.5 *
         (dot(work->xbuff, work->x, work->data->n) - regularization_correction);

  solver->sol->obj = obj;
}
void iprm_eval_adjust(IPRMSolver* solver){

  IPRMWorkspace* work = solver->work;
  IPRMProblemData* data = solver->work->data;
  iprm_compute_zy(work);
  
  iprm_compute_psi_phi(solver);
  ////////
  while (work->mu > solver->settings->eta * work->kkt->phi && work->mu > solver->settings->epsilon){
    work->mu /= solver->settings->eta;
    iprm_compute_zy(work);
    iprm_compute_psi_phi(solver);
  }
  work->gamma = iprm_min(solver->settings->gamma0, work->mu / work->kkt->phi);
  iprm_compute_obj(solver);
}
void iprm_linesearch_update(IPRMSolver* solver){
  // 线搜索的同时更新x, mu, t, s, xi
  IPRMWorkspace* work = solver->work;
  IPRMProblemData* data = work->data;
  IPRMKKT* kkt = work->kkt;
  IPRMSettings* settings = solver->settings;
  IPRMFloat* Dx = kkt->w;
  IPRMFloat* Dt = &kkt->w[data->n + data->nsoc + data->nsoc];
  IPRMFloat* Ds = &kkt->w[data->n + data->nsoc + data->nsoc + data->p];
  IPRMFloat* x0 = kkt->kktbuff1;
  IPRMFloat* t0 = &kkt->kktbuff1[data->n];
  IPRMFloat* s0 = &kkt->kktbuff1[data->n + data->nsoc + data->nsoc + data->p];
  IPRMFloat* xi0 = &kkt->kktbuff2[data->n + data->nsoc + data->nsoc + data->p]; 
  IPRMFloat mu0 = work->mu;
  copy_arrayf(work->x, x0, data->n);
  copy_arrayf(work->t, t0, data->p);
  copy_arrayf(work->s, s0, data->m);
  copy_arrayf(work->xi, xi0, data->m);
  IPRMFloat a = 1.0;
  IPRMFloat phi0 = work->kkt->phi; // 更新前的phi值
  for(IPRMInt i = 0;i < settings->linesearch_iters;i++){
    work->mu = mu0 + a * work->Dmu;
    axpy(Dx, x0, work->x, a, data->n);
    axpy(Dt, t0, work->t, a, data->p);
    axpy(Ds, s0, work->s, a, data->m);
    axpy(work->Dxi, xi0, work->xi, a, data->m);
	  iprm_compute_zy(work);
    iprm_compute_psi_phi(solver);
    //printf("%.12e %.12e\n", work->kkt->phi, phi0);
    if (work->kkt->phi <= (1 - 2 * settings->tau * a) * phi0){
	    printf("alpha = %lf\n", a);
	    iprm_compute_obj(solver);
      return; // 达到更新标准
    }
    a *= settings->delta;
  }
  
  // 没有找到，重启步长
  
  a = 1;
  printf("alpha = %lf\n", a);
  //work->mu = mu0 + a * work->Dmu; // mu不更新，只更新其他的
  axpy(Dx, x0, work->x, a, data->n);
  axpy(Dt, t0, work->t, a, data->p);
  axpy(Ds, s0, work->s, a, data->m);
  axpy(work->Dxi, xi0, work->xi, a, data->m);
  iprm_compute_zy(work);
  iprm_compute_psi_phi(solver);
  iprm_compute_obj(solver);
  
}