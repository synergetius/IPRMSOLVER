
#ifndef CONE_H
#define CONE_H

#include "linalg.h"

void iprm_compute_zy(IPRMWorkspace* work);
void iprm_compute_obj(IPRMSolver* solver);
void iprm_eval_adjust(IPRMSolver* solver);
void iprm_compute_psi_phi(IPRMSolver* solver);
IPRMFloat soc_det(const IPRMFloat* u, IPRMInt n);
void soc_product(const IPRMFloat* u, const IPRMFloat* v, IPRMFloat* p,
                 IPRMInt n);
void cone_product(const IPRMFloat* u, const IPRMFloat* v, IPRMFloat* p,
                  IPRMInt l, IPRMInt nsoc, const IPRMInt* q);
void soc_sqrt(const IPRMFloat* u, IPRMFloat* p, IPRMInt n);
void cone_sqrt(const IPRMFloat* u, IPRMFloat* p,
                  IPRMInt l, IPRMInt nsoc, const IPRMInt* q);
void iprm_linesearch_update(IPRMSolver* solver);

#endif