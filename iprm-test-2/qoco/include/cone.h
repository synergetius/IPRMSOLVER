
#ifndef CONE_H
#define CONE_H

#include "linalg.h"

void iprm_compute_zy(IPRMWorkspace* work);
void iprm_compute_obj(IPRMSolver* solver);
void iprm_eval_adjust(IPRMSolver* solver);
void iprm_compute_psi_phi(IPRMSolver* solver);
void iprm_linesearch_update_fast(IPRMSolver* solver);
void iprm_linesearch_update(IPRMSolver* solver);

#endif