

#ifndef QOCO_API_H
#define QOCO_API_H

#include "definitions.h"
#include "enums.h"
#include "input_validation.h"
#include "kkt.h"
#include "linalg.h"
#include "iprm_error.h"
#include "structs.h"
#include "utils.h"
IPRMInt iprm_setup(IPRMSolver* solver, IPRMInt n, IPRMInt m, IPRMInt p,
  IPRMCscMatrix* P, IPRMFloat* c, IPRMCscMatrix* A,
  IPRMFloat* b, IPRMCscMatrix* G, IPRMFloat* h, IPRMInt l,
  IPRMInt nsoc, IPRMInt* q, IPRMSettings* settings);
void iprm_set_csc(IPRMCscMatrix* A, IPRMInt m, IPRMInt n, IPRMInt Annz,
                  IPRMFloat* Ax, IPRMInt* Ap, IPRMInt* Ai);


void iprm_set_default_settings(IPRMSettings* settings);
IPRMInt iprm_solve(IPRMSolver* solver);
IPRMInt iprm_cleanup(IPRMSolver* solver);
#endif /* #ifndef QOCO_API_H */