

#ifndef QOCO_API_H
#define QOCO_API_H

#include "definitions.h"
#include "enums.h"
#include "input_validation.h"
#include "kkt.h"
#include "linalg.h"
#include "qoco_error.h"
#include "structs.h"
#include "utils.h"
QOCOInt iprm_setup(IPRMSolver* solver, QOCOInt n, QOCOInt m, QOCOInt p,
    QOCOCscMatrix* P, QOCOFloat* c, QOCOCscMatrix* A,
    QOCOFloat* b, QOCOCscMatrix* G, QOCOFloat* h,
    IPRMSettings* settings);
void qoco_set_csc(QOCOCscMatrix* A, QOCOInt m, QOCOInt n, QOCOInt Annz,
                  QOCOFloat* Ax, QOCOInt* Ap, QOCOInt* Ai);


void iprm_set_default_settings(IPRMSettings* settings);
QOCOInt iprm_solve(IPRMSolver* solver);
QOCOInt iprm_cleanup(IPRMSolver* solver);
#endif /* #ifndef QOCO_API_H */