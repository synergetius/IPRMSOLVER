#ifndef KKT_H
#define KKT_H

#include "cone.h"
#include "linalg.h"
#include "qdldl.h"
#include "structs.h"
void print_K(IPRMCscMatrix* K);
void iprm_allocate_kkt(IPRMWorkspace* work);
void iprm_construct_kkt(IPRMSolver* solver);
void iprm_initialize(IPRMSolver* solver);
void iprm_update_blocks(IPRMSolver* solver);
void iprm_construct_kkt_rhs(IPRMWorkspace* work);
void iprm_kkt_multiply(IPRMSolver* solver, IPRMFloat* src, IPRMFloat* tar);
void iprm_kkt_solve(IPRMSolver* solver);


#endif /* #ifndef KKT_H */