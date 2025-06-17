/**
 * @file kkt.h
 * @author Govind M. Chari <govindchari1@gmail.com>
 *
 * @section LICENSE
 *
 * Copyright (c) 2024, Govind M. Chari
 * This source code is licensed under the BSD 3-Clause License
 *
 * @section DESCRIPTION
 *
 * Provides various functions for solving, constructing and updating KKT
 * systems.
 */

#ifndef KKT_H
#define KKT_H

#include "cone.h"
#include "linalg.h"
#include "qdldl.h"
#include "structs.h"

void iprm_allocate_kkt(IPRMWorkspace* work);
void iprm_construct_kkt(IPRMSolver* solver);
void iprm_initialize(IPRMSolver* solver);
void iprm_update_blocks(IPRMSolver* solver);
void iprm_construct_kkt_rhs(IPRMWorkspace* work);
void iprm_kkt_solve(IPRMSolver* solver);


#endif /* #ifndef KKT_H */