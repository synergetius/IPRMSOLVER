#ifndef STRUCTS_H
#define STRUCTS_H

#include "definitions.h"
#include "timer.h"


typedef struct {
  /** Number of rows. */
  IPRMInt m;
  /** Number of columns. */
  IPRMInt n;
  /** Number of nonzero elements. */
  IPRMInt nnz;
  /** Row indices (length: nnz). */
  IPRMInt* i;
  /** Column pointers (length: n+1). */
  IPRMInt* p;
  /** Data (length: nnz). */
  IPRMFloat* x;
} IPRMCscMatrix;
typedef struct{
  IPRMCscMatrix* P; 
  IPRMFloat* c;
  IPRMCscMatrix* A;
  IPRMCscMatrix* At;
  IPRMFloat* b;
  IPRMCscMatrix* G;
  IPRMCscMatrix* Gt;
  IPRMFloat* h;
  IPRMInt n;
  IPRMInt m;
  IPRMInt p;
} IPRMProblemData;
typedef struct{
  IPRMInt max_iters;
  IPRMInt iter_ref_iters;
  IPRMInt linesearch_iters;
  IPRMFloat kkt_static_reg;
  IPRMFloat kkt_dynamic_reg;
  
  IPRMFloat mu0;
  IPRMFloat rho0;
  IPRMFloat eta;
  IPRMFloat gamma0;
  IPRMFloat delta;
  IPRMFloat tau;
  IPRMFloat sigma;
  IPRMFloat epsilon;
  unsigned char verbose;
} IPRMSettings;
typedef struct{
  IPRMCscMatrix* K;
  IPRMInt* Pnzadded_idx;
  IPRMInt Pnum_nzadded;
  IPRMInt* PregtoKKT;
  IPRMInt* AtoKKT;
  IPRMInt* GtoKKT;
  IPRMInt* ZtoKKT; // 右下角对角阵
  /** Permutation vector. */
  IPRMInt* p;
  /** Inverse of permutation vector. */
  IPRMInt* pinv;
  /** Elimination tree for LDL factorization of K. */
  IPRMInt* etree;
  IPRMInt* Lnz;
  IPRMFloat* Lx;
  IPRMInt* Lp;
  IPRMInt* Li;
  IPRMFloat* D;
  IPRMFloat* Dinv;
  IPRMInt* iwork;
  unsigned char* bwork;
  IPRMFloat* fwork;

  IPRMFloat* psi;
  IPRMFloat phi;

  IPRMFloat* xts; // x, t, s memory
  IPRMFloat* rhs;
  IPRMFloat* npm_buff; // buffer (length = n + p + m)
  IPRMFloat* npm_buff2;
  IPRMFloat* psi_buff;
  IPRMFloat* psi_buff2;
} IPRMKKT;
typedef struct{
  IPRMProblemData* data;
  IPRMTimer solve_timer;
  IPRMKKT* kkt;
  IPRMFloat* x;
  IPRMFloat* t;
  IPRMFloat* s;
  IPRMFloat* xi;
  IPRMFloat mu;
  IPRMFloat rho;
  IPRMFloat* z;
  IPRMFloat* y;
  IPRMFloat* xbuff;
  IPRMFloat* tbuff;
  IPRMFloat* ubuff1;
  IPRMFloat gamma;
  IPRMFloat Dmu; // delta mu
  IPRMFloat* Dxi;
} IPRMWorkspace;
typedef struct{
  IPRMFloat* x;
  IPRMFloat* xi;
  IPRMFloat* t;
  IPRMFloat* s;
  IPRMInt iters;
  IPRMFloat setup_time_sec;
  IPRMFloat solve_time_sec;
  IPRMFloat obj;
  IPRMFloat phi;
  IPRMFloat pres;
  IPRMFloat dres;
  IPRMFloat gap;
  IPRMInt status;
  
} IPRMSolution;
typedef struct{
  IPRMSettings* settings;
  IPRMWorkspace* work;
  IPRMSolution* sol;
} IPRMSolver;
#endif 