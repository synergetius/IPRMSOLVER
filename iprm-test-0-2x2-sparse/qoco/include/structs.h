#ifndef STRUCTS_H
#define STRUCTS_H

#include "definitions.h"
#include "timer.h"


typedef struct {
  /** Number of rows. */
  QOCOInt m;
  /** Number of columns. */
  QOCOInt n;
  /** Number of nonzero elements. */
  QOCOInt nnz;
  /** Row indices (length: nnz). */
  QOCOInt* i;
  /** Column pointers (length: n+1). */
  QOCOInt* p;
  /** Data (length: nnz). */
  QOCOFloat* x;
} QOCOCscMatrix;
typedef struct{
  QOCOCscMatrix* M;
  QOCOCscMatrix* MtM;
  QOCOInt** ie;
  QOCOInt** je;
  QOCOInt* cnt;
} MtMPipe;
typedef struct{
  QOCOCscMatrix* A;
  QOCOCscMatrix* B;
  QOCOCscMatrix* ApkB;
  QOCOInt* AtoApkB;
  QOCOInt* BtoApkB;
} ApkBPipe;
typedef struct{
  QOCOCscMatrix* P; 
  QOCOFloat* c;
  QOCOCscMatrix* A;
  QOCOCscMatrix* At;
  QOCOFloat* b;
  QOCOCscMatrix* G;
  QOCOCscMatrix* Gt;
  QOCOFloat* h;
  QOCOInt n;
  QOCOInt m;
  QOCOInt p;
} IPRMProblemData;
typedef struct{
  QOCOInt max_iters;
  QOCOInt iter_ref_iters;
  QOCOInt linesearch_iters;
  QOCOFloat kkt_static_reg;
  QOCOFloat kkt_dynamic_reg;
  
  QOCOFloat mu0;
  QOCOFloat rho0;
  QOCOFloat eta;
  QOCOFloat gamma0;
  QOCOFloat delta;
  QOCOFloat tau;
  QOCOFloat sigma;
  QOCOFloat epsilon;
  unsigned char verbose;
} IPRMSettings;
typedef struct{
  QOCOCscMatrix* K;
  QOCOInt* Pnzadded_idx;
  QOCOInt Pnum_nzadded;
  QOCOInt* AtoKKT;
  /** Permutation vector. */
  QOCOInt* p;
  /** Inverse of permutation vector. */
  QOCOInt* pinv;
  /** Elimination tree for LDL factorization of K. */
  QOCOInt* etree;
  QOCOInt* Lnz;
  QOCOFloat* Lx;
  QOCOInt* Lp;
  QOCOInt* Li;
  QOCOFloat* D;
  QOCOFloat* Dinv;
  QOCOInt* iwork;
  unsigned char* bwork;
  QOCOFloat* fwork;

  QOCOFloat* psi;
  QOCOFloat phi;

  QOCOFloat* rhs;

  ////////////// 2x2版本添加的变量
  //QOCOCscMatrix* GTMG;
  //QOCOInt* GTMGtoKKT;
  QOCOFloat* xt;
  QOCOFloat* xtbuff1;
  QOCOFloat* xtbuff2;
  // 2x2sparse添加的变量
  //QOCOInt* MtMtoKKT;
  QOCOInt* PpkMtMtoKKT;
} IPRMKKT;

typedef struct{
  IPRMProblemData* data;
  QOCOTimer solve_timer;
  IPRMKKT* kkt;
  QOCOFloat* x;
  QOCOFloat* t;
  QOCOFloat* s;
  QOCOFloat* xi;
  QOCOFloat mu;
  QOCOFloat rho;
  QOCOFloat* z;
  QOCOFloat* y;
  QOCOFloat* xbuff;
  QOCOFloat* tbuff;
  //QOCOFloat* ubuff1;
  QOCOFloat* sbuff;
  QOCOFloat* xibuff;
  QOCOFloat gamma;
  QOCOFloat Dmu; // delta mu
  QOCOFloat* Dxi;
  ////////////// 2x2版本添加的变量
  QOCOFloat* rhoy;
  QOCOFloat* Ds;
  // 2x2-sparse
  MtMPipe* mtmpipe;
  ApkBPipe* apkbpipe;
} IPRMWorkspace;
typedef struct{
  QOCOFloat* x;
  QOCOFloat* xi;
  QOCOFloat* t;
  QOCOFloat* s;
  QOCOInt iters;
  QOCOFloat setup_time_sec;
  QOCOFloat solve_time_sec;
  QOCOFloat obj;
  QOCOFloat phi;
  QOCOFloat pres;
  QOCOFloat dres;
  QOCOFloat gap;
  QOCOInt status;
  
} IPRMSolution;
typedef struct{
  IPRMSettings* settings;
  IPRMWorkspace* work;
  IPRMSolution* sol;
} IPRMSolver;
#endif 