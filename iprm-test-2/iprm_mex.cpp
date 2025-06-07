/* Copyright (c) 2025, Govind M. Chari <govindchari1@gmail.com> */
/* Copyright (c) 2017, Bartolomeo Stellato */
/* This source code is licensed under the BSD 3-Clause License  */

#include "iprm_mex.hpp"
#include <iostream>

const char *IPRM_SETTINGS_FIELDS[] = {
  "max_iters",
  "iter_ref_iters",
  "linesearch_iters",
  "kkt_static_reg",
  "kkt_dynamic_reg",
  "mu0",
  "rho0",
  "eta",
  "gamma0",
  "delta",
  "tau",
  "sigma",
  "epsilon",
  "verbose"
};
const char *IPRM_RESULT_FIELDS[] = {
  "x",
  "xi",
  "t",
  "s",
  "iters",
  "setup_time_sec",
  "solve_time_sec",
  "obj",
  "phi",
  "status"
  
};


const char *CSC_FIELDS[] = {"m",
                            "n",
                            "nnz",
                            "i",
                            "p",
                            "x"};

// wrapper class for all qoco data and settings
class IprmData
{
public:
  IprmData()
  {
    solver = NULL;
  }
  IPRMSolver *solver; // Solver
};

// internal utility functions
void castToDoubleArr(QOCOFloat *arr, double *arr_out, QOCOInt len);
void setToNaN(double *arr_out, QOCOInt len);
void copyMxStructToSettings(const mxArray *, IPRMSettings *);
void copyUpdatedSettingsToWork(const mxArray *, IPRMSettings *);
void castCintToDoubleArr(QOCOInt *arr, double *arr_out, QOCOInt len);
QOCOInt *copyToCintVector(mwIndex *vecData, QOCOInt numel);
QOCOInt *copyDoubleToCintVector(double *vecData, QOCOInt numel);
QOCOInt *roundDoubleIntoInt(double *vecData, QOCOInt numel);
QOCOFloat *copyToCfloatVector(double *vecData, QOCOInt numel);
// mxArray *copyInfoToMxStruct(OSQPInfo *info);
mxArray *copySettingsToMxStruct(IPRMSettings *settings);
// mxArray *copyCscMatrixToMxStruct(csc *M);
// mxArray *copyDataToMxStruct(OSQPWorkspace *work);
// mxArray *copyLinsysSolverToMxStruct(OSQPWorkspace *work);
// mxArray *copyScalingToMxStruct(OSQPWorkspace *work);
// mxArray *copyWorkToMxStruct(OSQPWorkspace *work);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  IprmData *iprmData; 
  
  // Exitflag
  QOCOInt exitflag = 1;

  // Static string for static methods
  char stat_string[64];

  // Get the command string
  char cmd[64];
  if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
    mexErrMsgTxt("First input should be a command string less than 64 characters long.");

  // new object
  if (!strcmp("new", cmd))
  {
    // Check parameters
    if (nlhs != 1)
    {
      mexErrMsgTxt("New: One output expected.");
    }
    // Return a handle to a new C++ wrapper instance
    plhs[0] = convertPtr2Mat<IprmData>(new IprmData);
    return;
  }

  // Check for a second input, which should be the class instance handle or string 'static'
  if (nrhs < 2)
    mexErrMsgTxt("Second input should be a class instance handle or the string 'static'.");

  if (mxGetString(prhs[1], stat_string, sizeof(stat_string)))
  {
    // If we are dealing with non-static methods, get the class instance pointer from the second input
    iprmData = convertMat2Ptr<IprmData>(prhs[1]);
  }
  else
  {
    if (strcmp("static", stat_string))
    {
      mexErrMsgTxt("Second argument for static functions is string 'static'");
    }
  }

  // delete the object and its data
  if (!strcmp("delete", cmd))
  {

    // clean up the problem workspace
    if (iprmData->solver)
    {
      iprm_cleanup(iprmData->solver);
    }

    // clean up the handle object
    destroyObject<IprmData>(prhs[1]);
    // Warn if other commands were ignored
    if (nlhs != 0 || nrhs != 2)
      mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
    return;
  }

  // report the current settings
  if (!strcmp("current_settings", cmd))
  {

    // throw an error if this is called before solver is configured
    if (!iprmData->solver)
    {
      mexErrMsgTxt("Solver is uninitialized.  No settings have been configured.");
    }
    // report the current settings
    plhs[0] = copySettingsToMxStruct(iprmData->solver->settings);
    return;
  }

  // write_settings
  if (!strcmp("update_settings", cmd))
  {

    // overwrite the current settings.  Mex function is responsible
    // for disallowing overwrite of selected settings after initialization,
    // and for all error checking
    // throw an error if this is called before solver is configured
    if (!iprmData->solver)
    {
      mexErrMsgTxt("Solver is uninitialized.  No settings have been configured.");
    }

    copyUpdatedSettingsToWork(prhs[2], iprmData->solver->settings);
    return;
  }

  // report the default settings
  if (!strcmp("default_settings", cmd))
  {

    // Warn if other commands were ignored
    if (nrhs > 2)
      mexWarnMsgTxt("Default settings: unexpected number of arguments.");

    // Create a Settings structure in default form and report the results
    // Useful for external solver packages (e.g. Yalmip) that want to
    // know which solver settings are supported
    IPRMSettings *defaults = (IPRMSettings *)mxCalloc(1, sizeof(IPRMSettings));
    iprm_set_default_settings(defaults);
    plhs[0] = copySettingsToMxStruct(defaults);
    mxFree(defaults);
    return;
  }

  // setup
  if (!strcmp("setup", cmd))
  {

    // throw an error if this is called more than once
    if (iprmData->solver)
    {
      mexErrMsgTxt("Solver is already initialized with problem data.");
    }

    // Create data and settings containers
    IPRMSettings *settings = new IPRMSettings;

    // handle the problem data first.  Matlab-side
    // class wrapper is responsible for ensuring that
    // P and A are sparse matrices,  everything
    // else is a dense vector and all inputs are
    // of compatible dimension

    const QOCOInt n = (QOCOInt)mxGetScalar(prhs[2]);
    const QOCOInt m = (QOCOInt)mxGetScalar(prhs[3]);
    const QOCOInt p = (QOCOInt)mxGetScalar(prhs[4]);
    const mxArray *P = prhs[5];
    const mxArray *c = prhs[6];
    const mxArray *A = prhs[7];
    const mxArray *b = prhs[8];
    const mxArray *G = prhs[9];
    const mxArray *h = prhs[10];


    // Create Data Structure
    QOCOFloat *cvec = copyToCfloatVector(mxGetPr(c), n);
    QOCOFloat *bvec = copyToCfloatVector(mxGetPr(b), p);
    QOCOFloat *hvec = copyToCfloatVector(mxGetPr(h), m);


    // Matrix P:  nnz = P->p[n]
    QOCOCscMatrix *Pcsc = new QOCOCscMatrix;
    QOCOInt *Pp = copyToCintVector(mxGetJc(P), n + 1);
    QOCOInt *Pi = copyToCintVector(mxGetIr(P), Pp[n]);
    QOCOFloat *Px = copyToCfloatVector(mxGetPr(P), Pp[n]);
    qoco_set_csc(Pcsc, n, n, Pp[n], Px, Pp, Pi);

    // // Matrix A: nnz = A->p[n]
    QOCOCscMatrix *Acsc = new QOCOCscMatrix;
    QOCOInt *Ap = copyToCintVector(mxGetJc(A), n + 1);
    QOCOInt *Ai = copyToCintVector(mxGetIr(A), Ap[n]);
    QOCOFloat *Ax = copyToCfloatVector(mxGetPr(A), Ap[n]);
    qoco_set_csc(Acsc, p, n, Ap[n], Ax, Ap, Ai);

    // Matrix G: nnz = G->p[n]
    QOCOCscMatrix *Gcsc = new QOCOCscMatrix;
    QOCOInt *Gp = copyToCintVector(mxGetJc(G), n + 1);
    QOCOInt *Gi = copyToCintVector(mxGetIr(G), Gp[n]);
    QOCOFloat *Gx = copyToCfloatVector(mxGetPr(G), Gp[n]);
    qoco_set_csc(Gcsc, m, n, Gp[n], Gx, Gp, Gi);

    // Create Settings
    const mxArray *mxSettings = prhs[11]; ///////////
    if (mxIsEmpty(mxSettings))
    {
      // use defaults
      iprm_set_default_settings(settings);
    }
    else
    {
      // populate settings structure from mxArray input
      copyMxStructToSettings(mxSettings, settings);
    }

    // Setup solver.
    iprmData->solver = new IPRMSolver;
    exitflag = iprm_setup(iprmData->solver, n, m, p, Pcsc, cvec, Acsc, bvec, Gcsc, hvec, settings);

    // cleanup temporary structures
    if (cvec)
      qoco_free(cvec);
    if (bvec)
      qoco_free(bvec);
    if (hvec)
      qoco_free(hvec);
    // if (qvec)
    //   qoco_free(qvec);
    if (Px)
      qoco_free(Px);
    if (Pi)
      qoco_free(Pi);
    if (Pp)
      qoco_free(Pp);
    if (Pcsc)
      delete Pcsc;
    if (Ax)
      qoco_free(Ax);
    if (Ai)
      qoco_free(Ai);
    if (Ap)
      qoco_free(Ap);
    if (Acsc)
      delete Acsc;
    if (Gx)
      qoco_free(Gx);
    if (Gi)
      qoco_free(Gi);
    if (Gp)
      qoco_free(Gp);
    if (Gcsc)
      delete Gcsc;

    // Settings
    if (settings)
      delete settings;

    // Report error (if any)
    if (exitflag != 0)
    {
      mexErrMsgTxt("Invalid problem setup");
    }
    return;
  }

  if (!strcmp("get_dimensions", cmd))
  {

    // throw an error if this is called before solver is configured
    if (!iprmData->solver)
    {
      mexErrMsgTxt("Solver has not been initialized.");
    }

    plhs[0] = mxCreateDoubleScalar(iprmData->solver->work->data->m);
    plhs[1] = mxCreateDoubleScalar(iprmData->solver->work->data->n);
    plhs[2] = mxCreateDoubleScalar(iprmData->solver->work->data->p);

    return;
  }

  // solve
  if (!strcmp("solve", cmd))
  {
    if (nlhs != 10 || nrhs != 2) ////////////////////
    {
      mexErrMsgTxt("Solve : wrong number of inputs / outputs");
    }
    if (!iprmData->solver)
    {
      mexErrMsgTxt("No problem data has been given.");
    }

    // solve the problem
    iprm_solve(iprmData->solver);

    // Allocate space for solution
    // primal variables
    plhs[0] = mxCreateDoubleMatrix(iprmData->solver->work->data->n, 1, mxREAL);
    // slack variables
    plhs[1] = mxCreateDoubleMatrix(iprmData->solver->work->data->m, 1, mxREAL);
    // dual variables for equality constraints
    plhs[2] = mxCreateDoubleMatrix(iprmData->solver->work->data->p, 1, mxREAL);
    // dual variables for inequality constraints
    plhs[3] = mxCreateDoubleMatrix(iprmData->solver->work->data->m, 1, mxREAL);

    // copy results to mxArray outputs
    castToDoubleArr(iprmData->solver->sol->x, mxGetPr(plhs[0]), iprmData->solver->work->data->n);
    castToDoubleArr(iprmData->solver->sol->xi, mxGetPr(plhs[1]), iprmData->solver->work->data->m);
    castToDoubleArr(iprmData->solver->sol->t, mxGetPr(plhs[2]), iprmData->solver->work->data->p);
    castToDoubleArr(iprmData->solver->sol->s, mxGetPr(plhs[3]), iprmData->solver->work->data->m);

    // copy the rest of the information.
    plhs[4] = mxCreateDoubleScalar(iprmData->solver->sol->iters);
    plhs[5] = mxCreateDoubleScalar(iprmData->solver->sol->setup_time_sec);
    plhs[6] = mxCreateDoubleScalar(iprmData->solver->sol->solve_time_sec);
    
    plhs[7] = mxCreateDoubleScalar(iprmData->solver->sol->obj);
    plhs[8] = mxCreateDoubleScalar(iprmData->solver->sol->phi);
    plhs[9] = mxCreateDoubleScalar(iprmData->solver->sol->status);
    
    return;
  }

  if (!strcmp("constant", cmd))
  { // Return solver constants

    char constant[32];
    mxGetString(prhs[2], constant, sizeof(constant));

    if (!strcmp("QOCO_UNSOLVED", constant))
    {
      plhs[0] = mxCreateDoubleScalar(QOCO_UNSOLVED);
      return;
    }

    if (!strcmp("QOCO_SOLVED", constant))
    {
      plhs[0] = mxCreateDoubleScalar(QOCO_SOLVED);
      return;
    }

    if (!strcmp("QOCO_SOLVED_INACCURATE", constant))
    {
      plhs[0] = mxCreateDoubleScalar(QOCO_SOLVED_INACCURATE);
      return;
    }

    if (!strcmp("QOCO_NUMERICAL_ERROR", constant))
    {
      plhs[0] = mxCreateDoubleScalar(QOCO_NUMERICAL_ERROR);
      return;
    }

    if (!strcmp("QOCO_MAX_ITER", constant))
    {
      plhs[0] = mxCreateDoubleScalar(QOCO_MAX_ITER);
      return;
    }

    mexErrMsgTxt("Constant not recognized.");

    return;
  }

  // Got here, so command not recognized
  mexErrMsgTxt("Command not recognized.");
}

QOCOFloat *copyToCfloatVector(double *vecData, QOCOInt numel)
{
  // This memory needs to be freed!
  QOCOFloat *out = (QOCOFloat *)qoco_malloc(numel * sizeof(QOCOFloat));

  // copy data
  for (QOCOInt i = 0; i < numel; i++)
  {
    out[i] = (QOCOFloat)vecData[i];
  }
  return out;
}

QOCOInt *roundDoubleIntoInt(double *vecData, QOCOInt numel)
{
  // This memory needs to be freed!
  QOCOInt *out = (QOCOInt *)qoco_malloc(numel * sizeof(QOCOInt));

  // round data
  for (QOCOInt i = 0; i < numel; i++)
  {
    out[i] = round(vecData[i]);
  }
  return out;
}

QOCOInt *copyToCintVector(mwIndex *vecData, QOCOInt numel)
{
  // This memory needs to be freed!
  QOCOInt *out = (QOCOInt *)qoco_malloc(numel * sizeof(QOCOInt));

  // copy data
  for (QOCOInt i = 0; i < numel; i++)
  {
    out[i] = (QOCOInt)vecData[i];
  }
  return out;
}

QOCOInt *copyDoubleToCintVector(double *vecData, QOCOInt numel)
{
  // This memory needs to be freed!
  QOCOInt *out = (QOCOInt *)qoco_malloc(numel * sizeof(QOCOInt));

  // copy data
  for (QOCOInt i = 0; i < numel; i++)
  {
    out[i] = (QOCOInt)vecData[i];
  }
  return out;
}

void castCintToDoubleArr(QOCOInt *arr, double *arr_out, QOCOInt len)
{
  for (QOCOInt i = 0; i < len; i++)
  {
    arr_out[i] = (double)arr[i];
  }
}

void castToDoubleArr(QOCOFloat *arr, double *arr_out, QOCOInt len)
{
  for (QOCOInt i = 0; i < len; i++)
  {
    arr_out[i] = (double)arr[i];
  }
}

void setToNaN(double *arr_out, QOCOInt len)
{
  QOCOInt i;
  for (i = 0; i < len; i++)
  {
    arr_out[i] = mxGetNaN();
  }
}

mxArray *copySettingsToMxStruct(IPRMSettings *settings)
{
  //////////////////////////////////////
  int nfields = sizeof(IPRM_SETTINGS_FIELDS) / sizeof(IPRM_SETTINGS_FIELDS[0]);
  mxArray *mxPtr = mxCreateStructMatrix(1, 1, nfields, IPRM_SETTINGS_FIELDS);


  mxSetField(mxPtr, 0, "max_iters", mxCreateDoubleScalar(settings->max_iters));
  mxSetField(mxPtr, 0, "iter_ref_iters", mxCreateDoubleScalar(settings->iter_ref_iters));
  mxSetField(mxPtr, 0, "linesearch_iters", mxCreateDoubleScalar(settings->linesearch_iters));
  mxSetField(mxPtr, 0, "kkt_static_reg", mxCreateDoubleScalar(settings->kkt_static_reg));
  mxSetField(mxPtr, 0, "kkt_dynamic_reg", mxCreateDoubleScalar(settings->kkt_dynamic_reg));
  mxSetField(mxPtr, 0, "mu0", mxCreateDoubleScalar(settings->mu0));
  mxSetField(mxPtr, 0, "rho0", mxCreateDoubleScalar(settings->rho0));
  mxSetField(mxPtr, 0, "eta", mxCreateDoubleScalar(settings->eta));
  mxSetField(mxPtr, 0, "gamma0", mxCreateDoubleScalar(settings->gamma0));
  mxSetField(mxPtr, 0, "delta", mxCreateDoubleScalar(settings->delta));
  mxSetField(mxPtr, 0, "tau", mxCreateDoubleScalar(settings->tau));
  mxSetField(mxPtr, 0, "sigma", mxCreateDoubleScalar(settings->sigma));
  mxSetField(mxPtr, 0, "epsilon", mxCreateDoubleScalar(settings->epsilon));
  mxSetField(mxPtr, 0, "verbose", mxCreateDoubleScalar(settings->verbose));

  return mxPtr;
}

// ======================================================================

void copyMxStructToSettings(const mxArray *mxPtr, IPRMSettings *settings)
{

  // this function assumes that only a complete and validated structure
  // will be passed.  matlab mex-side function is responsible for checking
  // structure validity

  settings->max_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "max_iters"));
  settings->iter_ref_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "iter_ref_iters"));
  settings->linesearch_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "linesearch_iters"));
  settings->kkt_static_reg = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "kkt_static_reg"));
  settings->kkt_dynamic_reg = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "kkt_dynamic_reg"));
  settings->mu0 = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "mu0"));
  settings->rho0 = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "rho0"));
  settings->eta = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "eta"));
  settings->gamma0 = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "gamma0"));
  settings->delta = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "delta"));
  settings->tau = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "tau"));
  settings->sigma = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "sigma"));
  settings->epsilon = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "epsilon"));
  settings->verbose = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "verbose"));
}

void copyUpdatedSettingsToWork(const mxArray *mxPtr, IPRMSettings *settings)
{
  // 此处QOCO的变量命名有些奇怪，等完整理解后再重新检查
  settings->max_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "max_iters"));
  settings->iter_ref_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "iter_ref_iters"));
  settings->linesearch_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "linesearch_iters"));
  settings->kkt_static_reg = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "kkt_static_reg"));
  settings->kkt_dynamic_reg = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "kkt_dynamic_reg"));
  settings->mu0 = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "mu0"));
  settings->rho0 = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "rho0"));
  settings->eta = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "eta"));
  settings->gamma0 = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "gamma0"));
  settings->delta = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "delta"));
  settings->tau = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "tau"));
  settings->sigma = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "sigma"));
  settings->epsilon = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "epsilon"));
  settings->verbose = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "verbose"));
  /*
  settings->max_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "max_iters"));
  settings->bisect_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "sigma"));
  settings->ruiz_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "scaling"));
  settings->iter_ref_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "adaptive_rho"));
  settings->kkt_static_reg = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "adaptive_rho_interval"));
  settings->kkt_dynamic_reg = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "adaptive_rho_tolerance"));
  settings->abstol = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "adaptive_rho_fraction"));
  settings->reltol = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "max_iter"));
  settings->abstol_inacc = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "eps_abs"));
  settings->reltol_inacc = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "eps_rel"));
  settings->verbose = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "eps_dual_inf"));
  */
}