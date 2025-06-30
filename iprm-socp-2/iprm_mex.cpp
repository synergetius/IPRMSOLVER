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

// wrapper class for all iprm data and settings
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
void castToDoubleArr(IPRMFloat *arr, double *arr_out, IPRMInt len);
void setToNaN(double *arr_out, IPRMInt len);
void copyMxStructToSettings(const mxArray *, IPRMSettings *);
void copyUpdatedSettingsToWork(const mxArray *, IPRMSettings *);
void castCintToDoubleArr(IPRMInt *arr, double *arr_out, IPRMInt len);
IPRMInt *copyToCintVector(mwIndex *vecData, IPRMInt numel);
IPRMInt *copyDoubleToCintVector(double *vecData, IPRMInt numel);
IPRMInt *roundDoubleIntoInt(double *vecData, IPRMInt numel);
IPRMFloat *copyToCfloatVector(double *vecData, IPRMInt numel);
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
  IPRMInt exitflag = 1;

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

    const IPRMInt n = (IPRMInt)mxGetScalar(prhs[2]);
    const IPRMInt m = (IPRMInt)mxGetScalar(prhs[3]);
    const IPRMInt p = (IPRMInt)mxGetScalar(prhs[4]);
    const mxArray *P = prhs[5];
    const mxArray *c = prhs[6];
    const mxArray *A = prhs[7];
    const mxArray *b = prhs[8];
    const mxArray *G = prhs[9];
    const mxArray *h = prhs[10];
    const IPRMInt l = (IPRMInt)mxGetScalar(prhs[11]);
    const IPRMInt nsoc = (IPRMInt)mxGetScalar(prhs[12]);
    const mxArray *q = prhs[13];

    // double *data = mxGetPr(b); // 获取数据指针
    // mwSize n0 = mxGetM(b);      // 获取行数
    // double max_val = 0.0;
    // 
    // for (mwSize i = 0; i < n0; i++) {
    //     double abs_val = fabs(data[i]);
    //     if (abs_val > max_val) {
    //         max_val = abs_val;
    //     }
    // }
    // printf("b inf_norm: %.12e\n", max_val);
    // printf("p0 = %d, p = %d\n", n0, p);
    // data = mxGetPr(c); // 获取数据指针
    // n0 = mxGetM(c);      // 获取行数
    // max_val = 0.0;
    // 
    // for (mwSize i = 0; i < n0; i++) {
    //     double abs_val = fabs(data[i]);
    //     if (abs_val > max_val) {
    //         max_val = abs_val;
    //     }
    // }
    // printf("c inf_norm: %.12e\n", max_val);
    // printf("n0 = %d, n = %d\n", n0, n);

    // Create Data Structure
    IPRMFloat *cvec = copyToCfloatVector(mxGetPr(c), n);
      printf("cvec: %.12e\n", inf_norm(cvec, n));
    IPRMFloat *bvec = copyToCfloatVector(mxGetPr(b), p);
    printf("bvec: %.12e\n", inf_norm(bvec, p));
    IPRMFloat *hvec = copyToCfloatVector(mxGetPr(h), m);
    IPRMInt *qvec = roundDoubleIntoInt(mxGetPr(q), nsoc); // Warning: This only works with IPRMInt is int32

    // Matrix P:  nnz = P->p[n]
    IPRMCscMatrix *Pcsc = new IPRMCscMatrix;
    IPRMInt *Pp = copyToCintVector(mxGetJc(P), n + 1);
    IPRMInt *Pi = copyToCintVector(mxGetIr(P), Pp[n]);
    IPRMFloat *Px = copyToCfloatVector(mxGetPr(P), Pp[n]);
    iprm_set_csc(Pcsc, n, n, Pp[n], Px, Pp, Pi);

    // // Matrix A: nnz = A->p[n]
    IPRMCscMatrix *Acsc = new IPRMCscMatrix;
    IPRMInt *Ap = copyToCintVector(mxGetJc(A), n + 1);
    IPRMInt *Ai = copyToCintVector(mxGetIr(A), Ap[n]);
    IPRMFloat *Ax = copyToCfloatVector(mxGetPr(A), Ap[n]);
    iprm_set_csc(Acsc, p, n, Ap[n], Ax, Ap, Ai);

    // Matrix G: nnz = G->p[n]
    IPRMCscMatrix *Gcsc = new IPRMCscMatrix;
    IPRMInt *Gp = copyToCintVector(mxGetJc(G), n + 1);
    IPRMInt *Gi = copyToCintVector(mxGetIr(G), Gp[n]);
    IPRMFloat *Gx = copyToCfloatVector(mxGetPr(G), Gp[n]);
    iprm_set_csc(Gcsc, m, n, Gp[n], Gx, Gp, Gi);

    // Create Settings
    const mxArray *mxSettings = prhs[14]; ///////////
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
    exitflag = iprm_setup(iprmData->solver, n, m, p, Pcsc, cvec, Acsc, bvec, Gcsc, hvec, l, nsoc, qvec, settings);

    // cleanup temporary structures
    if (cvec)
      iprm_free(cvec);
    if (bvec)
      iprm_free(bvec);
    if (hvec)
      iprm_free(hvec);
    if (qvec)
      iprm_free(qvec);
    if (Px)
      iprm_free(Px);
    if (Pi)
      iprm_free(Pi);
    if (Pp)
      iprm_free(Pp);
    if (Pcsc)
      delete Pcsc;
    if (Ax)
      iprm_free(Ax);
    if (Ai)
      iprm_free(Ai);
    if (Ap)
      iprm_free(Ap);
    if (Acsc)
      delete Acsc;
    if (Gx)
      iprm_free(Gx);
    if (Gi)
      iprm_free(Gi);
    if (Gp)
      iprm_free(Gp);
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

    if (!strcmp("IPRM_UNSOLVED", constant))
    {
      plhs[0] = mxCreateDoubleScalar(IPRM_UNSOLVED);
      return;
    }

    if (!strcmp("IPRM_SOLVED", constant))
    {
      plhs[0] = mxCreateDoubleScalar(IPRM_SOLVED);
      return;
    }

    if (!strcmp("IPRM_SOLVED_INACCURATE", constant))
    {
      plhs[0] = mxCreateDoubleScalar(IPRM_SOLVED_INACCURATE);
      return;
    }

    if (!strcmp("IPRM_NUMERICAL_ERROR", constant))
    {
      plhs[0] = mxCreateDoubleScalar(IPRM_NUMERICAL_ERROR);
      return;
    }

    if (!strcmp("IPRM_MAX_ITER", constant))
    {
      plhs[0] = mxCreateDoubleScalar(IPRM_MAX_ITER);
      return;
    }

    mexErrMsgTxt("Constant not recognized.");

    return;
  }

  // Got here, so command not recognized
  mexErrMsgTxt("Command not recognized.");
}

IPRMFloat *copyToCfloatVector(double *vecData, IPRMInt numel)
{
  // This memory needs to be freed!
  IPRMFloat *out = (IPRMFloat *)iprm_malloc(numel * sizeof(IPRMFloat));

  // copy data
  for (IPRMInt i = 0; i < numel; i++)
  {
    out[i] = (IPRMFloat)vecData[i];
  }
  return out;
}

IPRMInt *roundDoubleIntoInt(double *vecData, IPRMInt numel)
{
  // This memory needs to be freed!
  IPRMInt *out = (IPRMInt *)iprm_malloc(numel * sizeof(IPRMInt));

  // round data
  for (IPRMInt i = 0; i < numel; i++)
  {
    out[i] = round(vecData[i]);
  }
  return out;
}

IPRMInt *copyToCintVector(mwIndex *vecData, IPRMInt numel)
{
  // This memory needs to be freed!
  IPRMInt *out = (IPRMInt *)iprm_malloc(numel * sizeof(IPRMInt));

  // copy data
  for (IPRMInt i = 0; i < numel; i++)
  {
    out[i] = (IPRMInt)vecData[i];
  }
  return out;
}

IPRMInt *copyDoubleToCintVector(double *vecData, IPRMInt numel)
{
  // This memory needs to be freed!
  IPRMInt *out = (IPRMInt *)iprm_malloc(numel * sizeof(IPRMInt));

  // copy data
  for (IPRMInt i = 0; i < numel; i++)
  {
    out[i] = (IPRMInt)vecData[i];
  }
  return out;
}

void castCintToDoubleArr(IPRMInt *arr, double *arr_out, IPRMInt len)
{
  for (IPRMInt i = 0; i < len; i++)
  {
    arr_out[i] = (double)arr[i];
  }
}

void castToDoubleArr(IPRMFloat *arr, double *arr_out, IPRMInt len)
{
  for (IPRMInt i = 0; i < len; i++)
  {
    arr_out[i] = (double)arr[i];
  }
}

void setToNaN(double *arr_out, IPRMInt len)
{
  IPRMInt i;
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

  settings->max_iters = (IPRMInt)mxGetScalar(mxGetField(mxPtr, 0, "max_iters"));
  settings->iter_ref_iters = (IPRMInt)mxGetScalar(mxGetField(mxPtr, 0, "iter_ref_iters"));
  settings->linesearch_iters = (IPRMInt)mxGetScalar(mxGetField(mxPtr, 0, "linesearch_iters"));
  settings->kkt_static_reg = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "kkt_static_reg"));
  settings->kkt_dynamic_reg = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "kkt_dynamic_reg"));
  settings->mu0 = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "mu0"));
  settings->rho0 = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "rho0"));
  settings->eta = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "eta"));
  settings->gamma0 = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "gamma0"));
  settings->delta = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "delta"));
  settings->tau = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "tau"));
  settings->sigma = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "sigma"));
  settings->epsilon = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "epsilon"));
  settings->verbose = (IPRMInt)mxGetScalar(mxGetField(mxPtr, 0, "verbose"));
}

void copyUpdatedSettingsToWork(const mxArray *mxPtr, IPRMSettings *settings)
{
  settings->max_iters = (IPRMInt)mxGetScalar(mxGetField(mxPtr, 0, "max_iters"));
  settings->iter_ref_iters = (IPRMInt)mxGetScalar(mxGetField(mxPtr, 0, "iter_ref_iters"));
  settings->linesearch_iters = (IPRMInt)mxGetScalar(mxGetField(mxPtr, 0, "linesearch_iters"));
  settings->kkt_static_reg = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "kkt_static_reg"));
  settings->kkt_dynamic_reg = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "kkt_dynamic_reg"));
  settings->mu0 = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "mu0"));
  settings->rho0 = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "rho0"));
  settings->eta = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "eta"));
  settings->gamma0 = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "gamma0"));
  settings->delta = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "delta"));
  settings->tau = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "tau"));
  settings->sigma = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "sigma"));
  settings->epsilon = (IPRMFloat)mxGetScalar(mxGetField(mxPtr, 0, "epsilon"));
  settings->verbose = (IPRMInt)mxGetScalar(mxGetField(mxPtr, 0, "verbose"));
}