/*=========================================================
 * mldivide_lapack.c - Replacement for the mldivide MATLAB
 * function by calling LAPACK dgels directly.
 *
 * X = mldivide_lapack(A,B) computes the solution to a 
 * system of linear equations A * X = B
 * using LAPACK routine DGELS, where 
 * A is a real M-by-N matrix.
 * X and B are real M-by-1 matrices.
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2009-2018 The MathWorks, Inc.
 *=======================================================*/

#if !defined(_WIN32)
#define dgels dgels_
#endif

#include <string.h> /* needed for memcpy() */
#include "mex.h"
#include "lapack.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *A, *B, *work;
    double *A2, *B2;
    size_t m,n,nrhsB;
    ptrdiff_t info, lwork;
    char *chn = "N";

    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    nrhsB = mxGetN(prhs[1]);

#ifdef LAPACK_VERBOSE
    printf("m = %d n = %d nrhsB = %d\n", m, n, nrhsB);
#endif

    A = mxGetPr(prhs[0]); /* pointer to first input matrix */
    B = mxGetPr(prhs[1]); /* pointer to second input matrix */

    /* DGELS works in-place, so we copy the inputs first. */
    mxArray *Awork;
    Awork = mxCreateDoubleMatrix(m, n, mxREAL);
    A2 = mxGetPr(Awork);
    memcpy(A2, A, m*n*mxGetElementSize(prhs[0]));

    plhs[0] = mxCreateDoubleMatrix(m, nrhsB, mxREAL);
    B2 = mxGetPr(plhs[0]);
    memcpy(B2, B, m*nrhsB*mxGetElementSize(prhs[1]));

    /* Call LAPACK to find work array size*/
    double wkopt;
    lwork = -1;
    dgels(chn, &m, &n, &nrhsB, A2, &m, B2, &m, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
#ifdef LAPACK_VERBOSE
    printf("info = %d lwork = %d\n", info, lwork);
#endif

    /* create work matrix */
    work = (double *)mxCalloc(lwork,sizeof(double));

    /* Call LAPACK to solve system */
    dgels(chn, &m, &n, &nrhsB, A2, &m, B2, &m, work, &lwork, &info);
#ifdef LAPACK_VERBOSE
    printf("info = %d\n", info);
#endif
    /* plhs[0] now holds X */

    mxDestroyArray(Awork);
    mxFree(work);
}
