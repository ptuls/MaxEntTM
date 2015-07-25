/* ============================================================= */
/* === MATLAB/floyd mexFunction ============================= */
/* ============================================================= */

/* ----------------------------------------------------------------
 * Compute all-pair shortest paths using Floyd's algorithm.
 *
 * Copyright (c) 2005-2006 Yin Zhang <yzhang@cs.utexas.edu>
 * ----------------------------------------------------------------
 */

/* $Header: /u/yzhang/PrivacyTE/src/Package/RCS/floyd_apsp.c,v 1.1 2005/12/29 05:24:18 yzhang Exp $ */

#include "mex.h"
#include "matrix.h"

void FloydAPSP (int N, double** C, double** D, double** P)
{
  int i,j,k;
#ifndef USE_ROW_ORDER /* USE_COLUMN_ORDER */
  /*
   * When C, D, P are stored in column order (as in MATLAB), we have
   *
   *   C[i][j] gives the weight from j to i
   *   D[i][j] gives the distance from j to i
   *   P[i][j] gives the predecessor of i (from j to i)
   */

  /* initialization */
  for (j = 0; j < N; j++) {
    for (i = 0; i < N; i++) {
      D[j][i] = C[j][i];
      P[j][i] = i;
    }
    D[j][j] = 0.0;
    P[j][j] = -1;
  }
  for (k = 0; k < N; k++) { /* k-> is the intermediate point */
    for (j = 0; j < N; j++) { /* reaching j */
      for (i = 0; i < N; i++) { /* start from i */
        /* if i-->k + k-->j is smaller than original i-->j */
        if (mxIsInf(D[j][k])) continue;
        if (D[k][i] + D[j][k] < D[j][i]) {
          /* reduce the i-->j distance to the smaller one i->k->j */
          D[j][i] = D[k][i] + D[j][k];
          /* and update the predecessor matrix */
          P[j][i] = P[j][k];
    	}
      }
    }
  }

#else /* USE_ROW_ORDER */
  /* initialization */
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      D[i][j] = C[i][j];
      P[i][j] = i;
    }
    D[i][i] = 0.0;
    P[i][i] = -1;
  }
  for (k = 0; k < N; k++) { /* k-> is the intermediate point */
    for (i = 0; i < N; i++) { /* start from i */
      if (mxIsInf(D[i][k])) continue;
      for (j = 0; j < N; j++) { /* reaching j */
        /* if i-->k + k-->j is smaller than original i-->j */
        if (D[i][k] + D[k][j] < D[i][j]) {
          /* reduce the i-->j distance to the smaller one i->k->j */
          D[i][j] = D[i][k] + D[k][j];
          /* and update the predecessor matrix */
          P[i][j] = P[k][j];
    	}
      }
    }
  }
#endif
} /* FloydAPSP */


void mexFunction(
        int nlhs,       mxArray *plhs[],
        int nrhs, const mxArray *prhs[]
        )
{
  /* Declare variables */
  int m,n,i,j;
  double *Cx,*Dx,*Px,**C,**D,**P;

  /* Check for proper number of input and output arguments. */    
  if (nrhs != 1) {
    mexErrMsgTxt("One input arguments required: [D,P] = floyd(C)");
  }
  if (nlhs > 2) {
    mexErrMsgTxt("Too many output arguments.");
  }

  /* Check data type of input argument. */
  if (mxGetNumberOfDimensions(prhs[0]) != 2) {
    mexErrMsgTxt("Input argument must be two dimensional\n");
  }
  if (mxIsChar(prhs[0]) || mxIsSparse(prhs[0]) || mxIsComplex(prhs[0])) {
    mexErrMsgTxt("Input argument must be a full real matrix.");
  }
 
  /* Get the size and pointers to input data. */
  m  = mxGetM(prhs[0]);
  n  = mxGetN(prhs[0]);
  if (m != n) {
    mexErrMsgTxt("Input argument must be square.");
  }
  Cx  = mxGetPr(prhs[0]);

  /* create output matrices */
  plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  Dx = mxGetPr(plhs[0]);

  plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
  Px = mxGetPr(plhs[1]);

  /* set up the 2-d arrays */
  C  = (double**) mxMalloc(n*sizeof(double*));
  D  = (double**) mxMalloc(n*sizeof(double*));
  P  = (double**) mxMalloc(n*sizeof(double*));
  for (i = 0, j = 0; i < n; i++, j+=n) {
    C[i] = Cx + j;
    D[i] = Dx + j;
    P[i] = Px + j;
  }

  FloydAPSP(n, C, D, P);

  /* translate base 0 into base 1 */
  for (i = 0; i < n*n; i++) {
    Px[i]++;
  }

  /* garbage collection */
  mxFree(C);
  mxFree(D);
  mxFree(P);
}
