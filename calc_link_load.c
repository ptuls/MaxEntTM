/* ============================================================= */
/* === MATLAB/calc_link_load mexFunction ======================= */
/* ============================================================= */

/* ----------------------------------------------------------------
 * [L] = calc_link_load(D,W,P,T,mode)
 *
 * where:
 *
 *    D is the distance matrix
 *    P is the predecessor matrix
 *    T is the traffic matrix
 *
 *    L is the link load matrix
 *
 * Copyright (c) 2005-2006 Yin Zhang <yzhang@cs.utexas.edu>
 * ----------------------------------------------------------------
 */

/* $Header: /u/yzhang/PrivacyTE/src/Package/RCS/calc_link_load.c,v 1.2 2005/12/29 13:27:00 yzhang Exp $ */

#include <math.h>
#include <stdlib.h>  /* needed for qsort() */
#include "mex.h"
#include "matrix.h"

/* Compute link load solely based on the predecessor matrix */
void compute_load_from_pred(int N, double**P, double**D, double**T, double **L)
{
  int i,j,k,p;
  double t;
  /*
   * compute the link loads
   * XXX be careful! C, D, P, L are stored in column order
   */
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      if (!mxIsInf(D[j][i]) && D[j][i] != 0) { /* path from node i to node j */
        t = T[j][i];           /* traffic from i to j */
        k = j;
        p = (int) P[k][i] - 1;
        while (p >= 0) {
          /* add t to link load k --> p */
          L[k][p] += t;
          /* update k and p */
          k = p;
          p = (int) P[k][i] - 1;
        }
/* 	if (p != i) { */
/* 	  mexErrMsgTxt("pred didn't get back to start i"); */
/* 	} */
      }
    }
  }
}
  
typedef struct dist_s {
  int    node;
  double dist;
} dist_t;

static int compare_dist(const void *d1, const void *d2)
{
  double v1, v2;
  v1 = ((dist_t*) d1)->dist;
  v2 = ((dist_t*) d2)->dist;
  if (v1 < v2)
    return 1;
  else if (v1 > v2)
    return -1;
  else
    return 0;
}

/* Compute link load by simulating OSPF routing */
#define EQUAL_TOL 1.0e-12  
void compute_load_with_ospf(int N, double**W, double**D, double**T, double **L)
{
  int deg,i,j,k,m,n,*node_lst;
  double *traffic,w,t;
  dist_t *dist_lst;

  traffic  = (double*) mxMalloc(N*sizeof(double));
  dist_lst = (dist_t*) mxMalloc(N*sizeof(dist_t));
  node_lst = (int*)    mxMalloc(N*sizeof(int));
  
  for (i = 0; i < N; i++) { /* reaching i */
    for (j = 0; j < N; j++) { /* from j */
      traffic[j]       = T[i][j];       /* traffic(j->i) */
      dist_lst[j].node = j;
      dist_lst[j].dist = D[i][j];
    }

    /* sort all j in decreasing order of dist(j->i) */
    qsort(dist_lst, N, sizeof(dist_t), compare_dist);

    for (m = 0; m < N; m++) {
      j = dist_lst[m].node;
      if (mxIsInf(D[i][j]) || (D[i][j] == 0)) continue;
    
      /* compute the splitting factor */
      deg = 0;
      for (k = 0; k < N; k++) {
        /* there is an edge j->k and j->k + k->i == j->i */
        w = W[k][j];
        if (!mxIsInf(w) && (w != 0) && 
            fabs(w + D[i][k] - D[i][j]) < EQUAL_TOL) {
          node_lst[deg++] = k;
        }
      }

      /* add traffic[j] to all the links */
      t = traffic[j]/deg;
      for (n = 0; n < deg; n++) {
        k = node_lst[n];
        traffic[k] += t;
        L[k][j]    += t;
      }
    }
  }

  mxFree(traffic);
  mxFree(dist_lst);
  mxFree(node_lst);
}

void mexFunction(
        int nlhs,       mxArray *plhs[],
        int nrhs, const mxArray *prhs[]
        )
{
  /* Declare variables */
  int m,n,i,j,ospf;
  double *Dx,*Px,*Tx,*Wx,*Lx,*Ox,**D,**P,**T,**W,**L;

  /* Check for proper number of input and output arguments. */    
  if (nrhs != 4 && nrhs != 5) {
    mexErrMsgTxt("Five input arguments required: L = calc_link_load(D,W,P,T,[ospf])");
  }
  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments.");
  }

  /* Check data type of input argument. */
  if (mxGetNumberOfDimensions(prhs[0]) != 2 ||
      mxGetNumberOfDimensions(prhs[1]) != 2 ||
      mxGetNumberOfDimensions(prhs[2]) != 2 ||
      mxGetNumberOfDimensions(prhs[3]) != 2) {
    mexErrMsgTxt("First 4 input arguments must be two dimensional\n");
  }
  if (mxIsChar(prhs[0]) || mxIsSparse(prhs[0]) || mxIsComplex(prhs[0]) ||
      mxIsChar(prhs[1]) || mxIsSparse(prhs[1]) || mxIsComplex(prhs[1]) ||
      mxIsChar(prhs[2]) || mxIsSparse(prhs[2]) || mxIsComplex(prhs[2]) ||
      mxIsChar(prhs[3]) || mxIsSparse(prhs[3]) || mxIsComplex(prhs[3])) {
    mexErrMsgTxt("First 4 input arguments must be full real matrices.");
  }
 
  /* Get the size and pointers to input data. */
  m  = mxGetM(prhs[0]);
  n  = mxGetN(prhs[0]);
  if (m != n ||
      mxGetM(prhs[1]) != n || mxGetN(prhs[1]) != n ||
      mxGetM(prhs[2]) != n || mxGetN(prhs[2]) != n ||
      mxGetM(prhs[3]) != n || mxGetN(prhs[3]) != n) {
    mexErrMsgTxt("First 4 input arguments must be square and have the same size.");
  }
  Dx  = mxGetPr(prhs[0]);
  Wx  = mxGetPr(prhs[1]);
  Px  = mxGetPr(prhs[2]);
  Tx  = mxGetPr(prhs[3]);

  ospf = 1;
  if (nrhs == 5) {
    if (mxGetM(prhs[4]) != 1 || mxGetN(prhs[4]) != 1) {
      mexErrMsgTxt("Last input arguments must be a single number.");
    }
    Ox = mxGetPr(prhs[4]);
    ospf = (Ox[0] == 0) ? 0 : 1;
  }

  /* create output matrices */
  plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
  Lx = mxGetPr(plhs[0]);

  /* set up the 2-d arrays */
  D  = (double**) mxMalloc(n*sizeof(double*));
  W  = (double**) mxMalloc(n*sizeof(double*));
  P  = (double**) mxMalloc(n*sizeof(double*));
  T  = (double**) mxMalloc(n*sizeof(double*));
  L  = (double**) mxMalloc(n*sizeof(double*));
  for (i = 0, j = 0; i < n; i++, j+=n) {
    D[i] = Dx + j;
    W[i] = Wx + j;
    P[i] = Px + j;
    T[i] = Tx + j;
    L[i] = Lx + j;
  }

  for (i = 0; i < n*n; i++) {
    Lx[i] = 0.0;
  }
  
  if (ospf) {
    compute_load_with_ospf(n,W,D,T,L);
  } else {
    compute_load_from_pred(n,P,D,T,L);
  }
  
  /* garbage collection */
  mxFree(D);
  mxFree(W);
  mxFree(P);
  mxFree(T);
  mxFree(L);
}
