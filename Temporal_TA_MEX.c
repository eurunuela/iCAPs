/*
  Author : E. Orliac, SCITAS, EPFL
  Date   : 27.10.2017
  Purpose: 
  Remarks: 
*/

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "filter_boundary_MEX.h"
#include "daubechies_6d.h"


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *TCIN;            /* Input signal */
  double *num;             /* param.f_Analyze.num */
  double noise_estimate;   /* ParametersIn.LambdaTemp(ParametersIn.VxlInd) */
  double maxeig;           /* ParametersIn.MaxEig */
  int    N;                /* ParametersIn.Dimension(4) */ 
  int    Nit;              /* */
  int    den_length;       /* Length of param.f_Analyze.den */
  double *den1;            /* denominator, causal part */
  double *den2;            /* denominator, non-causal part */
  double lambda;           /* */
  bool   cost_save;        /* */
  double *TCOUT;           /* */
  double *paramOut;        /* */

  double *C;               /* Debauchies DB6 wavelet coefficients */

  int i, k, ncd1, nd1, nd2;
  double med, precision, *out, *z, *z2, *z3, *s, t, t_l, acc, nv;
  double *z_l, *test;

  const int DB = 3;

  TCIN           = mxGetPr(prhs[0]);
  num            = mxGetPr(prhs[1]);
  noise_estimate = mxGetScalar(prhs[2]);
  maxeig         = mxGetScalar(prhs[3]);
  N              = (int) mxGetScalar(prhs[4]);
  Nit            = (int) mxGetScalar(prhs[5]);
  den_length     = (int) mxGetScalar(prhs[6]);
  den1           = mxGetPr(prhs[7]);
  den2           = mxGetPr(prhs[8]);
  lambda         = mxGetScalar(prhs[9]);
  cost_save      = (int) mxGetScalar(prhs[10]);
  TCOUT          = mxGetPr(prhs[11]);
  paramOut       = mxGetPr(prhs[12]);

  const mwSize *dimsTCIN = mxGetDimensions(prhs[0]);
  const mwSize *dimsNum  = mxGetDimensions(prhs[1]);
  const mwSize *dimsDen1 = mxGetDimensions(prhs[7]);
  const mwSize *dimsDen2 = mxGetDimensions(prhs[8]);

  /*
  printf("Length of X   = %d x %d\n", dimsTCIN[0], dimsTCIN[1]);
  printf("Length of num = %d x %d\n", dimsNum[0], dimsNum[1]);
  printf("num           = [%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f]\n",
	 num[0], num[1], num[2], num[3], num[4], num[5]);
  printf("noise_estimate = %.10f\n", noise_estimate);
  printf("maxeig         = %.10f\n", maxeig);
  printf("N              = %d\n", N);
  printf("Nit            = %d\n", Nit);
  printf("den1           = %d x %d\n", dimsDen1[0], dimsDen1[1]);
  printf("den2           = %d x %d\n", dimsDen2[0], dimsDen2[1]);
  printf("lambda         = %.10f\n", lambda);
  printf("cost_save      = %d\n", cost_save);
  */

  precision = noise_estimate/100000.d;

  z      = (double*) calloc(sizeof(double), N);
  s      = (double*) calloc(sizeof(double), N);
  z2     = (double*) calloc(sizeof(double), N);
  z3     = (double*) calloc(sizeof(double), N);
  z_l    = (double*) calloc(sizeof(double), N);
  t      = 1.d;


  for(i=0; i<dimsTCIN[0]; i++)
    z[i] = TCIN[i];

  for(k=0; k<Nit; k++) {

    for(i=0; i<dimsTCIN[0]; i++)
      z_l[i] = z[i];

    filter_boundary_MEX(z,  TCIN , N, num, dimsNum[1], den_length, den1, dimsDen1[1], den2, dimsDen2[0], 0);
    filter_boundary_MEX(z2, s,     N, num, dimsNum[1], den_length, den1, dimsDen1[1], den2, dimsDen2[0], 1);
    filter_boundary_MEX(z3, z2,    N, num, dimsNum[1], den_length, den1, dimsDen1[1], den2, dimsDen2[0], 0);

    for(i=0; i<N; i++) {
      z[i] = 1.d/(lambda*maxeig)*z[i] + s[i] - z3[i]/maxeig;
      /* Clipping */
      if (z[i] > 1.d)  z[i] =  1.d;
      if (z[i] < -1.d) z[i] = -1.d;
    }

    t_l = t;
    t   = (1.d+sqrt(1.d+4*(t*t)))*0.5d;
    for(i=0; i<dimsTCIN[0]; i++)
      s[i] = z[i] + (t_l - 1.d)/t * (z[i]-z_l[i]);

    if (cost_save) {
      fprintf(stderr, "Fatal. Missing implementation if cost_save == 1 in Temporal_TA_MEX.\n");
      exit(EXIT_FAILURE);
    } else {
      filter_boundary_MEX(z2, z, N, num, dimsNum[1], den_length, den1, dimsDen1[1], den2, dimsDen2[0], 1);
      acc = 0.d;
      for(i=0; i<N; i++)
	acc += z2[i]*z2[i];
      nv = sqrt(lambda*lambda*acc/N);
    }

    if (fabs(nv-noise_estimate) > precision) {
      lambda = lambda*noise_estimate/nv;
    }
  }

  filter_boundary_MEX(z2, z, N, num, dimsNum[1], den_length, den1, dimsDen1[1], den2, dimsDen2[0], 1);

  for(i=0; i<N; i++)
    TCOUT[i] = TCIN[i] - lambda*z2[i];
 
  paramOut[0] = lambda;
  paramOut[1] = nv;

  free(z);
  free(s);
  free(z2);
  free(z3);
  free(z_l);
}

