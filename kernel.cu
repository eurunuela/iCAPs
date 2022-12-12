/*
  Author : E. Orliac, SCITAS, EPFL
  Date   : 08.11.2017
  Purpose: 
  Remarks: 
*/

#include "cuda_runtime.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <device_launch_parameters.h>
#include <device_functions.h>
#include <math.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>


#define DB6   6
#define DB3   3

__device__ __constant__ double d_g6[DB6];
__device__ __constant__ double d_h6[DB6];
__device__ __constant__ double d_unity[1];


/*
__inline__ __device__ void reverse_vector(double *data, const int size) {

  unsigned   j;
  double   tmp;

  for(unsigned i=0; i<(size/2); i++) {
    j   = size - i - 1;
    tmp = data[j];
    data[j] = data[i];
    data[i] = tmp;
  }
} 
*/

__device__ void filter(double *out,
		       double *in,  const int n,
		       double *num, const int lnum,
		       double *den, const int lden,
		       const bool reverse) {
  int    i, j;
  double acc;

  if (lden == 1 && den[0] == 1.0) {
    if (reverse) {
      for(i=0; i<n; i++) {
	out[n-1-i] = num[0]*in[n-1-i];
	for(j=1; j<lnum; j++) {
	  if (j<=i) {
	    out[n-1-i] += num[j]*in[n-1-i+j];
	  }
	}
      }
    } else {
      /*
      for(i=0; i<n; i++) {
	out[i] = 0.0;
	for(j=0; j<lnum; j++) {
	  if (j<=i) {
	    out[i] += num[j]*in[i-j];
	  }
	}
      }
      */
      for(i=0; i<n; i+=2) {
	out[i]   = num[0]*in[i];
	out[i+1] = num[0]*in[i+1];
	for(j=1; j<lnum; j++) {
	  if (j<=i) {
	    out[i] += num[j]*in[i-j];
	  }
	  if (j<=i+1) {
	    out[i+1] += num[j]*in[i-j+1];
	  }
	}
      }
    }
  } else if (lnum == 1 && num[0] == 1.0) {

    for(i=0; i<n; i++)
      out[i]   = 0.0;

    if (reverse) {
      for(i=0; i<n; i++) {
	acc = in[n-1-i];
	for(j=1; j<lden; j++) {
	  if (j<=i) {
	    acc -= den[j]*out[i-j];
	  }
	}
	out[i] = acc;
      }
    } else {
      for(i=0; i<n; i++) {
	acc = in[i];
	for(j=1; j<lden; j++) {
	  if (j<=i) {
	    acc -= den[j]*out[i-j];
	  }
	}
	out[i] = acc;
      }
    }
  } else {
    printf("Fatal. Unknown case in kernel:filter.\n");	
    assert(false);
  }
}

__device__ void filter_boundary(double *tmp,
				double *out,
				double *data, const int n,
				double *num,  const int lnum,
				const int lden, 
				double *den1, const int lden1,
				double *den2, const int lden2,
				const int condition) {

  unsigned   j  = 0;
  double   dtmp = 0.0;

  if (condition == 1) {
    filter(out, data, n, num, lnum, d_unity, 1, true);
  } else if (condition == 0) {
    filter(out, data, n, num, lnum, d_unity, 1, false);
  } else {
    printf("Fatal. Unknown case in kernel filter_boundary.\n");
    assert(false);
  }
 
  if (lden == 2) {
 
    /* den1 is causal; den2 is non-causal */
    if (lden1+lden2 > 2) {

      // shiftnc = lden2 - 1
      if (lden2 - 1 != 0) {
	printf("FATAL. Non-zero shiftnc! Need to add missing implementation.\n");
	assert(false);
      }

      if (condition == 0) {
	filter(tmp, out, n, d_unity, 1, den1, lden1, false);
	filter(out, tmp, n, d_unity, 1, den2, lden2, true);
	for (unsigned i=0; i<n; i+=2) {
	  out[i]   *= den2[lden2-1];
	  out[i+1] *= den2[lden2-1];
	}
      } else {
	filter(tmp, out, n, d_unity, 1, den1, lden1, true);
	filter(out, tmp, n, d_unity, 1, den2, lden2, true);

	//EO: combine revert and scaling
	for(unsigned i=0; i<(n/2); i++) {
	  j      = n - i - 1;
	  dtmp   = out[j];
	  out[j] = out[i] * den2[lden2-1];
	  out[i] = dtmp   * den2[lden2-1];
	}
      }
    } else {
      printf("fatal: lden1+lden2 <= 2.");
      assert(false);
    }
  } else {
    printf("fatal: lden not equal to 2.");
    assert(false);
  }
}

__global__ void my_kernel(double* __restrict__ tcIn,
			  const int tcLength,
			  const int voxelNb,
			  double *num,
			  const int lnum,
			  const int lden,
			  double *den1, const int lden1,
			  double *den2, const int lden2,
			  const double lambdaTemp,
			  const double maxeig,
			  const int cost_save,
			  const int nit,
			  const double *noiseFinIn,
			  double *tcOut,
			  double *noiseFinOut) {

  double acc, nv, noise_estimate, precision, t_l, lambda_t;

  double lambda   = 0.0;
  double t        = 1.0;
  double maxeig_  = 1.0 / maxeig;
  double tclsqrt_ = sqrt(1.0 / tcLength);

  double *z   = (double*) malloc(tcLength*sizeof(double));
  assert(z != NULL);
  double *z_l = (double*) malloc(tcLength*sizeof(double));
  assert(z_l != NULL);
  double *z3  = (double*) malloc(tcLength*sizeof(double));
  assert(z3 != NULL);
  double *s   = (double*) malloc(tcLength*sizeof(double));
  assert(s != NULL);
  double *tmp = (double*) malloc(tcLength*sizeof(double));
  assert(tmp != NULL);

  unsigned    tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= voxelNb)
    return;

  unsigned tcinid = tid * tcLength;

  /* 
     Compute Daubechies 6 level 1 detailed coefficients
     Replicates Matlab sequence (see MyTemporal.m):
         [coef,  len] = wavedec(s,1,'db3');
         coef(1:len(1)) = [];
   ! Assumes default dwtmode('sym') is set.
   ! Requires even length input signal.
  */

  // Requires even length input signal
  assert(tcLength%2==0);
  
  // Number of detailed coefficients to be computed
  unsigned nc = floor((tcLength-1)*0.5) + DB3;

  double *coeff = (double*) malloc(nc * sizeof(double));
  assert(coeff != NULL);

  // Number of coefficients that do not depend on padded values
  unsigned n1 = nc - 4;

  coeff[0]  = 0.0;
  coeff[0] += d_g6[0] * tcIn[tcinid + 3];
  coeff[0] += d_g6[1] * tcIn[tcinid + 2];
  coeff[0] += d_g6[2] * tcIn[tcinid + 1];
  coeff[0] += d_g6[3] * tcIn[tcinid];
  coeff[0] += d_g6[4] * tcIn[tcinid];
  coeff[0] += d_g6[5] * tcIn[tcinid + 1];

  coeff[1]  = 0.0;
  coeff[1] += d_g6[0] * tcIn[tcinid + 1];
  coeff[1] += d_g6[1] * tcIn[tcinid];
  coeff[1] += d_g6[2] * tcIn[tcinid];
  coeff[1] += d_g6[3] * tcIn[tcinid + 1];
  coeff[1] += d_g6[4] * tcIn[tcinid + 2];
  coeff[1] += d_g6[5] * tcIn[tcinid + 3];

  for (unsigned i=0; i < n1; i++) {
    coeff[i+2] = 0.0;
    for (unsigned j=0; j<DB6; j++) {
      coeff[i+2] += d_g6[j] * tcIn[tcinid + 2*i + j];
    }
  }

  coeff[n1+2]  = 0.0;
  coeff[n1+2] += d_g6[0] * tcIn[tcinid + tcLength - 4];
  coeff[n1+2] += d_g6[1] * tcIn[tcinid + tcLength - 3];
  coeff[n1+2] += d_g6[2] * tcIn[tcinid + tcLength - 2];
  coeff[n1+2] += d_g6[3] * tcIn[tcinid + tcLength - 1];
  coeff[n1+2] += d_g6[4] * tcIn[tcinid + tcLength - 1];
  coeff[n1+2] += d_g6[5] * tcIn[tcinid + tcLength - 2];

  coeff[n1+3]  = 0.0;
  coeff[n1+3] += d_g6[0] * tcIn[tcinid + tcLength - 2];
  coeff[n1+3] += d_g6[1] * tcIn[tcinid + tcLength - 1];
  coeff[n1+3] += d_g6[2] * tcIn[tcinid + tcLength - 1];
  coeff[n1+3] += d_g6[3] * tcIn[tcinid + tcLength - 2];
  coeff[n1+3] += d_g6[4] * tcIn[tcinid + tcLength - 3];
  coeff[n1+3] += d_g6[5] * tcIn[tcinid + tcLength - 4];

  thrust::sort(thrust::seq, coeff, coeff+nc);

  double median = 0.0;
  const size_t lhs = (nc - 1) * 0.5;
  const size_t rhs = nc * 0.5;
  if (lhs == rhs) {
    median = coeff[lhs];
  } else {
    median = (coeff[lhs] + coeff[rhs]) * 0.5;
  }

  for(unsigned i=0; i<nc; i++) {
    coeff[i] = fabs(coeff[i] - median);
  }

  thrust::sort(thrust::seq, coeff, coeff+nc);
  if (lhs == rhs) {
    median = coeff[lhs];
  } else {
    median = (coeff[lhs] + coeff[rhs]) * 0.5;
  }

  free(coeff);

  __syncthreads();


  lambda_t       = median * lambdaTemp;
  noise_estimate = lambda_t;
  precision      = noise_estimate * 1E-5;

  if (noiseFinIn[tid] == 0.0) {
    //printf("case 0.0 in cuda: lambda = %15.10f\n", lambda_t);
    lambda = lambda_t;
  } else {
    //printf("case != 0.0 in cuda: lambda = %15.10f (vox %i)\n", noiseFinIn[tid], tid);
    lambda = noiseFinIn[tid];
  }

  for (unsigned i=0; i<tcLength; i+=2) {
    s[i]   = 0.0;
    s[i+1] = 0.0;
    z[i]   = tcIn[tcinid + i];
    z[i+1] = tcIn[tcinid + i + 1];
  }

  filter_boundary(tmp, &tcOut[tcinid], &tcIn[tcinid], tcLength, num, lnum, lden, den1, lden1, den2, lden2, 0);

  nv = 0.0;

  for (unsigned k=0; k<nit; k++) {
    
    for(unsigned i=0; i<tcLength; i++)
      z_l[i] = z[i];

    filter_boundary(tmp, z,  s, tcLength, num, lnum, lden, den1, lden1, den2, lden2, 1);
    filter_boundary(tmp, z3, z, tcLength, num, lnum, lden, den1, lden1, den2, lden2, 0);

    for (unsigned i=0; i<tcLength; i++) {
      z[i] = maxeig_ * (tcOut[tcinid+i]/lambda - z3[i]) + s[i];
      if (z[i] >  1.0) z[i] =  1.0;
      if (z[i] < -1.0) z[i] = -1.0;
    }
    
    t_l = t;
    t   = 0.5 + sqrt(0.25 + t*t);
    double ct = (t_l - 1.0)/t;

    for (unsigned i=0; i<tcLength; i++)
      s[i] = z[i] + ct * (z[i]-z_l[i]);

    nv = 0.0;
    if (cost_save) {
      printf("Fatal. Missing implementation if cost_save == 1 in Temporal_TA_MEX.\n");
      assert(false);
    } else {
      filter_boundary(tmp, z3, z, tcLength, num, lnum, lden, den1, lden1, den2, lden2, 1);
      acc = 0.0;
      for(unsigned i=0; i<tcLength; i++)
	acc += z3[i]*z3[i];
      nv = lambda * sqrt(acc) * tclsqrt_;
    }
    
    if (fabs(nv-noise_estimate) > precision) {
      lambda *= noise_estimate/nv;
    }
  }

  for(unsigned i=0; i<tcLength; i++)
    tcOut[tcinid + i] = tcIn[tcinid + i] - lambda*z3[i];

  noiseFinOut[tid] = nv;

  free(s);
  free(z3);
  free(z_l);
  free(z);
  free(tmp);
}


// Function to set constant memory for the Daubechies db6 coefficients
void set_constant_memory(const double *g6, const double *h6, const double *one) {
  cudaMemcpyToSymbol(d_g6,    g6,  DB6  * sizeof(double));
  cudaMemcpyToSymbol(d_h6,    h6,  DB6  * sizeof(double));
  cudaMemcpyToSymbol(d_unity, one,        sizeof(double));
}

