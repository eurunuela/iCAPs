/*
  Author : E. Orliac, SCITAS, EPFL
  Date   : 08.11.2017
  Purpose: CUDA C MEX implementation of MyTemporal.m loop over voxels.
  Remarks: 
 */

#include "mex.h"
#include "cuda_runtime.h"
#include "kernel.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <device_launch_parameters.h>
#include <device_functions.h>
#include <math.h>

#define CHECK(call) { \
cudaError_t err; \
if ( (err = (call)) != cudaSuccess) { \
  fprintf(stderr, "Got error: %s at %s:%d\n", cudaGetErrorString(err), __FILE__, __LINE__); \
  exit(0); \
 } \
}

static const double h_6[6] = { 
   0.33267055295008261599851158914,
   0.80689150931109257649449360409,
   0.45987750211849157009515194215,
  -0.13501102001025458869638990670,
  -0.08544127388202666169281916918,
   0.03522629188570953660274066472
};

static const double g_6[6] = { 
   0.03522629188570953660274066472,
   0.08544127388202666169281916918,
  -0.13501102001025458869638990670,
  -0.45987750211849157009515194215,
   0.80689150931109257649449360409,
  -0.33267055295008261599851158914
};

static const double unity[1] = {1.0};

extern void my_kernel_wrapper(dim3          dimGrid, 
			      dim3          dimBlock, 
			      double       *tcIn,     const int  tcLength,
			      const int     voxelNb,
			      double       *num,      const int  lnum,
			      const int     lden,
			      double       *den1,     const int  lden1,
			      double       *den2,     const int  lden2,
			      const double lambdaTemp,
			      double        maxeig,
			      const int     cost_save,
			      const int     nit,
			      const double *noiseFinIn,
			      double       *tcOut,
			      double       *noiseFinOut);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double      *tcin,   *num,   *den1,   *den2,   *tcout,   *noiseFinIn,   *noiseFinOut;
  double      *d_tcin, *d_num, *d_den1, *d_den2, *d_tcout, *d_noiseFinIn, *d_noiseFinOut;
  double       lambdaTemp, /*lambda,*/ maxeig, noise_est;
  int          lden, cost_save, nit;
  const double *g6, *h6, *p_unity;

  g6      = g_6;
  h6      = h_6;
  p_unity = unity;

  tcin       = mxGetPr(prhs[0]);           /* Input data                     */
  num        = mxGetPr(prhs[1]);           /* Numerator of filter function   */
  lden       = (int) mxGetScalar(prhs[2]);
  den1       = mxGetPr(prhs[3]);           /* Causal part of denominator     */
  den2       = mxGetPr(prhs[4]);           /* Non-causal part of denominator */
  lambdaTemp = mxGetScalar(prhs[5]);
  maxeig     = mxGetScalar(prhs[6]);
  cost_save  = (int) mxGetScalar(prhs[7]);
  nit        = (int) mxGetScalar(prhs[8]);
  noiseFinIn = mxGetPr(prhs[9]);

  const mwSize *dimsTcin = mxGetDimensions(prhs[0]); /* tc length x number of voxels */
  const mwSize *dimsNum  = mxGetDimensions(prhs[1]); /* Numerator dimensions*/
  const mwSize *dimsDen1 = mxGetDimensions(prhs[3]); 
  const mwSize *dimsDen2 = mxGetDimensions(prhs[4]); 

  printf("tcin dimensions = %d x %d\n", dimsTcin[0], dimsTcin[1]);
  printf("num  dimensions = %d x %d\n", dimsNum[0],  dimsNum[1]);
  printf("lden            = %d\n",      lden);
  printf("den1 dimensions = %d x %d\n", dimsDen1[0], dimsDen1[1]);
  printf("den2 dimensions = %d x %d\n", dimsDen2[0], dimsDen2[1]);
  printf("lambdaTemp      = %15.10f\n", lambdaTemp);
  printf("maxeig          = %15.10f\n", maxeig);
  printf("cost_save       = %d\n",      cost_save);
  printf("numb. it. nit   = %d\n",      nit);

  // See if that should be flexible or not
  if (dimsNum[0] != 1 || dimsNum[1] != 6)
    mexErrMsgIdAndTxt("MyTemporal_MEX:num:dims", "Expected dims 1x6");

  uint tcinBytes     = dimsTcin[0] * dimsTcin[1] * sizeof(double);
  uint numBytes      = dimsNum[0]  * dimsNum[1]  * sizeof(double);
  uint den1Bytes     = dimsDen1[0] * dimsDen1[1] * sizeof(double);
  uint den2Bytes     = dimsDen2[0] * dimsDen2[1] * sizeof(double);
  uint noiseFinBytes = dimsTcin[1]               * sizeof(double);

  // Set up output for Matlab parallel process
  plhs[0]     = mxCreateDoubleMatrix(dimsTcin[0], dimsTcin[1], mxREAL);
  tcout       = mxGetPr(plhs[0]);
  plhs[1]     = mxCreateDoubleMatrix(dimsTcin[1], 1, mxREAL);
  noiseFinOut = mxGetPr(plhs[1]);

  // Set device
  int dev = 0;
  int deviceCount = 0;
  CHECK(cudaSetDevice(dev));
  CHECK(cudaGetDeviceCount(&deviceCount));
  printf("There are %d GPUs available.\n", deviceCount);

  size_t size;
  CHECK(cudaDeviceGetLimit(&size, cudaLimitMallocHeapSize));
  printf("GPU cudaLimitMallocHeapSize = %d\n", size);
  CHECK(cudaDeviceSetLimit(cudaLimitMallocHeapSize, size*50));
  CHECK(cudaDeviceGetLimit(&size, cudaLimitMallocHeapSize));
  printf("GPU cudaLimitMallocHeapSize = %d\n", size);


  // Force the creation of the CUDA context so context
  // creation overhead is readily distinguishable when profiling
  //CHECK(cudaFree(0));

  // Allocate device global memory
  CHECK(cudaMalloc((double**)&d_tcin,        tcinBytes));
  CHECK(cudaMalloc((double**)&d_num,         numBytes));
  CHECK(cudaMalloc((double**)&d_den1,        den1Bytes));
  CHECK(cudaMalloc((double**)&d_den2,        den2Bytes));
  CHECK(cudaMalloc((double**)&d_noiseFinIn,  noiseFinBytes));
  CHECK(cudaMalloc((double**)&d_tcout,       tcinBytes));
  CHECK(cudaMalloc((double**)&d_noiseFinOut, noiseFinBytes));


  // Transfer data from host to device
  CHECK(cudaMemcpy(d_tcin,     tcin,     tcinBytes,     cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(d_num,      num,      numBytes,      cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(d_den1,     den1,     den1Bytes,     cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(d_den2,     den2,     den2Bytes,     cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(d_noiseFinIn, noiseFinIn, noiseFinBytes, cudaMemcpyHostToDevice));

  //CHECK(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));
  CHECK(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));

  // Transfer common constant information in constant memory:
  //   - Daubechies D6 coefficients
  //   - Numerator and denominator filter coefficients
  set_constant_memory(g6, h6, p_unity);


  // Call CUDA kernel
  int  blockSize = 5;
  int  numBlocks = (dimsTcin[1] + blockSize -1) / blockSize;
  printf("blocksize is %d and numBlocks is %d\n", blockSize, numBlocks);
  dim3 dimGrid(numBlocks, 1);
  dim3 dimBlock(blockSize, 1);

  my_kernel_wrapper(dimGrid, dimBlock,
		    d_tcin, dimsTcin[0], dimsTcin[1],
		    d_num,  dimsNum[1],
		    lden,
		    d_den1, dimsDen1[1],
		    d_den2, dimsDen2[1],
		    lambdaTemp,
		    maxeig,
		    cost_save,
		    nit,
		    d_noiseFinIn,
		    d_tcout,
		    d_noiseFinOut);

  // Wait for GPU to finish before accessing on host
  CHECK(cudaDeviceSynchronize());
  
  CHECK(cudaMemcpy(tcout,       d_tcout,       tcinBytes,     cudaMemcpyDeviceToHost));
  CHECK(cudaMemcpy(noiseFinOut, d_noiseFinOut, noiseFinBytes, cudaMemcpyDeviceToHost));


  CHECK(cudaFree(d_tcin));
  CHECK(cudaFree(d_num));
  CHECK(cudaFree(d_den1));
  CHECK(cudaFree(d_den2));
  CHECK(cudaFree(d_tcout));  
  CHECK(cudaFree(d_noiseFinIn));
  CHECK(cudaFree(d_noiseFinOut));

  // Reset device
  CHECK(cudaDeviceReset());

  return;
}
