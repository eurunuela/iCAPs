#include <device_launch_parameters.h>
#include <device_functions.h>
#include <stdio.h>


__global__ void my_kernel(    double *tcin,     const int  tcLength, const int  voxelNb,
			      double *num,      const int  lnum,
			      const int lden,
			      double *den1,     const int  lden1,
			      double *den2,     const int  lden2,
			      double lambdaTemp,
			      double maxeig,
			      const int cost_save,
			      const int nit,
			      const double *noiseFinIn,
			      double *tcout,
			      double *noiseFinOut);

extern void my_kernel_wrapper(dim3   dimGrid,
			      dim3   dimBlock,
			      double *tcin,     const int  tcLength, const int  voxelNb,
			      double *num,      const int  lnum,
			      const int lden,
			      double *den1,     const int  lden1,
			      double *den2,     const int  lden2,
			      double lambdaTemp,
			      double maxeig,
			      const int cost_save,
			      const int nit,
			      const double *noiseFinIn,
			      double *tcout,
			      double *noiseFinOut) {
  my_kernel<<<dimGrid, dimBlock>>>(tcin, tcLength, voxelNb,
				   num,  lnum,
				   lden,
				   den1, lden1,
				   den2, lden2,
				   lambdaTemp,
				   maxeig,
				   cost_save,
				   nit,
				   noiseFinIn,
				   tcout,
				   noiseFinOut);

  return;
}
