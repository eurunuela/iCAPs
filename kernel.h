#include "cuda_runtime.h"
#ifndef KERNEL_H
#define KERNEL_H

// Set constant memory (Daubechies D6 coefficients)
void set_constant_memory(const double *g6,
			 const double *h6,
			 const double *one);

#endif
