/*
  Author : E. Orliac, SCITAS, EPFL
  Date   : 27.10.2017
  Purpose: 
  Remarks: 
*/

/* 
   num      : numerator
   den1     : causal
   den2     : non-causal
   lnum     : length numerator
   lden1    : length causal
   lden2    : length non-causal
   condition: 0 = 'normal'
              1 = 'transpose' 
 */

#include <stdio.h>
#include <stdlib.h>
#include "filter_MEX.h"
#include "reverse_vector.h"

void filter_boundary_MEX(double *out,
			 double *in,  int lin, 
			 double *num, int lnum, 
			 int lden, 
			 double *den1, int lden1, 
			 double *den2, int lden2, 
			 int condition) {
  
  int    i, shiftnc;
  double unity[1] = {1.d};
  double *den;
  double *tmp = malloc(lin * sizeof(double));

  /* Check input condition */
  if (condition != 0 && condition != 1) {
    fprintf(stderr, "FATAL. Unknown input 'condition'. Allowed are 0 (normal) and 1 (traspose).\n");
    exit(EXIT_FAILURE);
  }

  /* Copy in in tmp */
  for (i=0; i<lin; i++)
    tmp[i] = in[i];

  /*printf("mex in out[0] = %.10f, out[end] = %.10f\n", out[0], out[lin-1]);*/
  
  if (condition == 1) {
    reverse_vector(tmp, lin);
    den = unity;
    filter_MEX(out, tmp, lin, num, lnum, den, 1);
    reverse_vector(out, lin);
  } else if (condition == 0) {
    den = unity;
    filter_MEX(out, tmp, lin, num, lnum, den, 1);
  } else {
    fprintf(stderr, "Fatal. Unknown case in filter_boundary_MEX.\n");
    exit(EXIT_FAILURE);
  }

  if (lden == 2) {
    /* den1 is causal; den2 is non-causal. */
    if (lden1+lden2 > 2) {
      shiftnc = lden2 - 1;
      if (shiftnc != 0) {
	fprintf(stderr, "FATAL. Non-zero shiftnc! Need to add missing implementation.\n");
        exit(EXIT_FAILURE); 
      }
      if (condition == 0) {
	num = unity;
	filter_MEX(tmp, out, lin, num, 1, den1, lden1);
	reverse_vector(tmp, lin);
	filter_MEX(out, tmp, lin, num, 1, den2, lden2);
	for (i=0; i<lin; i++)
	  out[i] *= den2[lden2-1];
	reverse_vector(out, lin);
	/*
	  TO IMPLEMENT IF shiftnc is not zero
	  out = out(2*shiftnc+1:end);
	*/
	/*printf("#normal# out1 out[0] = %.10f, out[end] = %.10f\n", out[0], out[lin-1]);*/
      } else {
	num = unity;
	reverse_vector(out, lin);
	filter_MEX(tmp, out, lin, num, 1, den1, lden1);
	reverse_vector(tmp, lin);
	filter_MEX(out, tmp, lin, num, 1, den2, lden2);
	for (i=0; i<lin; i++)
	  out[i] *= den2[lden2-1];
	/*
	  TO IMPLEMENT IF shiftnc is not zero
	  out = out(1:end-2*shiftnc);
	*/
	/*printf("#transp# out1 out[0] = %.10f, out[end] = %.10f\n", out[0], out[lin-1]);*/
      }
    }
  }

  free(tmp);
}
