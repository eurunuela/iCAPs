/*
  Author : E. Orliac, SCITAS, EPFL
  Date   : 27.10.2017
  Purpose: 
  Remarks: 
*/

/*
  num      : numerator
  lnum     : length of numerator
  den1     : causal
  den2     : non-causal
  ld1      : length causal
  ld2      : length non-causal
  condition: 0 = 'normal'
             1 = 'transpose' 
 */

#include <stdio.h>
#include <stdlib.h>
#include "filter_MEX.h"

void filter_MEX(double *out, double *in, int lin, double *num, int lnum, double *den, int lden) {

  int    i, j;
  double acc;

  for(i=0; i<lin; i++)
    out[i] = 0.0;

  if (lden == 1 && den[0] == 1.d) {
    for(i=0; i<lin; i++) {
      acc = 0.0d;
      for(j=0; j<lnum; j++) {
	if (j<=i) {
	  acc += num[j]*in[i-j];
	}
      }
      out[i] = acc;
    }
  } else if (lnum == 1 && num[0] == 1.d) {
    for(i=0; i<lin; i++) {
      acc = in[i];
      for(j=0; j<lden; j++) {
	if (j<=i) {
	  acc -= den[j]*out[i-j];
	}
      }
      out[i] = acc;
    }
  } else {
    fprintf (stderr, "Fatal. Unknown case in filter_MEX.\n");
    exit (EXIT_FAILURE);
  }
}
