#include <stdlib.h>
#include <stdio.h>
#define nc 6


double* daubechies_6d(double *a,  size_t n) {

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

  int i, ii, j, jf, k, nmod, n1, ni, nh;
  

  double *out = (double *) malloc(n * sizeof(double));
  if (out == NULL) {
    fprintf (stderr, "calloc failed\n");
    exit(0);
  }

  double *scratch = (double *) calloc(sizeof(double), n);
  if (scratch == NULL) {
    fprintf (stderr, "calloc failed\n");
    exit(0);
  }

  /* Init to input */
  for (i=0; i<n; i++)
    out[i] = a[i];

  for (j = n; j >= 2; j >>= 1) {

    nmod = nc * j;
    n1   = j - 1;
    nh   = j >> 1;

    for (i=0; i<n; i++)
      scratch[i] = 0.0;

    for (ii = 0, i = 0; i < j; i += 2, ii++) {

      double h = 0, g = 0;
      ni = i + nmod;

      for (k = 0; k < nc; k++) {
	jf = n1 & (ni + k);
	h += h_6[k] * out[jf];
	g += g_6[k] * out[jf];
      }
      
      scratch[ii]      += h;
      scratch[ii + nh] += g;
    }

    for (i = 0; i < j; i++) {
      out[i] = scratch[i];
    }
  }

  free(scratch);

  return out;
}
