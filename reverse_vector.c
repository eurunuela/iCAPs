/*

  Purpose: Reverse order of a vector of type double

*/
#include <stdlib.h>

double* reverse_vector(double *data,  int size) {

  size_t i, j;
  double tmp;

  for(i=0; i<(size/2); i++) {
    j   = size - i - 1;
    tmp = data[j];
    data[j] = data[i];
    data[i] = tmp;
  }
}
