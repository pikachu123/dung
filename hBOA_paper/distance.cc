#include <math.h>
#include <stdlib.h>

#include "mymath.h"

int hamming(char *x, char *y, int n)
{
  int result=0;
  
  for (int i=0; i<n; i++)
      result += abs(x[i]-y[i]);

  return result;
};

double euclidBinaryDouble(char *x, double *y, int n)
{
  double result;
  
  result=0;

  for (int i=0; i<n; i++)
    result += fsqr(x[i]-y[i]);

  return sqrt(result);
};

double euclid(double *x, double *y, int n)
{
  double result;
  
  result=0;

  for (int i=0; i<n; i++)
    result += fsqr(x[i]-y[i]);

  return sqrt(result);
};
