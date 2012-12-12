#include "mymath.h"
#include "memalloc.h"

#define invLog2 1.44269504088896338700465094007086008787155151367188
#define invLog3 0.9102392266268373166582250632927753031254
#define MIN_XLOG 1E-50

double *precomputedCummulativeLogarithm;
long   *precomputedCombine;

//======================================================================

//  float fabs(float x)
//  {
//    if (x<0)
//      return -x;
//    else
//      return x;
//  }

//======================================================================

//  long labs(long x)
//  {
//    if (x<0)
//      return -x;
//    else
//      return x;
//  }

//======================================================================

float fsqr(float x)
{
  return x*x;
}

//======================================================================

double dsqr(double x)
{
  return x*x;
}

//======================================================================

long lsqr(long x)
{
  return x*x;
}

//======================================================================

float fmax(float a, float b)
{
  if (a>b)
    return a;
  else
    return b;
}

//======================================================================

int iabs(int x)
{
  if (x<0)
    return -x;
  else
    return x;
}

//======================================================================

int iPower(int x, int a)
{
  if (a==0)
    return 1;
  else
  if (a%2==0)
    return iPower(x*x,a>>1);
  else
    return x*iPower(x,a-1);
}

//======================================================================

double dFactorial(int n)
{
  if (n)
    return n*dFactorial(n-1);
  else
    return 1;
}

//======================================================================

double sumLog(long i, long j)
{
  double result;
  register long k;

  result = 0;

  for (k=i; k<=j; k++)
      result += log2(k);
  
  return result;
};

//======================================================================

double getPrecomputedCummulativeLog(long i, long j)
{
  return ((double) (precomputedCummulativeLogarithm[j]-precomputedCummulativeLogarithm[i-1]));
};

//======================================================================

int precomputeCummulativeLogarithms(long n)
{
  long i;
  double sum;
  
  precomputedCummulativeLogarithm = (double*) Calloc(n+1,sizeof(double));
  sum = 0;
  precomputedCummulativeLogarithm[0]=0;
  
  for (i=1; i<=n; i++)
    {
      sum += log2(i);
      precomputedCummulativeLogarithm[i] = sum;
    };

  return 0;
};

//======================================================================

int freePrecomputedCummulativeLogarithms()
{
  Free(precomputedCummulativeLogarithm);

  return 0;
};

//======================================================================

double log2(double x)
{
  if (x<=MIN_XLOG)
    x=MIN_XLOG;
  return invLog2*log(x);
};

//======================================================================

double log3(double x)
{
  return invLog3*log(x);
};
