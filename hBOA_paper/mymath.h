#ifndef _mymath_h_
#define _mymath_h_

#include <math.h>

//float fabs(float x);
//long labs(long x);
float fsqr(float x);
double dsqr(double x);
long lsqr(long x);
float fmax(float a, float b);
int iabs(int x);
int iPower(int x, int a);  
double dFactorial(int n);
double sumLog(long i, long j);
double getPrecomputedCummulativeLog(long i, long j);
int precomputeCummulativeLogarithms(long n);
int freePrecomputedCummulativeLogarithms();
double log2(double x);
double log3(double x);

#endif
