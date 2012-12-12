#include <stdlib.h>

#include "arrays.h"

int CoupleArray::add(int a, int b)
{
  n++;
  if (x)
    x = (int*) realloc(x,sizeof(int)*n);
  else
    x = (int*) malloc(sizeof(int));

  if (y)
    y = (int*) realloc(y,sizeof(int)*n);
  else
    y = (int*) malloc(sizeof(int));

  x[n-1]=a;
  y[n-1]=b;

  return n;
}

CoupleArray::~CoupleArray()
{
  free(x);
  free(y);
}

CoupleArray::CoupleArray()
{
  x=y=NULL;
  clear();
}

void CoupleArray::clear()
{
  n=0;
  if (x!=NULL) free(x);
  if (y!=NULL) free(y);
  x=y=NULL;
}

int CoupleArray::getSize()
{
  return n;
}

int CoupleArray::get(int *a, int *b, int i)
{
  *a = x[i];
  *b = y[i];

  return 0;
}
