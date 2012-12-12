#include <stdio.h>
#include <stdlib.h>

#include "memalloc.h"
#include "reordering.h"

#define numReordering 2

Reordering reordering[numReordering] = {
  {"Normal order",&normalOrdering},
  {"Disorder-k",&disorderK}
};

int *order;

// =======================================================

void *getReorderingOperator(int n)
{
    if (n>=numReordering)
	{
	    fprintf(stderr,"ERROR: Reordering operator not available!");
	    exit(-1);
	    return 0;
	}
    else
	return &(reordering[n]);
};

// =======================================================

char *getReorderingDescription(int n)
{
  return reordering[n].description;
};

// =======================================================

int reorder(char *y, char *x, int n)
{
  register int i;
  
  for (i=0; i<n; i++)
    {
      y[i] = x[order[i]];
      ///      printf("Order of %u bit is %u\n",i,order[i]);
    }

  return 0;
}

// =======================================================

int normalOrdering(int n, double *params)
{
  register int i;

//  printf("Normal ordering initialization entered with n=%i.\n",n);

  if (n>0)
    order = (int*) Calloc(n,sizeof(int));
  else
    order = NULL;

  for (i=0; i<n; i++)
    order[i]=i;

  return 0;
};

// =======================================================

int disorderK(int n, double *params)
{
  register int i;
  int k,nk;

  order = (int*) Calloc(n,sizeof(int));

  k = (int) params[0];

  if (n%k!=0)
    {
      fprintf(stderr,"ERROR: For disorderK ordering operator n has to be dividable by K!");
      exit(-1);
    }

  nk = n/k;

  for (i=0; i<n; i++)
    order[i] = (i%k)*nk+i/k;

  return 0;
};

// =======================================================

int doneReordering()
{
  if (order)
    Free(order);

  return 0;
};
