#include "population.h"
#include "replace.h"
#include "random.h"
#include "select.h"
#include "distance.h"

#include <math.h>
#include <stdlib.h>

char *replacementDesc[6] = { "Replace Worst",
			     "Replace Any",
			     "lambda,mju",
			     "lambda+mju",
			     "RTS",
			     "Crowding"
};

// ============================================================

int replaceWorst(Population *population, Population *offspring)
{
  long i,j;
  long M,N,NM;
  int numDiscrete,numContinuous;

  // initialize variables

  N  = population->N;
  M  = offspring->N;
  numDiscrete  = population->numDiscrete;
  numContinuous = population->numContinuous;

  NM = N-M;
      
  // shuffle the individuals a little
  
  for (i=0; i<N; i++)
    {
      j = (long) ((double) drand()*N);
      
      if (i!=j)
	swapIndividuals(population,i,j);
    };
  
  // separate worst M individuals by divide and conquer 
  
  divideWorst(population,0,N-1,numDiscrete,numContinuous,N-M);
  
  /*
    for (i=0; i<N; i++)
    {
    if (i==NM)
    printf("------------------\n");
    
    printf("fitness[%3lu]= %1.2f\n",i,population->individual[i].f);
    };
    getchar();
  */
  
  // replace the worst M
  
  for (i=NM; i<N; i++)
    replaceIndividual(population,i,&(offspring->individual[i-NM]));
  
  return 0;
}

// ============================================================

int replaceAny(Population *population, Population *offspring)
{
  long i,j;
  long M,N;
  int n;

  // initialize variables

  N  = population->N;
  M  = offspring->N;
  n  = population->n;
  
  // shuffle the individuals a little
  
  for (i=0; i<N; i++)
    {
      j = (long) ((double) drand()*N);
      
      if (i!=j)
	swapIndividuals(population,i,j);
    };
  
  // replace the first individuals
  
  for (i=0; i<M; i++)
    replaceIndividual(population,i,&(offspring->individual[i]));
  
  // get back
  
  return 0;
}

// ===============================================================

int lambdaCommaMju(Population *population, Population *offspring)
{
  int numDiscrete,numContinuous;

  // assign some helper variables

  numDiscrete  = population->numDiscrete;
  numContinuous = population->numContinuous;

  // get the best from offspring

  divideBest(offspring, 0, offspring->N-1, numDiscrete, numContinuous, population->N);

  // put them into the original population

  for (long i=0; i<population->N; i++)
    copyIndividual(population->individual+i,offspring->individual+i,numDiscrete, numContinuous);

  // get back

  return 0;
};

// ===============================================================

int lambdaPlusMju(Population *population, Population *offspring)
{
  long i;
  Population big;
  int numDiscrete,numContinuous;

  // assign some helper variables

  numDiscrete  = population->numDiscrete;
  numContinuous = population->numContinuous;

  // create a big population to put the both in

  allocatePopulation(&big,offspring->N+population->N,numDiscrete,numContinuous);

  // copy them all inside

  for (i=0; i<population->N; i++)
    copyIndividual(big.individual+i,population->individual+i,numDiscrete,numContinuous);

  for (i=0; i<offspring->N; i++)
    copyIndividual(big.individual+i+population->N,offspring->individual+i,numDiscrete,numContinuous);

  // get the best from the big population

  divideBest(&big, 0, big.N-1, numDiscrete,numContinuous, population->N);

  // put them into the original population

  for (i=0; i<population->N; i++)
    copyIndividual(population->individual+i,big.individual+i,numDiscrete,numContinuous);

  // get back

  return 0;
};

// ===============================================================

char *getReplacementDesc(int n)
{
  return replacementDesc[n];
};

//=====================================================
// do divide and conquer until the [NM..N] individuals
// are worst or equal to the rest (NM = N-M, last M 
// should be worst)

int divideWorst(Population *population, long left, long right, int numDiscrete, int numContinuous, long NM)
{
  long l,r;
  double pivot;

  l  = left;
  r  = right;

  pivot = (population->individual[l].f+population->individual[r].f)/2;

  while (l<=r)
    {
      while ((l<right)&&(population->individual[l].f>pivot)) l++;
      while ((r>left)&&(population->individual[r].f<pivot)) r--;

      if (l<=r)
	{
	  if (l!=r)
	    swapIndividuals(population,l,r);

	  l++;
	  r--;
	}
    };

  if ((l==NM)||(r==(NM-1)))
    return 0;
  
  if ((r>=NM)&&(left<r))
    divideWorst(population,left,r,numDiscrete,numContinuous,NM);

  if ((l<NM)&&(l<right))
    divideWorst(population,l,right,numDiscrete,numContinuous,NM);
  
  return 0;
};

// ============================================================

int restrictedTournament(Population *population, Population *offspring)
{
  long i,j;
  long M,N;
  int numDiscrete,numContinuous;
  double dist=0;
  int windowSize;

  // initialize variables

  N  = population->N;
  M  = offspring->N;
  numDiscrete  = population->numDiscrete;
  numContinuous = population->numContinuous;

  windowSize = (numDiscrete>numContinuous)? numDiscrete:(numContinuous*5);

  if ((numDiscrete==0)&&(numContinuous==0))
    {
      fprintf(stderr,"ERROR: Restricted tournament can't work with no variables whatsoever.\n");
      exit(-1);
    };
  
  // shuffle the individuals a little
  
  for (i=0; i<N; i++)
    {
      j = (long) ((double) drand()*N);
      
      if (i!=j)
	swapIndividuals(population,i,j);
    };

  // replace the closest if better

  for (i=0; i<M; i++)
    { 
      double min=0;
      long which;

      which = longRand(N);

      if (numContinuous==0)
	min   = hammingDistance(&(population->individual[which]),&(offspring->individual[i]),numDiscrete);
      else
	if (numDiscrete==0)
	  min = euclid(population->individual[which].continuous,population->individual[i].continuous,numContinuous);
      
      for (long jj=1; jj<windowSize; jj++)
	{
	  j=longRand(N);
	  
	  if (numContinuous==0)
	    dist = hammingDistance(&(population->individual[j]),&(offspring->individual[i]),numDiscrete);
	  else
	    if (numDiscrete==0)
	      min = euclid(population->individual[j].continuous,population->individual[i].continuous,numContinuous);
	  if (dist<min)
	    {
	      which = j;
	      min=dist;
	    };
	};
      
      if (population->individual[which]<offspring->individual[i])
	replaceIndividual(population,which,&(offspring->individual[i]));
    };

  // get back

  return 0;
}

// ============================================================

int crowding(Population *population, Population *offspring)
{
  long i,j;
  long M,N;
  int n;
  double dist;

  // initialize variables

  N  = population->N;
  M  = offspring->N;
  n  = population->n;
  
  // shuffle the individuals a little
  
  for (i=0; i<N; i++)
    {
      j = (long) ((double) drand()*N);
      
      if (i!=j)
	swapIndividuals(population,i,j);
    };

  // replace the closest

  for (i=0; i<M; i++)
    { 
      double min;
      long which;

      which = longRand(N);
      min   = hammingDistance(&(population->individual[which]),&(offspring->individual[i]),n);
      
      for (long jj=1; jj<n; jj++)
	{
	  j=longRand(N);
	  
	  dist = hammingDistance(&(population->individual[j]),&(offspring->individual[i]),n);
	  
	  if (dist<min)
	    {
	      which = j;
	      min=dist;
	    };
	};
      
      replaceIndividual(population,which,&(offspring->individual[i]));
    };

  // get back

  return 0;
}
