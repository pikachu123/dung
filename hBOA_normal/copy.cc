#include <stdio.h>
#include <stdlib.h>

#include "mymath.h"
#include "memalloc.h"
#include "random.h"
#include "population.h"
#include "recombination.h"
#include "copy.h"

// ===========================================================================

int copyRecombination(Population *parents, Population *children, long M, RecombinationParams *params)
{
  register long i;
  long N;
  int  n;

  // initialize the variables

  N = parents->N;
  n = parents->n;

  // size of parent population has to be equal to the size of offspring

  if (M!=N)
    {
      fprintf(stderr,"ERROR: With copy recombination the size of offspring has to be the same as the size of the set of selected parents.\n");
      exit(-1);
    }

  // allocate memory for offspring

  allocatePopulation(children,N,parents->numDiscrete,parents->numContinuous);

  // generate offspring successively

  for (i=0; i<N; i++)
    {
      // copy the parent into the offspring pop.

      copyIndividual(&(children->individual[i]),&(parents->individual[i]),parents->numDiscrete,parents->numContinuous);
    }

  children->evaluated=0;
  
  if (children->numContinuous>0) // do sigma-self-adaptive mutation
    {
/*       printf("Mutating children\n");getchar(); */
      double tau=4;
      double coef = exp(tau*gaussianRandom(0,1)/sqrt(children->numContinuous));
      
      for (i=0; i<N; i++)
	{
	  for (int k=0; k<children->numContinuous; k++)
	    {
	      children->individual[i].mutation[k]*=coef;
	      /* 	    printf("Before %f\n",children->individual[i].continuous[k]); */
	      children->individual[i].continuous[k] += gaussianRandom(0,children->individual[i].mutation[k]);
	      /* 	    printf("After %f\n",children->individual[i].continuous[k]); */
	      /* 	    getchar(); */
	      
	      if (children->individual[i].continuous[k]>1)
		children->individual[i].continuous[k]=1;
	      if (children->individual[i].continuous[k]<0)
		children->individual[i].continuous[k]=0;
	    };

	  children->individual[i].fCalculated=0;
	};
    };

  // get back

  return 0;
}
