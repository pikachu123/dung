#include "counts.h"
#include "population.h"
#include "binary.h"

// ==============================================================================

int computeCounts(int *pos, int n, Population *P, long *count)
{
  long i,N;
  long numConfigurations;
  long configuration;
  char *x;

  // initialize the variables

  N = P->N;
  numConfigurations = (long) 1<<n;
  
  // set counts to 0

  for (i=0; i<numConfigurations; i++)
    count[i]=0;

  // go through all strings and do the job

  for (i=0; i<N; i++)
    {
      x = P->individual[i].chromosome;
//       x = getChromosome(&(P->individual[i]),d);
      
      configuration = indexedBinaryToInt(x,pos,n);

      count[configuration]++;
    };

//   printf("Counts for ");
//   for (i=0; i<n; i++)
//     printf("%u ",pos[i]);
//   printf(" : ");
//   for (i=0; i<numConfigurations; i++)
//     printf("%u ",count[i]);
 
  // get back

  return 0;
};

// ==============================================================================

int computeCountsForList(int node, int *list, int numList, int *pos, int n, Population *P, long **count)
{
  register long i;
  long N;
  long numConfigurations;
  long numConfigurations4;
  long configuration;
  char *x;
  int l;

  // initialize the variables

  N = P->N;
  numConfigurations = (long) 1<<n;
  numConfigurations4 = (long) numConfigurations<<2;
  
  // set counts to 0

  for (i=0; i<numConfigurations4; i++)
    for (l=0; l<numList; l++)
      count[l][i]=0;

  // go through all strings and do the job

  for (i=0; i<N; i++)
    {
      x = P->individual[i].chromosome;
//       x = getChromosome(&(P->individual[i]),d);
      
      configuration   = (long) indexedBinaryToInt(x,pos,n);
      configuration <<= 2;
      configuration  += x[node];
      
      for (l=0; l<numList; l++)
	count[l][configuration+(x[list[l]]<<1)]++;
    };
  
  // get back
  
  return 0;
};
