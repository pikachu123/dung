#include "bic.h"

#include "counts.h"
#include "population.h"
#include "memalloc.h"
#include "mymath.h"
#include "utils.h"
#include "recombination.h"
#include "index.h"
#include "frequencyTree.h"

//#define DISPLAY_THE_INSTANCES
//#define DEBUG

int bicComputeAdditionGains( int i, 
			     double oldContribution,
			     double **additionGain, 
			     int *updateIdx, 
			     int numUpdated, 
			     int *parentList,
			     int numParents, 
			     Population *P,
			     RecombinationParams *params)
{
  long  j;
  int   k,l;
  long  **count;
  double **pwi,*tpwi;
  long  numParentConfigurations;
  long  numParentConfigurations2;
  long  numParentConfigurations4;
  double result;
  double structureDescriptionLength;
  double tableAndDataDescriptionLength;
  double entropyN;

#ifdef DEBUG
  printf("in the bic (additions)!\n");
#endif

  // initialize some variables
  
  numParentConfigurations  = (long) 1<<numParents;
  numParentConfigurations2 = (long) numParentConfigurations<<1;
  numParentConfigurations4 = (long) numParentConfigurations2<<1;

  // allocate memory for counts of all n-ths

  if (numUpdated>0)
    count = (long**) Calloc(numUpdated,sizeof(long*));
  else
    count = NULL;
  
  for (l=0; l<numUpdated; l++)
    count[l] = (long*) Calloc(numParentConfigurations4,sizeof(long));

  // allocate memory for p's and so (need only one of these)

  tpwi  = (double*)  Calloc(numParentConfigurations2,sizeof(double));
  pwi   = (double**) Calloc(numParentConfigurations2,sizeof(double*));
  
  for (j=0; j<numParentConfigurations2; j++)
    pwi[j] = (double*) Calloc(2,sizeof(double));
  
  // compute the counts -----------------------------------------------
  
  computeCountsForList(i,updateIdx,numUpdated,parentList,numParents,P,count);
  
  // compute the frequencies for the edge and its parents before the addition
  

  structureDescriptionLength = 0; // BIC doesn't care about the structure like MDL 
  ///  structureDescriptionLength  = log2(P->n)+getPrecomputedCummulativeLog(P->n-numParents,P->n)-getPrecomputedCummulativeLog(1,numParents+1);
  
  // for each element of the nodes to be updated update the gain
  
  for (l=0; l<numUpdated; l++)
    {      
      ///printf("Computing gain for (%u,%u)\n",updateIdx[l],i);
      
      // compute the N's
      
      for (j=0; j<numParentConfigurations2; j++)
	{
	  tpwi[j] = 0;
	  
	  for (k=0; k<2; k++)
	    {
	      pwi[j][k] = (double) (count[l][(j<<1)+k])/P->N;
	      tpwi[j]  += pwi[j][k];
	    };
	  
	  ///	  if (i==1)
	  ///	    printf("parent = %u, tpwi[%u] = %f\n",updateIdx[l],j,tpwi[j]);
	};
      

      // compute the new entropy
      
      entropyN = 0;
      
      for (j=0; j<numParentConfigurations2; j++)
	for (k=0; k<2; k++) 
	  if (pwi[j][k]>0)
	    entropyN -= pwi[j][k]*log2(pwi[j][k]/tpwi[j]);
      
      entropyN *= P->N;
      
      tableAndDataDescriptionLength = entropyN+double(numParentConfigurations)*log2(P->N);
      
      // compute the resulting gain for the addition of an edge from updateIdx[l] to i
      
      result = -oldContribution-structureDescriptionLength-tableAndDataDescriptionLength;
      
      // update the gain
       
      ///printf("Addition gain for (%u,%u)...%f (%f, %f, %f) with %u parents\n",updateIdx[l],i,result,oldContribution,structureDescriptionLength,tableAndDataDescriptionLength,numParents);
      ///getchar();
      
      additionGain[updateIdx[l]][i]=result;
    };
  
  // free the memory
  
  for (l=0; l<numUpdated; l++)
    Free(count[l]);
  
  if (numUpdated>0)
    Free(count);
  
  for (j=0; j<numParentConfigurations2; j++)
    Free(pwi[j]);
  
  Free(pwi);
  Free(tpwi);
  
  // get back
  
#ifdef DEBUG
  printf("Exiting bic...\n");
#endif

  return 0;
};

// ====================================================================

int bicComputeRemovalGains( int i, 
			    double oldContribution,
			    double **removalGain, 
			    int *parentList,
			    int numParents, 
			    Population *P,
			    RecombinationParams *params)
{
  long  j;
  int   k,l;
  long  *count;
  double **poi,*tpoi;
  long  numParentConfigurations;
  long  numParentConfigurations2;
  long  numParentConfigurationsD2;
  double result;
  double structureDescriptionLength;
  double tableAndDataDescriptionLength;
  double entropyN;
  long lExp;
  long lTmp;

  ///  printf("in the bic (removals)!\n");

  // any parents to remove?
  
  if (numParents==0)
    return 0;

  // initialize some variables
  
  numParentConfigurations   = (long) 1<<numParents;
  numParentConfigurations2  = (long) numParentConfigurations<<1;
  numParentConfigurationsD2 = (long) numParentConfigurations>>1;

  // allocate memory for counts of all n-ths

  count = (long*) Calloc(numParentConfigurations2,sizeof(long));

  // allocate memory for p's and so (need only one of these)

  tpoi  = (double*)  Calloc(numParentConfigurationsD2,sizeof(double));
  poi   = (double**) Calloc(numParentConfigurationsD2,sizeof(double*));

  for (j=0; j<numParentConfigurationsD2; j++)
    poi[j] = (double*) Calloc(2,sizeof(double));

  // add the current node into the list of parents to compute the counts (at the end)

  parentList[numParents]=i;

  // compute the counts -----------------------------------------------

  computeCounts(parentList,numParents+1,P,count);

  // compute the model description length decrease (the same for all) as 
  // a sum of the increase in graph size in bits + size of tables in bits

  structureDescriptionLength  = 0; // BIC doesn't care about the structure
  ///structureDescriptionLength  = log2(P->n)+getPrecomputedCummulativeLog(P->n-numParents+2,P->n)-getPrecomputedCummulativeLog(1,numParents-1);

  // for each element of the nodes to be removed update the gain

  for (l=0; l<numParents; l++)
    {      
      ///      printf("Compute removal gain for (%u, %u)\n",parentList[l],i);
      
      lExp = (long) 1<<(numParents-l);

      // compute the po's
      
      for (j=0; j<numParentConfigurationsD2; j++)
	{
	  tpoi[j] = 0;
	  
	  lTmp = (((j>>(numParents-l-1))<<(numParents-l))+(j&(((long)1<<(numParents-l-1))-1)))<<1;

	  ///	  printf("lExp = %lu, lTmp = %lu\n",lExp,lTmp);

	  for (k=0; k<2; k++)
	    {
	      poi[j][k] = (double) (count[lTmp+k]+count[lTmp+lExp+k])/P->N;
	      tpoi[j]  += poi[j][k];
	      ///	      if (i==1)
	      ///		printf("parent = %u, poi[%u][%u] = %f\n",parentList[l],j,k,poi[j][k]);
	    };
	};

      // compute the new entropy

      entropyN = 0;
      
      for (j=0; j<numParentConfigurationsD2; j++)
	for (k=0; k<2; k++) 
	  if (poi[j][k]>0)
	    entropyN -= poi[j][k]*log2(poi[j][k]/tpoi[j]);
      
      entropyN *= P->N;
      
      tableAndDataDescriptionLength = entropyN+double(0.25)*numParentConfigurations*log2(P->N);

      // compute the resulting gain for the removal of an edge from updateIdx[l] to i

      result = -oldContribution-structureDescriptionLength-tableAndDataDescriptionLength;
      
      // update the gain

      ///      printf("Removal gain for (%u,%u)...%f (%f + %f + %f)\n",parentList[l],i,result,oldContribution,structureDescriptionLength,tableAndDataDescriptionLength);
      ///      getchar();
     
      removalGain[parentList[l]][i]=result;
    };

  // free the memory

  Free(count);

  for (j=0; j<numParentConfigurationsD2; j++)
    Free(poi[j]);
  
  Free(poi);
  Free(tpoi);

  // get back

  ///  printf("Exiting bic...\n");

  return 0;
};

// =======================================================================================

double bicComputeIsolatedNodeContribution(int i, Population *P)
{
  double result;
  int    idx[1];
  long   counts[2];

  // compute the frequency

  idx[0]=i;
  computeCounts(idx,1,P,counts);

  // network complexity is simply log n

  result = log2(P->n);

  // table complexity is log N / 2

  result += log2(P->N)/2;

  // entropy is like this

  result -= counts[0]*log2(double(counts[0])/P->N);
  result -= counts[1]*log2(double(counts[1])/P->N);

  ///  printf("Contribution[%u] = %lf (%lf,%lf)\n",i,result,log2(P->n)+log2(P->N)/2,counts[0]*log2(double(counts[0])/P->N)+counts[1]*log2(double(counts[1])/P->N));

  // get back with the result

  return -result;  
};
