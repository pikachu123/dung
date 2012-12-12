#include "BDe.h"

#include "graph.h"
#include "mymath.h"
#include "memalloc.h"
#include "binary.h"
#include "counts.h"
#include "random.h"
#include "recombination.h"

//#define DEBUG

// ==============================================================================

double logS(int i, Population *P, AcyclicOrientedGraph *G)
{
  double result;
  int numIn;
  int numParentConfigurations;
  int k;
  long j;
  int *parentIdx;
  int n;
  long N,l;
  long *tNi,**Ni;
  long *tNPi,**NPi;
  char *x=NULL;
  int parentConfiguration;

  // initialize the variables

  n                       = G->size();
  N                       = P->N;
  numIn                   = G->getNumIn(i);
  numParentConfigurations = (long) 1<<numIn;
  
  // allocate the memory for the index of parents of the node i

  if (numIn>0)
      parentIdx = (int*) Calloc(numIn,sizeof(int));
  else
      parentIdx = NULL;

  tNi  = (long*) Calloc(numParentConfigurations,sizeof(long));
  tNPi = (long*) Calloc(numParentConfigurations,sizeof(long));
  Ni   = (long**) Calloc(numParentConfigurations,sizeof(long*));
  NPi  = (long**) Calloc(numParentConfigurations,sizeof(long*));  

  for (j=0; j<numParentConfigurations; j++)
    {
      Ni[j] = (long*) Calloc(2,sizeof(long));
      NPi[j] = (long*) Calloc(2,sizeof(long));
    };

  // set parent index

  k=0;
  for (j=0; j<n; j++)          // optimize!!! (through the graph library..the list of parents..)
    if ((j!=i)&&(G->connected(j,i)))
      parentIdx[k++]=j;

  // calculate the N and NP (N-prime) coefficients

  for (j=0; j<numParentConfigurations; j++)
    {
      tNi[j]  = 0;
      tNPi[j] = 0;
      
      for (k=0; k<2; k++)
	{
	  Ni[j][k]  = 0;
	  NPi[j][k] = 0;
	}
    };

  // compute the N's

  for (l=0; l<N; l++)
    {
      x = P->individual[l].chromosome;
      
      parentConfiguration = indexedBinaryToInt(x,parentIdx,numIn);
      
#ifdef DEBUG      
      if (parentConfiguration>=numParentConfigurations)
	printf("ERROR: parent configuration out of range!\n");
#endif
      
      Ni[parentConfiguration][x[i]]++;
    };

  // compute the N-prime's (a simple uniformative K2 metric)

  for (j=0; j<numParentConfigurations; j++)
    for (k=0; k<2; k++)
      NPi[j][k]=1;

  // compute the N and N-prime totals
  
  for (j=0; j<numParentConfigurations; j++)
    for (k=0; k<2; k++)
      {
	tNi[j]  += Ni[j][k];
	tNPi[j] += NPi[j][k];
      };

  // compute the result

  result = 0;

  for (j=0; j<numParentConfigurations; j++)
    {
      //    result += sumLog(1,tNPi[j]);        // optimize!!! these two lines!!!! (unnecessary stuff here...)
      //    result -= sumLog(1,tNPi[j]+tNi[j]);
      
      result += -sumLog(tNPi[j]+1,tNPi[j]+tNi[j]); // better but can be still slightly optimize...+1 can be precomputed

      for (k=0; k<2; k++)
	{
	  //	  result += sumLog(1,NPi[j][k]+Ni[j][k]);  // optimize!!! these two lines!!! (unnecessary stuff here...)
	  //	  result -= sumLog(1,NPi[j][k]);

	  result += sumLog(NPi[j][k]+1,NPi[j][k]+Ni[j][k]);
	};
    };
  

  // free the allocated memory

  if (numIn>0)
    Free(parentIdx);

  Free(tNi);
  Free(tNPi);

  for (j=0; j<numParentConfigurations; j++)
    {
      Free(Ni[j]);
      Free(NPi[j]);
    };
  
  Free(Ni);
  Free(NPi);

  // get back with the result

  return result;
};

// ==============================================================================

double logGain(int i, int *parents, int n, Population *P)
{
  int k;
  long j;
  double gain;
  long  *count;
  long *tNwi,**Nwi;
  long *tNPwi,**NPwi;
  long *tNoi,**Noi;
  long *tNPoi,**NPoi;
  long numParentConfigurations;
  long numParentConfigurations2;
  long numParentConfigurations4;

  // initialize some variables
  
  numParentConfigurations  = (long) 1<<n;
  numParentConfigurations2 = (long) numParentConfigurations<<1;
  numParentConfigurations4 = (long) numParentConfigurations2<<1;

  // allocate memory for count

  count = (long*) Calloc(numParentConfigurations4,sizeof(long));

  // allocate memory for N's and so

  tNoi  = (long*) Calloc(numParentConfigurations,sizeof(long));
  tNPoi = (long*) Calloc(numParentConfigurations,sizeof(long));
  Noi   = (long**) Calloc(numParentConfigurations,sizeof(long*));
  NPoi  = (long**) Calloc(numParentConfigurations,sizeof(long*));  

  for (j=0; j<numParentConfigurations; j++)
    {
      Noi[j] = (long*) Calloc(2,sizeof(long));
      NPoi[j] = (long*) Calloc(2,sizeof(long));
    };

  tNwi  = (long*) Calloc(numParentConfigurations2,sizeof(long));
  tNPwi = (long*) Calloc(numParentConfigurations2,sizeof(long));
  Nwi   = (long**) Calloc(numParentConfigurations2,sizeof(long*));
  NPwi  = (long**) Calloc(numParentConfigurations2,sizeof(long*));  

  for (j=0; j<numParentConfigurations2; j++)
    {
      Nwi[j] = (long*) Calloc(2,sizeof(long));
      NPwi[j] = (long*) Calloc(2,sizeof(long));
    };

  // add i to the list of positions to compute frequencies for -------

  parents[n+1]=i;

  // compute the counts -----------------------------------------------

  computeCounts(parents,n+2,P,count);

#ifdef DEBUG  
  printPopulation(stdout,P);
  printf("Positions: ");
  for (j=0; j<n+2; j++)
    printf("%u ",parents[j]);
  printf("\n");
  for (j=0; j<1<<(n+2); j++)
    {
      char y[n+2];
      
      printf("Count for ");
      longToBinary(j,y,n+2);
      printBinary(stdout,y,n+2);
      printf(" is %lu\n",count[j]);
      }
#endif

  // compute the things we need ---------------------------------------

  // compute the N's

  for (j=0; j<numParentConfigurations2; j++)
    {
      tNwi[j]=0;
      tNPwi[j]=0;

      for (k=0; k<2; k++)
	{
	  Nwi[j][k] = count[(j<<1)+k];
	  NPwi[j][k]=1;

	  tNwi[j]  += Nwi[j][k];
	  tNPwi[j] += NPwi[j][k];

	  ///printf("Nwi[%u][%u] = %lu\n",j,k,Nwi[j][k]);
	}
    };
			       
  for (j=0; j<numParentConfigurations; j++)
    {
      tNoi[j]=0;
      tNPoi[j]=0;

      for (k=0; k<2; k++)
	{
	  Noi[j][k] = Nwi[j<<1][k]+Nwi[(j<<1)+1][k];
	  NPoi[j][k]=1;

	  tNoi[j]  += Noi[j][k];
	  tNPoi[j] += NPoi[j][k];

	  ///printf("Noi[%u][%u] = %lu\n",j,k,Noi[j][k]);
	}
    };

  // compute the gain -------------------------------------------------

  gain = 0;

  for (j=0; j<numParentConfigurations; j++)
    {
      gain += getPrecomputedCummulativeLog(tNPoi[j],tNPoi[j]+tNoi[j]-1);

      for (k=0; k<2; k++)
	gain -= getPrecomputedCummulativeLog(NPoi[j][k],NPoi[j][k]+Noi[j][k]-1);
    };

  for (j=0; j<numParentConfigurations2; j++)
    {
      gain -= getPrecomputedCummulativeLog(tNPwi[j],tNPwi[j]+tNwi[j]-1);

      for (k=0; k<2; k++)
	gain += getPrecomputedCummulativeLog(NPwi[j][k],NPwi[j][k]+Nwi[j][k]-1);
    };

#ifdef DEBUG 

  printf("Gain for (%u,%u)...%f\n",parents[n+1],parents[n],gain);
  getchar();

#endif

  // free the memory --------------------------------------------------

  Free(count);

  for (j=0; j<numParentConfigurations; j++)
    {
      Free(Noi[j]);
      Free(NPoi[j]);
    };
  
  Free(Noi);
  Free(NPoi);
  Free(tNoi);
  Free(tNPoi);

  for (j=0; j<numParentConfigurations2; j++)
    {
      Free(Nwi[j]);
      Free(NPwi[j]);
    };
  
  Free(Nwi);
  Free(NPwi);  
  Free(tNwi);
  Free(tNPwi);
  
  // return the gain --------------------------------------------------

  return gain;
}

// ==============================================================================

int K2ComputeLogGains( int i, 
		       double oldContribution,
		       double **gain, 
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
  long  *tNwi,**Nwi;
  long  *tNPwi,**NPwi;
  long  *tNoi,**Noi;
  long  *tNPoi,**NPoi;
  long  numParentConfigurations;
  long  numParentConfigurations2;
  long  numParentConfigurations4;
  double result;

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

  // allocate memory for N's and so (need only one of these)

  tNoi  = (long*) Calloc(numParentConfigurations,sizeof(long));
  tNPoi = (long*) Calloc(numParentConfigurations,sizeof(long));
  Noi   = (long**) Calloc(numParentConfigurations,sizeof(long*));
  NPoi  = (long**) Calloc(numParentConfigurations,sizeof(long*));  

  for (j=0; j<numParentConfigurations; j++)
    {
      Noi[j] = (long*) Calloc(2,sizeof(long));
      NPoi[j] = (long*) Calloc(2,sizeof(long));
    };

  tNwi  = (long*) Calloc(numParentConfigurations2,sizeof(long));
  tNPwi = (long*) Calloc(numParentConfigurations2,sizeof(long));
  Nwi   = (long**) Calloc(numParentConfigurations2,sizeof(long*));
  NPwi  = (long**) Calloc(numParentConfigurations2,sizeof(long*));  

  for (j=0; j<numParentConfigurations2; j++)
    {
      Nwi[j] = (long*) Calloc(2,sizeof(long));
      NPwi[j] = (long*) Calloc(2,sizeof(long));
    };

  // compute the counts -----------------------------------------------

  computeCountsForList(i,updateIdx,numUpdated,parentList,numParents,P,count);

  // for each element of the nodes to be updated update the gain

  for (l=0; l<numUpdated; l++)
    {      
      ///      printf("Computing gain for (%u,%u)\n",updateIdx[l],i);

      // compute the N's
      
      for (j=0; j<numParentConfigurations2; j++)
	{
	  tNwi[j]=0;
	  tNPwi[j]=0;
	  
	  for (k=0; k<2; k++)
	    {
	      Nwi[j][k] = count[l][(j<<1)+k];
	      NPwi[j][k]=1;
	      
	      tNwi[j]  += Nwi[j][k];
	      tNPwi[j] += NPwi[j][k];

	      ///	      printf("Nwi[%u][%u] = %lu\n",j,k,Nwi[j][k]);
	    }
	};

      for (j=0; j<numParentConfigurations; j++)
	{
	  tNoi[j]=0;
	  tNPoi[j]=0;
	  
	  for (k=0; k<2; k++)
	    {
	      Noi[j][k] = Nwi[j<<1][k]+Nwi[(j<<1)+1][k];
	      NPoi[j][k]=1;
	      
	      tNoi[j]  += Noi[j][k];
	      tNPoi[j] += NPoi[j][k];

	      ///	      printf("Noi[%u][%u] = %lu\n",j,k,Noi[j][k]);
	    }
	};

      ///            getchar();
      
      // compute the resulting gain for the addition of an edge from updateIdx[l] to i

      result = 0;
      
      for (j=0; j<numParentConfigurations; j++)
	{
	  result += getPrecomputedCummulativeLog(tNPoi[j]+1,tNPoi[j]+tNoi[j]);
	  
	  for (k=0; k<2; k++)
	    result -= getPrecomputedCummulativeLog(NPoi[j][k]+1,NPoi[j][k]+Noi[j][k]);
	};
      
      for (j=0; j<numParentConfigurations2; j++)
	{
	  result -= getPrecomputedCummulativeLog(tNPwi[j]+1,tNPwi[j]+tNwi[j]);
	  
	  for (k=0; k<2; k++)
	    result += getPrecomputedCummulativeLog(NPwi[j][k]+1,NPwi[j][k]+Nwi[j][k]);
	};
      
      // update the gain

      //printf("Gain for (%u,%u)...%f\n",updateIdx[l],i,result);
      ///      getchar();

      gain[updateIdx[l]][i]=result;
    };

  // free the memory

  for (l=0; l<numUpdated; l++)
    Free(count[l]);

  if (numUpdated>0)
      Free(count);

  for (j=0; j<numParentConfigurations; j++)
    {
      Free(Noi[j]);
      Free(NPoi[j]);
    };
  
  Free(Noi);
  Free(NPoi);
  Free(tNoi);
  Free(tNPoi);

  for (j=0; j<numParentConfigurations2; j++)
    {
      Free(Nwi[j]);
      Free(NPwi[j]);
    };
  
  Free(Nwi);
  Free(NPwi);  
  Free(tNwi);
  Free(tNPwi);
  
  // get back

  return 0;
}

// =============================================================

double K2ComputeIsolatedNodeContribution(int i, Population *P)
{
  // !!!!!!!!!!!! not done !!!!!!!!!!!!!

  return 0;
};
