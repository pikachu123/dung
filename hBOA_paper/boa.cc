#include <stdio.h>
#include <string.h>

#include "boa.h"

#include "bic.h"
#include "mymath.h"
#include "memalloc.h"
#include "individual.h"
#include "population.h"
#include "recombination.h"
#include "graph.h"
#include "binary.h"
#include "utils.h"
#include "BDe.h"
//#include "mda.h"
#include "mdl.h"
#include "priors.h"
#include "random.h"
#include "counts.h"

//#define BOA_PRINT_PEARSON
//#define BOA_PRINT_PARENTS
//#define BOA_PRINT_COINCIDENCE
//#define BOA_PRINT_GROUPS 
//#define BOA_PRINT_MAXGAIN
//#define BOA_PRINT_TOTAL_SCORE_INCREASE

#define GRAPH_OPERATION_NONE     0
#define GRAPH_OPERATION_ADDITION 1
#define GRAPH_OPERATION_REMOVAL  2
#define GRAPH_OPERATION_REVERSAL 3

MetricDescription metricDescription[3] = {
  {"K2",&K2ComputeLogGains,NULL,NULL,&K2ComputeIsolatedNodeContribution},
  {"MDL",&mdlComputeAdditionGains,&mdlComputeRemovalGains,NULL,&mdlComputeIsolatedNodeContribution},
  {"BIC",&bicComputeAdditionGains,&bicComputeRemovalGains,NULL,&bicComputeIsolatedNodeContribution}
};

MetricDescription *metric;

int boaRecombination(Population *parents, Population *children, long M, RecombinationParams *params)
{
  long i,m;
  int k;
  long l;
  long totalSize;
  int  cluster;
  Population clusterParents;

  double *nodeContribution;

  long N;
  int  n,numDiscrete,numContinuous;
  int maxIncoming;
  int nMaxIncoming;

  int   **group, *groupSize;
  long   *numGroupInstances,*numGroupInstances0,maxNumGroupInstances;
  double **marginalAll, **marginalAllButOne;
  long   *count;

  char *added;
  int   numAdded;
  char  canAdd;
  char *x;
  double prob1;
  int   where;
  int   position;

  double **gainAddition;
  double **gainRemoval;
  double **gainReversal;

  char   **additionPossible;
  char   **removalPossible;
  char   **reversalPossible;

  double maxGain;
  int    maxType;
  int    maxFrom;
  int    maxTo;
  char  *full;
  double tmpGain;

  char finito;
  
  double totalScoreIncrease;

  long *normalizedSize;
  int numClusters;

  AcyclicOrientedGraph *G;

  Individual *child;

#ifdef BOA_PRINT_PARENTS

  // print the parent population

  printPopulation(stdout,parents);

#endif

  // initialize some variables (just to make it faster)

  N            = parents->N;
  n            = parents->n;
  numDiscrete  = parents->numDiscrete;
  numContinuous= parents->numContinuous;
  maxIncoming  = params->maxIncoming;
  nMaxIncoming = n*maxIncoming;
  numClusters  = params->numClusters;

///  printf("Check 0\n");

  // if maxIncoming is 0 then the umda recombination is the one we want

  //  if ((maxIncoming==0)||(!params->allowAdditions))
  //  umdaRecombination(parents,children,M,params);

  // allocate the children population

  allocatePopulation(children,M,parents->numDiscrete,parents->numContinuous);

  // allocate the array for network score contribution of each node

  nodeContribution = (double*) Calloc(n,sizeof(double));

///  printf("Check 0\n");

  // allocate the array for the gains with respect to the addition of each of the edges

  gainAddition     = (double**) Calloc(n,sizeof(double*));
  gainRemoval      = (double**) Calloc(n,sizeof(double*));
  gainReversal     = (double**) Calloc(n,sizeof(double*));
  additionPossible = (char**)   Calloc(n,sizeof(char*));
  removalPossible  = (char**)   Calloc(n,sizeof(char*));
  reversalPossible = (char**)   Calloc(n,sizeof(char*));

  for (k=0; k<n; k++)
    {
      gainAddition[k]     = (double*) Calloc(n,sizeof(double));
      gainRemoval[k]      = (double*) Calloc(n,sizeof(double));
      gainReversal[k]     = (double*) Calloc(n,sizeof(double));
      additionPossible[k] = (char*)   Malloc(n);
      removalPossible[k]  = (char*)   Malloc(n);
      reversalPossible[k] = (char*)   Malloc(n);
    };

  // allocate the memory for the matrix of groups, their sizes and marginal frequencies
  
  group              = (int**) Calloc(n,sizeof(int*));
  groupSize          = (int*) Calloc(n,sizeof(int));
  numGroupInstances  = (long*) Calloc(n,sizeof(long));
  numGroupInstances0 = (long*) Calloc(n,sizeof(long));
  marginalAll        = (double**) Calloc(n,sizeof(double*));
  marginalAllButOne  = (double**) Calloc(n,sizeof(double*));

  for (k=0; k<n; k++)
    group[k]=(int*) Calloc(n,sizeof(int));

  // allocate the memory for some auxilary array used for ordering vertices topologically

  added = (char*) Malloc(n);

  // allocate memory for the array containing information about filled nodes (no more edge can go in those)

  full = (char*) Malloc(n);

  // allocate memory for and initialize dependency graph

  G = new AcyclicOrientedGraph(n);

  // size of the clusters

  normalizedSize = (long*) Calloc(numClusters,sizeof(long));

  // for all coordinates do the same - create the network and generate the new individuals
  
  m=0;

  // cluster the population
    
  if (numClusters>0)
    kMeansClusterPopulation(parents,params->numClusters,params->numRestarts,params->phenotypicClustering);
  else
    for (i=0; i<parents->N; i++)
      parents->individual[i].cluster=0;
  
      // compute the normalized size of the clusters (to be used to generate the kids)
      
  if (params->fitnessProportionalClusterReproduction)
    {
      double *avgFitness;
      double total;
	  
      avgFitness = (double*) Calloc(numClusters,sizeof(double));
	  
      for (i=0; i<numClusters; i++)
	avgFitness[i]=0;
	  
      for (l=0; l<parents->N; l++)
	avgFitness[parents->individual[l].cluster] += parents->individual[l].f;
	  
      total=0;
      for (i=0; i<numClusters; i++)
	{
	  if (parents->clusterSize[i]>0)
	    avgFitness[i] /= parents->clusterSize[i];
	  else
	    avgFitness[i] = 0;
	      
	  if (avgFitness<0)
	    avgFitness=0;
	      
	  total += avgFitness[i];
	};
	  
      totalSize=0;
      for (i=0; i<numClusters; i++)
	{
	  normalizedSize[i] = long(M*avgFitness[i]/total);
	  totalSize += normalizedSize[i];
	}
	  
      normalizedSize[numClusters-1] += M-totalSize;
	  
      Free(avgFitness);
    }
  else
    {
      totalSize=0;
      for (i=0; i<numClusters; i++)
	{
	  normalizedSize[i] = long(double(parents->clusterSize[i])*M/parents->N);
	  totalSize += normalizedSize[i];
	}
      normalizedSize[numClusters-1] += M-totalSize;
    };

  for (cluster=0; cluster<numClusters; cluster++)
    if (parents->clusterSize[cluster]>0)
      {

	///	printf("Processing cluster %u of size %lu:\n",cluster,parents->clusterSize[cluster]);
///  printf("Check 1\n");

	// move the parents from this cluster to the clusterPopulation

	allocatePopulation(&clusterParents,parents->clusterSize[cluster],parents->numDiscrete,parents->numContinuous);
      
	l=0;
	for (i=0; i<parents->N; i++)
	  if (parents->individual[i].cluster==cluster)
	    copyIndividual(&(clusterParents.individual[l++]),&(parents->individual[i]),numDiscrete,numContinuous);
	N=parents->clusterSize[cluster];

	// allocate and compute univariate frequencies (we might need them in generating weird instances)

	double *p0;
	double *p1;

	allocateUMF(&p0,&p1,n);
	calculateUMF(&clusterParents,p0,p1);
	
///  printf("Check 2\n");

/* 	for (int ii=0; ii<n; ii++) */
/* 	  printf("%f ",p1[ii]); */
/* 	printf("\n"); */

	// 	  printf("   %lu members in the population found.\n",l);
	// 	  printf("   %lu members are to be generated.\n",normalizedSize[cluster]);

	// delete all edges
    
	G->removeAllEdges();

	// no edge is full so far

	for (k=0; k<n; k++)
	  full[k]=0;
	
	// all edge additions are possible, and all removals and reversals are not
	// (Note: if not starting with an empty network, we would have to change this)
	
	for (k=1; k<n; k++)
	  for (l=0; l<k; l++)
	    {
	      additionPossible[k][l] = additionPossible[l][k] = 1;
	      removalPossible[k][l]  = removalPossible[l][k]  = 0;
	      reversalPossible[k][l] = reversalPossible[l][k] = 0;
	    };

///  printf("Check 3\n");

	// compute the node contributions for all edges

	for (k=0; k<n; k++)
	  nodeContribution[k]=(*metric->boaComputeIsolatedNodeContribution)(k,&clusterParents);

	// compute the gains for all possible edge additions (no need to check most of the conditions, 
	// the first edge can be added anywhere (maxIncoming>=1, otherwise the umd recombination would 
	// be performed) (Note: If not starting with empty network, we should add here the reversal gains
	// and removal gains recompute too...)

	if (params->allowAdditions)	    
	  for (k=0; k<n; k++)
	    boaRecomputeAdditionGains(0,k,0,nodeContribution,gainAddition,gainRemoval,gainReversal,G,&clusterParents,params);
    
	// let's get greedy

	finito = 0;
	totalScoreIncrease = 0;

///  printf("Check 4\n");

	for (numAdded=0; (numAdded<nMaxIncoming)&&(!finito); numAdded++)
	  {
	    maxType = GRAPH_OPERATION_NONE;
	    maxGain = -1;
	    maxFrom = 0;
	    maxTo   = 0;

///            printf("Check A\n");

	    for (k=1; k<n; k++)
	      for (l=0; l<k; l++)
		{
		  // what about an edge from k to l?
	      
		  if ((params->allowAdditions)&&(additionPossible[k][l])&&(gainAddition[k][l]>maxGain))
		    {
		      maxType = GRAPH_OPERATION_ADDITION;
		      maxFrom = k;
		      maxTo   = l;
		      maxGain = gainAddition[k][l];
		    };

		  if ((params->allowRemovals)&&(removalPossible[k][l])&&(gainRemoval[k][l]>maxGain))
		    {
		      maxType = GRAPH_OPERATION_REMOVAL;
		      maxFrom = k;
		      maxTo   = l;
		      maxGain = gainRemoval[k][l];		  
		    };

		  ///	      printf("Checking reversal %u,%u with gain %f which status is %u\n",k,l,gainReversal[k][l],reversalPossible[k][l]);
		  if ((params->allowReversals)&&(reversalPossible[k][l])&&(gainReversal[k][l]>maxGain))
		    {
		      ///		  printf("updating reversal max...\n");
		      maxType = GRAPH_OPERATION_REVERSAL;
		      maxFrom = k;
		      maxTo   = l;
		      maxGain = gainReversal[k][l];		  
		    };
	      
		  // what about an edge from l to k?

		  if ((params->allowAdditions)&&(additionPossible[l][k])&&(gainAddition[l][k]>maxGain))
		    {
		      maxType = GRAPH_OPERATION_ADDITION;
		      maxFrom = l;
		      maxTo   = k;
		      maxGain = gainAddition[l][k];
		    };

		  if ((params->allowRemovals)&&(removalPossible[l][k])&&(gainRemoval[l][k]>maxGain))
		    {
		      maxType = GRAPH_OPERATION_REMOVAL;
		      maxFrom = l;
		      maxTo   = k;
		      maxGain = gainRemoval[l][k];		  
		    };

		  if ((params->allowReversals)&&(reversalPossible[l][k])&&(gainReversal[l][k]>maxGain))
		    {
		      ///		  printf("updating reversal max...\n");
		      maxType = GRAPH_OPERATION_REVERSAL;
		      maxFrom = l;
		      maxTo   = k;
		      maxGain = gainReversal[l][k];		  
		    };
		};

///            printf("Check B\n");

#ifdef BOA_PRINT_MAXGAIN

	    // print out a maximal gain for the best graph operation 

	    printf("\n\nMax gain -> %0.30lf (edge %u-%u) (and %0.30lf for the opposite addition)\n\n",maxGain,maxFrom,maxTo,gainAddition[maxTo][maxFrom]);
	    getchar();
#endif

	    // perform the operation that gives the biggest score improvement, if this is positive

	    if (maxGain>0)
	      {
		totalScoreIncrease += maxGain;
	    
		switch (maxType) 
		  {

		  case GRAPH_OPERATION_ADDITION:
		
                   //printf("Operator: Addition %u->%u (%f)\n",maxFrom,maxTo,maxGain);
		    // 		// !!!!!!!! don't swap really
		    // 		swapInt(&maxFrom,&maxTo);
		    // 		maxGain = -maxGain;
		    // 		// !!!!!!!!

		    // add the new edge
		
		    G->addEdge(maxFrom,maxTo);

		    // the contribution of node maxTo changes

		    nodeContribution[maxTo] += maxGain;

		    // recompute additions and, possibly, removals and reversals

		    boaRecomputeAdditionGains(maxFrom,maxTo,maxGain,nodeContribution,gainAddition,gainRemoval,gainReversal,G,&clusterParents,params);

		    if ((params->allowRemovals)||(params->allowReversals))
		      boaRecomputeRemovalGains(maxFrom,maxTo,maxGain,nodeContribution,gainAddition,gainRemoval,gainReversal,G,&clusterParents,params);

		    if (params->allowReversals)
		      boaRecomputeReversalGains(maxFrom,maxTo,maxGain,nodeContribution,gainAddition,gainRemoval,gainReversal,G,&clusterParents,params);
		
		    // update what is possible now

		    boaUpdateAfterAddition(maxFrom,maxTo,additionPossible,removalPossible,reversalPossible,G,full,params->maxIncoming);

		    break;
		
		  case GRAPH_OPERATION_REMOVAL:
	
///                    printf("Operator: Removal (%f)\n",maxGain);
	
		    ///		printf("Removing the edge %u->%u\n",maxFrom,maxTo);
		    ///		getchar();
		
		    // remove the chosen edge
		
		    G->removeEdge(maxFrom,maxTo);

		    // the contribution of node maxTo changes

		    nodeContribution[maxTo] += maxGain;

		    // update additions and, possibly, removals and reversals
		
		    boaRecomputeAdditionGains(maxFrom,maxTo,maxGain,nodeContribution,gainAddition,gainRemoval,gainReversal,G,&clusterParents,params);

		    if ((params->allowRemovals)||(params->allowReversals))
		      boaRecomputeRemovalGains(maxFrom,maxTo,maxGain,nodeContribution,gainAddition,gainRemoval,gainReversal,G,&clusterParents,params);

		    if (params->allowReversals)
		      boaRecomputeReversalGains(maxFrom,maxTo,maxGain,nodeContribution,gainAddition,gainRemoval,gainReversal,G,&clusterParents,params);
		
		    // update what is and what is not possible to do now

		    boaUpdateAfterRemoval(maxFrom,maxTo,additionPossible,removalPossible,reversalPossible,G,full,params->maxIncoming);

		    break;

		  case GRAPH_OPERATION_REVERSAL:

///                    printf("Operator: Reversal (%f)\n",maxGain);

		    ///		printf("Reversing the edge %u->%u\n",maxFrom,maxTo);
		    ///		getchar();

		    // reverse the chosen edge

		    G->reverseEdge(maxFrom,maxTo);

		    // need temporary gains

		    tmpGain = gainRemoval[maxFrom][maxTo];

		    // the contribution of nodes maxTo and maxFrom changes

		    nodeContribution[maxTo]   += tmpGain;
		    nodeContribution[maxFrom] += maxGain-tmpGain;

		    // update additions and so on, again, now we need to do more than is usual (in both affected nodes)

		    boaRecomputeAdditionGains(maxFrom,maxTo,tmpGain,nodeContribution,gainAddition,gainRemoval,gainReversal,G,&clusterParents,params);
		    boaRecomputeAdditionGains(maxTo,maxFrom,maxGain-tmpGain,nodeContribution,gainAddition,gainRemoval,gainReversal,G,&clusterParents,params);

		    if ((params->allowRemovals)||(params->allowReversals))
		      {
			boaRecomputeRemovalGains(maxFrom,maxTo,tmpGain,nodeContribution,gainAddition,gainRemoval,gainReversal,G,&clusterParents,params);
			boaRecomputeRemovalGains(maxTo,maxFrom,maxGain-tmpGain,nodeContribution,gainAddition,gainRemoval,gainReversal,G,&clusterParents,params);
		      };
		
		    if (params->allowReversals)
		      {
			boaRecomputeReversalGains(maxFrom,maxTo,maxGain,nodeContribution,gainAddition,gainRemoval,gainReversal,G,&clusterParents,params);
			boaRecomputeReversalGains(maxTo,maxFrom,maxGain-tmpGain,nodeContribution,gainAddition,gainRemoval,gainReversal,G,&clusterParents,params);
		      };

		    // update what is and what is not possible to do now

		    boaUpdateAfterReversal(maxFrom,maxTo,additionPossible,removalPossible,reversalPossible,G,full,params->maxIncoming);
		  }
	      }
	    else
	      finito = 1;
	  };

///            printf("Check C\n");

#ifdef BOA_PRINT_TOTAL_SCORE_INCREASE

	printf("Total score increase: %f\n",totalScoreIncrease);

#endif

#ifdef BOA_PRINT_COINCIDENCE
     
	// print the resulting coincidence matrix

	G->printCoincidenceMatrix(stdout);

#endif

	// create the array of incoming edges numbers for all vartices (with it itself
	// on the first position) and the number of them. 
	// + allocate the memory for marginal frequencies

	maxNumGroupInstances = 0;

	for (k=0; k<n; k++)
	  {
	    // create k-th group (corresponding to the k-th variable)

	    group[k][0]=k;
	    if (G->getNumIn(k)>0)
	      memcpy(&(group[k][1]),G->getParentList(k),G->getNumIn(k)*sizeof(int));
	    groupSize[k]=G->getNumIn(k)+1;

	    // compute the number of its instances

	    numGroupInstances[k] = 1<<groupSize[k];
	    numGroupInstances0[k] = numGroupInstances[k]>>1;

	    // allocate the memory for marginals and counts

	    marginalAll[k] = (double*) Calloc(numGroupInstances[k],sizeof(double));
	    marginalAllButOne[k] = (double*) Calloc(numGroupInstances0[k],sizeof(double));

	    // update maxNumGroupInstances

	    if (numGroupInstances[k]>maxNumGroupInstances)
	      maxNumGroupInstances=numGroupInstances[k];
	  }

	count = (long*) Calloc(maxNumGroupInstances,sizeof(long));

	// order groups of positions and the corresponding arrays topologically
	// (needed for generation of new instances)

	for (k=0; k<n; k++)
	  added[k]=0;
	numAdded=0;

	while (numAdded<n)
	  {
	    for (k=numAdded; k<n; k++)
	      {
		canAdd=1;
	    
		for (l=1; l<groupSize[k]; l++)
		  if (!added[group[k][l]])
		    canAdd=0;
	    
		if (canAdd)
		  {
		    added[group[k][0]]=1;
		    if (k!=numAdded)
		      {
			swapPointers((void**) &(group[k]),(void**) &(group[numAdded]));
			swapInt(&(groupSize[k]), &(groupSize[numAdded]));
			swapLong(&(numGroupInstances[k]),&(numGroupInstances[numAdded]));
			swapLong(&(numGroupInstances0[k]),&(numGroupInstances0[numAdded]));
			swapPointers((void**) &(marginalAll[k]), (void**) &(marginalAll[numAdded]));
			swapPointers((void**) &(marginalAllButOne[k]), (void**) &(marginalAllButOne[numAdded]));
		      }
		
		    numAdded++;
		  }
	      }
	  }

	// display the covered dependencies (if wanted)

	if (params->displayDependencies)
	  G->printCoincidenceMatrix(stdout);
    
#ifdef BOA_PRINT_GROUPS

	// output the info about groups 

	printf("Groups:\n");
	for (k=0; k<n; k++)
	  {
	    printf("( ");
	    for (l=0; l<groupSize[k]; l++)
	      printf("%u ",group[k][l]);
	    printf(")\n");
	  }
	printf("\n");
	getchar();
    
#endif

	// calculate the marginal frequencies for created groups of positions
   
	for (k=0; k<n; k++)
	  {
	    computeCounts(group[k],groupSize[k],&clusterParents,count);
	    
	    for (l=0; l<numGroupInstances[k]; l++)
	      {
		///	    printf("Count[%lu]=%lu\n",l,count[l]);
		marginalAll[k][l] = (double) count[l]/N;
		///	    getchar();
	      };

	    for (l=0; l<numGroupInstances0[k]; l++)
	      marginalAllButOne[k][l] = marginalAll[k][l]+marginalAll[k][l+numGroupInstances0[k]];
	  };

   
	//     for (k=0; k<n; k++)
	//       {
	// 	printf("Group (");
	// 	for (l=0; l<groupSize[k]; l++)
	// 	  printf("%u ",group[k][l]);
	// 	printf(") -> ");
	// 	for (l=0; l<1<<groupSize[k]; l++)
	// 	  printf("%1.5f ",marginalAll[k][l]);
	// 	printf(" # ");
	// 	for (l=0; l<1<<(groupSize[k]-1); l++)
	// 	  printf("%1.5f ",marginalAllButOne[k][l]);
	// 	printf("\n");
	//       };
	//     getchar();
	

	// generate the ith coordinate (chromosome) for all dhildren

	for (i=0; (i<normalizedSize[cluster])&&(m<M); i++)
	  {
	    child = &(children->individual[m]);
	    x = child->chromosome;

	    for (k=0; k<n; k++)
	      {
		position = group[k][0];

		if (groupSize[k]==1)
		  prob1 = marginalAll[k][1];
		else
		  {
		    x[position] = 0;
		    where       = indexedBinaryToInt(x,group[k],groupSize[k]);
		    if (marginalAllButOne[k][where]>0)
		      prob1       = 1-marginalAll[k][where]/marginalAllButOne[k][where];
		    else
		      {
/* 			printf("Here we go.\n"); */
/* 			printf("%u <- ",group[k][0]); */
/* 			for (int ii=1; ii<groupSize[k]; ii++) */
/* 			  printf("%u ",group[k][ii]); */
/* 			printf("\n"); */
			
/* 			int mm=0; */
/* 			for (int ii=0; ii<N; ii++) */
/* 			  mm+=parents->individual[ii].chromosome[position]; */
/* 			printf("m(position=1) = %u\n",mm); */
/* 			getchar(); */
			prob1       = p1[position]; // in this case the parents never appear and we use univariate frequency;
			///			printf("Filling in %f as opposed to %f\n",p1[position],0.5);
		      };
		  }

		if (drand()<prob1)
		  x[position]=1;
		else
		  x[position]=0;
	      }

	    m++;
	  }
	
	// free the memory used by marginal frequencies (for groups of bits) and the "count" array
	
	for (k=0; k<n; k++)
	  {
	    Free(marginalAll[k]);
	    Free(marginalAllButOne[k]);
	  };
	
	Free(count);
	
	freePopulation(&clusterParents);

	freeUMF(&p0,&p1,n);
	
      };

  // the fitness has to be evaluated for all children,

  for (i=0; i<M; i++)
    {
      child = &(children->individual[i]);
      child->fCalculated = 0;
    }    

  // free the memory used by the gain array

  for (k=0; k<n; k++)
    {
      Free(gainAddition[k]);
      Free(gainRemoval[k]);
      Free(gainReversal[k]);
      Free(additionPossible[k]);
      Free(removalPossible[k]);
      Free(reversalPossible[k]);
    };
  Free(gainAddition);
  Free(gainRemoval);
  Free(gainReversal);
  Free(additionPossible);
  Free(removalPossible);
  Free(reversalPossible);

  // free the memory used by the matrix of groups, their sizes, and marginal frequencies

  for (k=0; k<n; k++)
    Free(group[k]);
  Free(group);
  Free(groupSize);
  Free(marginalAll);
  Free(marginalAllButOne);

  // free memory occupied by dependency graph

  delete G;

  // free memory used by "numGroupInstances" and "numGroupInstances0"

  Free(numGroupInstances);
  Free(numGroupInstances0);

  // free memory used by an auxilary array "added"

  Free(added);

  // free memory used by an auxilary array "full"

  Free(full);

  // offspring has not been evaluated yet

  children->evaluated = 0;

  // free the network score contribution of each node array

  Free(nodeContribution);

  // free the array with sizes of the clusters

  Free(normalizedSize);

  // get back

  return 0;
};

//=============================================================================================

//==========================================================

int boaUpdateAfterAddition(int newFrom, 
			   int newTo, 
			   char **additionPossible, 
			   char **removalPossible, 
			   char **reversalPossible, 
			   AcyclicOrientedGraph *G, 
			   char *full,
			   int maxIncoming)
{
  int k,l,n;

  n = G->size();

  // the removal is possible and the reversal can be also

  removalPossible[newFrom][newTo] = 1;

  // carefully with this.....would get into trouble if updated addion possibilities first....
  // could be also ((G->canReverseEdge[newFrom][newTo])&&(!full[newTo])) - safer and slower

  ///  printf("New reversal possibility %u,%u = %u\n",newFrom,newTo,additionPossible[newTo][newFrom]);
  reversalPossible[newFrom][newTo] = additionPossible[newTo][newFrom];
  
  // can't add anymore

  additionPossible[newFrom][newTo] = 0;
  additionPossible[newTo][newFrom] = 0;

  // is the ending node of the new edge full yet?
  
  full[newTo]=(G->getNumIn(newTo)==maxIncoming);

  // if full no more edges can end here

  if (full[newTo])
    for (k=0; k<n; k++)
      {
	additionPossible[k][newTo]=0;
	reversalPossible[newTo][k]=0;
      };
 
  // check for cycles with the new edge...

  for (k=1; k<n; k++)
    for (l=0; l<k; l++)
      if (((k!=newFrom)||(l!=newTo))&&((k!=newTo)||(l!=newFrom)))
	{
	  // does the new edge forbid creating an edge k,l by means of a path that might create a cycle with this?
	
	  if ((additionPossible[k][l]>0)&&(G->existsPath(l,newFrom)&&(G->existsPath(newTo,k))))
	    {
	      additionPossible[k][l]=0;
	      reversalPossible[l][k]=0;
	    };
	    
	  // does the new edge forbid creating an edge l,k by means of a path that might create a cycle with this?
	
	  if ((additionPossible[l][k]>0)&&(G->existsPath(k,newFrom)&&(G->existsPath(newTo,l))))
	    {
	      additionPossible[l][k]=0;
	      reversalPossible[k][l]=0;
	    };
	}
     
  return 0;
};

//==========================================================

int boaUpdateAfterRemoval(int newFrom, 
			  int newTo, 
			  char **additionPossible, 
			  char **removalPossible, 
			  char **reversalPossible, 
			  AcyclicOrientedGraph *G, 
			  char *full,
			  int maxIncoming)
{
  int k,l,n;

  // initialize the variables

  n = G->size();

  // can't remove or reverse the thing anymore

  removalPossible[newFrom][newTo]  = 0;
  reversalPossible[newFrom][newTo] = 0;

  // can't be full after removing...
  
  full[newTo]=0;

  // allows the edges that got possible to be added after removal to be added

  for (k=1; k<n; k++)
    for (l=0; l<k; l++)
      {
	if ((!additionPossible[k][l])&&(!full[l])&&(!G->connected(k,l))&&(!G->existsPath(l,k)))
	  additionPossible[k][l]=1;
	
	if ((!additionPossible[l][k])&&(!full[k])&&(!G->connected(l,k))&&(!G->existsPath(k,l)))
	  additionPossible[l][k]=1;

      };

  return 0;
};

//==========================================================

int boaUpdateAfterReversal(int newFrom, 
			   int newTo, 
			   char **additionPossible, 
			   char **removalPossible, 
			   char **reversalPossible, 
			   AcyclicOrientedGraph *G, 
			   char *full,
			   int maxIncoming)
{
  //what has happened after we removed an edge newFrom->newTo

  boaUpdateAfterRemoval(newFrom,newTo,additionPossible,removalPossible,reversalPossible,G,full,maxIncoming);

  // what has happened after we added edge newTo->newFrom

  boaUpdateAfterAddition(newTo,newFrom,additionPossible,removalPossible,reversalPossible,G,full,maxIncoming);

  // get back

  return 0;
};

//==========================================================

int boaRecomputeAdditionGains(int from,
			      int to,
			      double gain,
			      double *nodeContribution,
			      double **gainAddition, 
			      double **gainRemoval, 
			      double **gainReversal, 
			      AcyclicOrientedGraph *G, 
			      Population *population,
			      RecombinationParams *params)
{
  int parentCount;
  int *parentList;
  int numUpdated;
  int *updateIdx;
  int n,k;

  // initialize the variables

  n = population->n;

  // initialize the variables
      
  parentCount = G->getNumIn(to);
  
  // initialize the parent list
  
  parentList = G->getParentList(to);
  
  // allocate memory for the index of nodes from which the edge is to be updated
  
  updateIdx = (int*) Calloc(n-parentCount,sizeof(int));
      
  // create the index of those nodes
      
  numUpdated = 0;
      
  for (k=0; k<n; k++)
    if ((k!=to)&&(!G->connected(k,to)))
      updateIdx[numUpdated++] = k;
      
  // compute the gains for edges from all updateIdx to i
  
  if (numUpdated>0)
    (*metric->boaComputeAdditionGains)(to,nodeContribution[to],gainAddition,updateIdx,numUpdated,parentList,parentCount,population,params);
  
  // if the prior network is used, penalize the additions when unmatched with the prior network
  
  for (k=0; k<numUpdated; k++)
    gainAddition[updateIdx[k]][to] -= boaGetPenalty(updateIdx[k],to);

  // free memory used by list of beginning nodes of the edges to be updated
  
  Free(updateIdx);

  // get back
  
  return 0;
};

//==========================================================

int boaRecomputeRemovalGains(int from,
			     int to,
			     double gain,
			     double *nodeContribution,
			     double **gainAddition, 
			     double **gainRemoval, 
			     double **gainReversal, 
			     AcyclicOrientedGraph *G, 
			     Population *population,
			     RecombinationParams *params)
{
  int parentCount;
  int *parentList;

  // initialize the variables
      
  parentCount = G->getNumIn(to);
      
  // initialize the parent list

  parentList = G->getParentList(to);
  
  // compute the gains for edges from all updateIdx to i

  if (parentCount>0)
    (*metric->boaComputeRemovalGains)(to,nodeContribution[to],gainRemoval,parentList,parentCount,population,params);
  
  // if the prior network is used, penalize the removals when unmatched with the prior network
  
///  for (k=0; k<parentCount; k++)
///    gainRemoval[parentList[k]][to] -= boaGetPenalty(parentList[k],to);
  
  // get back

  return 0;
};

//==========================================================

int boaRecomputeReversalGains(int from,
			      int to,
			      double gain,
			      double *nodeContribution,
			      double **gainAddition, 
			      double **gainRemoval, 
			      double **gainReversal, 
			      AcyclicOrientedGraph *G, 
			      Population *population,
			      RecombinationParams *params)
{
  int k,n;

  // initialize the variables

  n = G->size();

  // compute reversal gains as a function of removal and addition gains (that are all ready waiting)
  // this is a very delicate thing....intricate....not sure it works

  for (k=0; k<n; k++)
    if (G->connected(k,to))
      {
	gainReversal[k][to] -= gain;
	///	printf("New reversal gain %u,%u = %lf\n",k,to,gainReversal[k][to]);
      }

  for (k=0; k<n; k++)
    if (G->connected(to,k))
      {
	gainReversal[to][k] = gainRemoval[to][k]+gainAddition[k][to];
	///	printf("New reversal gain %u,%u = %lf\n",to,k,gainReversal[to][k]);
      }

  // get back

  return 0;
};

// ========================================================

double boaGetPenalty(int from, int to)
{
  // add penalty function and stuff !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  return 0;
};

//==========================================================

int boaSetMetric(int n)
{
  metric = &(metricDescription[n]);

  return 0;
};

//==========================================================

char *boaGetMetricDescription(int n)
{
  return metric->description;
};
