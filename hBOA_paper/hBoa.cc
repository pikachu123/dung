#include <stdio.h>
#include <string.h>

#include "hBoa.h"
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
#include "frequencyTree.h"
#include "index.h"

//#define HBOA_PRINT_PEARSON
//#define HBOA_PRINT_PARENTS
//#define HBOA_DISPLAY_MODEL
//#define HBOA_DISPLAY_TOPOLOGICAL_ORDERING
//#define HBOA_PRINT_GROUPS
//#define HBOA_PRINT_MAXGAIN
//#define HBOA_PRINT_TOTAL_SCORE_INCREASE

#define getPenalty(a,b) 0;                      // have to finish!!!!!!!!!!!

#define GRAPH_OPERATION_NONE     0
#define GRAPH_OPERATION_ADDITION 1
#define GRAPH_OPERATION_REMOVAL  2
#define GRAPH_OPERATION_REVERSAL 3
#define GRAPH_OPERATION_JOINT    4

#define assignMaxAddition(k_,l_) { maxType = GRAPH_OPERATION_ADDITION; maxFrom = k_; maxTo   = l_; maxGain = operatorGain.addition[k_][l_]; }
#define assignMaxRemoval(k_,l_)  { maxType = GRAPH_OPERATION_REMOVAL; maxFrom = k_; maxTo   = l_; maxGain = operatorGain.removal[k_][l_]; }
#define assignMaxReversal(k_,l_) { maxType = GRAPH_OPERATION_REVERSAL; maxFrom = k_; maxTo   = l_; maxGain = operatorGain.reversal[k_][l_]; }
#define assignMaxJoint(k_,l_)    { maxType = GRAPH_OPERATION_JOINT; maxFrom = k_; maxTo   = l_; maxGain = operatorGain.joint[k_][l_]; }

int hboaRecombination(Population *parents, Population *children, long M, RecombinationParams *params)
{
  long         i;
  int          k;
  long         l;
  int          m;

  GroupInformation groupInformation;

  int    *fixedIndex;
  char   *fixedValue;
  int     numFixed;
  int     numVars;

  int  n;
  int maxIncoming;
  int nMaxIncoming;

  double *p0,*p1;

  char  canAdd;
  char  *added;
  int   numAdded;
  int   *index;
  char   *x;

  OperatorGain          operatorGain;
  OperatorApplicability operatorApplicability;

  double maxGain;
  int    maxType;
  int    maxFrom;
  int    maxTo;

  char finito;

  double totalScoreIncrease;

  AcyclicOrientedGraph *G;
  AcyclicOrientedGraph *newG;

  Individual *child;

  int shiftedMaxFrom;

  printf("Recombination entered.\n");

#ifdef HBOA_PRINT_PARENTS

  // print the parent population

  printPopulation(stdout,parents);

#endif

  // initialize some variables (just to make it faster)

  n                             = parents->n;
  maxIncoming                   = params->maxIncoming;
  nMaxIncoming                  = n*maxIncoming;
  groupInformation.maxGroupSize = params->maxGroupSize;

  // if maxIncoming is 0 and the maximal size of the groups is 0
  // then the umda recombination is the one we want. Also if no
  // additions are allowed...
  /*
  if (((maxIncoming==0) && (groupInformation.maxGroupSize==0))||(!params->allowAdditions))
    umdaRecombination(parents,children,M,params);
  */
  // ----------------------------- allocate memory --------------------------------

  // allocate the children population

  allocatePopulation(children,M,parents->numDiscrete,parents->numContinuous);

  // and some more for information about each group of variables in the network

  groupInformation.groupSize         = (int*)    Calloc(n,sizeof(int));
  groupInformation.groupIndex        = (int**)   Calloc(n,sizeof(int*));
  groupInformation.frequencyTree     = (FrequencyTree**) Calloc(n,sizeof(FrequencyTree*));
  groupInformation.numInstances      = (int*)     Calloc(n,sizeof(int));
  groupInformation.instances         = (char***)  Calloc(n,sizeof(char**));
  groupInformation.frequencies       = (double**) Calloc(n,sizeof(double*));
  groupInformation.full              = (char*)    Malloc(n);
  groupInformation.indexSize         = (int*)     Calloc(n,sizeof(int));
  groupInformation.index             = (int**)    Calloc(n,sizeof(int*));
  groupInformation.parentIndexSize   = (int*)     Calloc(n,sizeof(int));
  groupInformation.contribution      = (double*)  Calloc(n,sizeof(double));
    
  for (k=0; k<n; k++)
    groupInformation.groupIndex[k]   = (int*)  Calloc(n,sizeof(int));
    
  // allocate memory for the array of fixed positions and their values

  fixedIndex = (int*)  Calloc(n,sizeof(int));
  fixedValue = (char*) Malloc(n);

  // allocate memory for the array of gains with respect to the addition,
  // removal, and reversal of each edge

  operatorGain.addition     = (double**) Calloc(n,sizeof(double*));
  operatorGain.removal      = (double**) Calloc(n,sizeof(double*));
  operatorGain.reversal     = (double**) Calloc(n,sizeof(double*));
  operatorGain.joint        = (double**) Calloc(n,sizeof(double*));

  for (k=0; k<n; k++)
    {
      operatorGain.addition[k]     = (double*) Calloc(n,sizeof(double));
      operatorGain.removal[k]      = (double*) Calloc(n,sizeof(double));
      operatorGain.reversal[k]     = (double*) Calloc(n,sizeof(double));
      operatorGain.joint[k]        = (double*) Calloc(n,sizeof(double));
    };

  // allocate memory for the array of operation-applicability

  operatorApplicability.addition = (char**) Calloc(n,sizeof(char*));
  operatorApplicability.removal  = (char**) Calloc(n,sizeof(char*));
  operatorApplicability.reversal = (char**) Calloc(n,sizeof(char*));
  operatorApplicability.joint    = (char**) Calloc(n,sizeof(char*));

  for (k=0; k<n; k++)
    {
      operatorApplicability.addition[k] = (char*) Malloc(n);
      operatorApplicability.removal[k]  = (char*) Malloc(n);
      operatorApplicability.reversal[k] = (char*) Malloc(n);
      operatorApplicability.joint[k]    = (char*) Malloc(n);
    };

  // allocate the memory for some auxilary array used for ordering the groups topologically

  added = (char*) Malloc(n);
  index = (int*)  Calloc(n,sizeof(int));

  // ---------------------------- initialize groups ------------------------------

  // all groups contain just one variable with both possible values...0 and 1
  // (if someone has converged forget about him totally and put him in the
  // array of fixed variables)

  allocateUMF(&p0,&p1,n);
  calculateUMF(parents,p0,p1);

  // --------------------- initialize the graph and related ------------------------
    
  groupInformation.numGroups = numFixed = 0;
	    
  for (k=0; k<n; k++)
    {
      // if position not converged add it to the groups (with its values 0,1)
		    
      if ((p0[k]!=0)&&(p1[k]!=0))
	{
	  groupInformation.groupSize[groupInformation.numGroups]         = 1;
	  groupInformation.groupIndex[groupInformation.numGroups][0]     = k;
	  groupInformation.numGroups++;
	}
      else
	{
	  fixedIndex[numFixed] = k;
	  fixedValue[numFixed] = (p0[k]>0)? 0:1;
			    
	  numFixed++;
	};
    };

  if (numFixed>0)
    printf("Number of fixed bits: %u\n",numFixed); 
	    
      // the number of variables we are taking into account

  numVars = n-numFixed;

  if (numVars>0)
    {
	  
      // allocate memory for and initialize dependency graph
	    
      G = new AcyclicOrientedGraph(groupInformation.numGroups);

      // delete all edges
	    
      G->removeAllEdges();

      // no edge is full so far if maxIncoming>0
  
      for (k=0; k<groupInformation.numGroups; k++)
	groupInformation.full[k]=(maxIncoming==0);

      // all edge additions are possible if maxIncoming>0, and all removals and
      // reversals are not
      // (Note: if not starting with an empty network, we would have to change this)

      for (k=1; k<groupInformation.numGroups; k++)
	{
	  for (l=0; l<k; l++)
	    {
	      operatorApplicability.addition[k][l] = operatorApplicability.addition[l][k] = (maxIncoming>0);
	      operatorApplicability.removal[k][l]  = operatorApplicability.removal[l][k]  = 0;
	      operatorApplicability.reversal[k][l] = operatorApplicability.reversal[l][k] = 0;
	      operatorApplicability.joint[k][l]    = operatorApplicability.joint[l][k]    = 0;
	    };
	};

      for (k=0; k<groupInformation.numGroups; k++)
	groupInformation.contribution[k] = mdlComputeGroupContribution(k,
								       numVars,
								       &groupInformation,
								       G,
								       parents,
								       params);

      // compute the gains for all possible edge additions (no need to check most of the conditions,
      // the first edge can be added anywhere (maxIncoming>=1, otherwise the umd recombination would
      // be performed) (Note: If not starting with an empty network, we should add here the reversal gains
      // and removal gains recompute too...)

      for (k=0; k<groupInformation.numGroups; k++)
	{
	  if (params->allowAdditions)
	    hboaRecomputeAdditionGains( 0,
					k,
					numVars,
					&groupInformation,
					&operatorGain,
					&operatorApplicability,
					G,
					parents,
					params);
 	
	  if (params->allowJoints)
	    hboaRecomputeJointGains( 0,
				     k,
				     numVars,
				     &groupInformation,
				     &operatorGain,
				     &operatorApplicability,
				     G,
				     parents,
				     params);
	};

      // ---------------------------- the greedy search -----------------------------
      // let's get greedy

      finito = 0;
      totalScoreIncrease = 0;

      for (numAdded=0; (numAdded<nMaxIncoming)&&(!finito); numAdded++)
	{
	  maxType = GRAPH_OPERATION_NONE;
	  maxGain = -1;
	  maxFrom = 0;
	  maxTo   = 0;
 	
	  for (k=1; k<groupInformation.numGroups; k++)
	    for (l=0; l<k; l++)
	      {
		// what about an edge from k to l?
 	
		if ((params->allowAdditions)&&
		    (operatorApplicability.addition[k][l])&&
		    (operatorGain.addition[k][l]>maxGain))
		  assignMaxAddition(k,l);

		if ((params->allowRemovals)&&
		    (operatorApplicability.removal[k][l])&&
		    (operatorGain.removal[k][l]>maxGain))
		  assignMaxRemoval(k,l);	

		if ((params->allowReversals)&&
		    (operatorApplicability.reversal[k][l])&&
		    (operatorGain.reversal[k][l]>maxGain))
		  assignMaxReversal(k,l);

		if ((params->allowJoints)&&
		    (operatorApplicability.joint[k][l])&&
		    (operatorGain.joint[k][l]>maxGain))
		  assignMaxJoint(k,l);
 	
		// what about an edge from l to k?

		if ((params->allowAdditions)&&
		    (operatorApplicability.addition[l][k])&&
		    (operatorGain.addition[l][k]>maxGain))
		  assignMaxAddition(l,k);

		if ((params->allowRemovals)&&
		    (operatorApplicability.removal[l][k])&&
		    (operatorGain.removal[l][k]>maxGain))
		  assignMaxRemoval(l,k);

		if ((params->allowReversals)&&
		    (operatorApplicability.reversal[l][k])&&
		    (operatorGain.reversal[l][k]>maxGain))
		  assignMaxReversal(l,k);
	      };

#ifdef HBOA_PRINT_MAXGAIN

	  // print out a maximal gain for the best graph operation

	  printf("\n\nMax gain -> %0.30lf (edge %u-%u) (and %0.30lf for the opposite addition)\n\n",maxGain,maxFrom,maxTo,gainAddition[maxTo][maxFrom]);
	  ///getchar();
#endif

	  //	  printf("Max operation is (%u,%u) of type %u with gain %f\n",maxFrom,maxTo,maxType,maxGain);
	  //	  getchar();

	  // perform the operation that gives the biggest score improvement, if this is positive

	  if (maxGain>0)
	    {
	      totalScoreIncrease += maxGain;
 	
	      switch (maxType)
		{

		  // -----------------------------------------------------------------------
				    
		case GRAPH_OPERATION_ADDITION:                                                       // checked completely

		  //		  printf("\n-----------------------------------------------------------------------\n\nDoing addition (%u,%u) with gain %f\n",maxFrom,maxTo,maxGain);
		  //		      printf("Adding the edge %u->%u with gain %f\n",maxFrom,maxTo,maxGain);
				    
		  // add the new edge
 		
		  G->addEdge(maxFrom,maxTo);

		  groupInformation.contribution[maxTo] += maxGain;

		  // update what is possible now

		  hboaUpdateAfterAddition( maxFrom,
					   maxTo,
					   &groupInformation,
					   &operatorApplicability,
					   G,
					   params);

		  // recompute additions and, possibly, removals and reversals

		  if ((params->allowAdditions)||(params->allowReversals))
		    hboaRecomputeAdditionGains( maxFrom,
						maxTo,
						numVars,
						&groupInformation,
						&operatorGain,
						&operatorApplicability,
						G,
						parents,
						params);

		  if ((params->allowRemovals)||(params->allowReversals))
		    hboaRecomputeRemovalGains( maxFrom,
					       maxTo,
					       numVars,
					       &groupInformation,
					       &operatorGain,
					       &operatorApplicability,
					       G,
					       parents,
					       params);

		  if (params->allowReversals)
		    hboaRecomputeReversalGains( maxFrom,
						maxTo,
						numVars,
						&groupInformation,
						&operatorGain,
						&operatorApplicability,
						G,
						parents,
						params);
 		
		  if (params->allowJoints)
		    hboaRecomputeJointGains( maxFrom,
					     maxTo,
					     numVars,
					     &groupInformation,
					     &operatorGain,
					     &operatorApplicability,
					     G,
					     parents,
					     params);

		  break;

		  // -----------------------------------------------------------------------
 		
		case GRAPH_OPERATION_REMOVAL:                                                   // completely checked
 		
		  printf("Removing the edge %u->%u with gain %f\n",maxFrom,maxTo,maxGain);
 		
		  // remove the chosen edge
 		
		  G->removeEdge(maxFrom,maxTo);

		  groupInformation.contribution[maxTo] += maxGain;

		  // update what is and what is not possible to do now

		  hboaUpdateAfterRemoval(maxFrom,
					 maxTo,
					 &groupInformation,
					 &operatorApplicability,
					 G,
					 params);

		  // update additions and, possibly, removals and reversals
 		
		  hboaRecomputeAdditionGains(maxFrom,
					     maxTo,
					     numVars,
					     &groupInformation,
					     &operatorGain,
					     &operatorApplicability,
					     G,parents,
					     params);

		  if ((params->allowRemovals)||(params->allowReversals))
		    hboaRecomputeRemovalGains(maxFrom,
					      maxTo,
					      numVars,
					      &groupInformation,
					      &operatorGain,
					      &operatorApplicability,
					      G,
					      parents,
					      params);

		  if (params->allowReversals)
		    hboaRecomputeReversalGains(maxFrom,
					       maxTo,
					       numVars,
					       &groupInformation,
					       &operatorGain,
					       &operatorApplicability,
					       G,
					       parents,
					       params);

		  if (params->allowJoints)
		    hboaRecomputeJointGains(maxFrom,
					    maxTo,
					    numVars,
					    &groupInformation,
					    &operatorGain,
					    &operatorApplicability,
					    G,
					    parents,
					    params);

		  break;

		  // -----------------------------------------------------------------------

		case GRAPH_OPERATION_REVERSAL:                                                  // checked completely

		  printf("Reversing the edge %u->%u with gain %f\n",maxFrom,maxTo,maxGain);
		  printf("The removal was %f and the addition %f\n",operatorGain.removal[maxFrom][maxTo],operatorGain.addition[maxTo][maxFrom]);

		  // reverse the chosen edge

		  G->reverseEdge(maxFrom,maxTo);

		  groupInformation.contribution[maxFrom] = mdlComputeGroupContribution(maxFrom,
										       numVars,
										       &groupInformation,
										       G,
										       parents,
										       params);

		  groupInformation.contribution[maxTo] = mdlComputeGroupContribution(maxTo,
										     numVars,
										     &groupInformation,
										     G,
										     parents,
										     params);
		      
		  // update what is and what is not possible to do now

		  hboaUpdateAfterReversal(maxFrom,
					  maxTo,
					  &groupInformation,
					  &operatorApplicability,
					  G,
					  params);

		  // update additions and so on, again, now we need to do more than is usual (in both affected nodes)
		  // (as if removed and added then instead of reversing, what is actually the same)

		  hboaRecomputeAdditionGains(maxTo,
					     maxFrom,
					     numVars,
					     &groupInformation,
					     &operatorGain,
					     &operatorApplicability,
					     G,
					     parents,
					     params);

		  hboaRecomputeAdditionGains(maxFrom,
					     maxTo,
					     numVars,
					     &groupInformation,
					     &operatorGain,
					     &operatorApplicability,
					     G,
					     parents,
					     params);

		  if ((params->allowRemovals)||(params->allowReversals))
		    {
		      hboaRecomputeRemovalGains(maxTo,
						maxFrom,
						numVars,
						&groupInformation,
						&operatorGain,
						&operatorApplicability,
						G,
						parents,
						params);

		      hboaRecomputeRemovalGains(maxFrom,
						maxTo,
						numVars,
						&groupInformation,
						&operatorGain,
						&operatorApplicability,
						G,
						parents,
						params);
		    };
 		
		  if (params->allowReversals)
		    {
		      hboaRecomputeReversalGains(maxTo,
						 maxFrom,
						 numVars,
						 &groupInformation,
						 &operatorGain,
						 &operatorApplicability,
						 G,
						 parents,
						 params);

		      hboaRecomputeReversalGains(maxFrom,
						 maxTo,
						 numVars,
						 &groupInformation,
						 &operatorGain,
						 &operatorApplicability,
						 G,
						 parents,
						 params);
		    };

		  if (params->allowJoints)
		    {
		      hboaRecomputeJointGains(maxFrom,
					      maxTo,
					      numVars,
					      &groupInformation,
					      &operatorGain,
					      &operatorApplicability,
					      G,
					      parents,
					      params);
		    };

		  break;

		  // -----------------------------------------------------------------------

		case GRAPH_OPERATION_JOINT:

		  ///		printf("Merging groups %u--%u\n",maxFrom,maxTo);
		  ///		getchar();

		  newG = new AcyclicOrientedGraph(G->size()-1);

		  // join the chosen groups

		  joinGroups(maxFrom,
			     maxTo,
			     &groupInformation,
			     &operatorApplicability,
			     G,
			     newG,
			     parents);

		  shiftedMaxFrom = (maxFrom>maxTo)? (maxFrom-1):maxFrom;
 	
		  // update what is and what is not possible to do now

		  hboaUpdateAfterJoint(shiftedMaxFrom,
				       0,
				       &groupInformation,
				       &operatorApplicability,
				       G,
				       params);

		  // recompute all. But can be optimized, not worth probably though. we'll see what future brings...

		  for (k=0; k<newG->size(); k++)
		    {
		      if (params->allowAdditions)
			hboaRecomputeAdditionGains(0,
						   k,
						   numVars,
						   &groupInformation,
						   &operatorGain,
						   &operatorApplicability,
						   G,parents,
						   params);
			
		      if ((params->allowRemovals)||(params->allowReversals))
			hboaRecomputeRemovalGains(0,
						  k,
						  numVars,
						  &groupInformation,
						  &operatorGain,
						  &operatorApplicability,
						  G,
						  parents,
						  params);
			
		      if (params->allowReversals)
			hboaRecomputeReversalGains(0,
						   k,
						   numVars,
						   &groupInformation,
						   &operatorGain,
						   &operatorApplicability,
						   G,
						   parents,
						   params);
			
		      if (params->allowJoints)
			hboaRecomputeJointGains(0,
						k,
						numVars,
						&groupInformation,
						&operatorGain,
						&operatorApplicability,
						G,
						parents,
						params);
		   
		    };
		  
		  // swap the graphs
		  
		  swapPointers((void**) &G,(void**) &newG);
		
		  // delete the old one (now located in newG)
		
		  delete newG;

		  break;
		}
	    }
	  else
	    finito = 1;
	};

#ifdef HBOA_PRINT_TOTAL_SCORE_INCREASE

      printf("Total score increase: %f\n",totalScoreIncrease);

#endif
	    
#ifdef HBOA_DISPLAY_MODEL

      printf("Groups:\n");
      for (k=0; k<groupInformation.numGroups; k++)
	{
	  printf("Group %2u: ( %3u",k,groupInformation.groupIndex[k][0]);
	  for (l=1; l<groupInformation.groupSize[k]; l++)
	    printf(" %3u",groupInformation.groupIndex[k][l]);
	  printf(" )\n");
	};

      printf("\nGroup Interactions:\n");
      for (k=0; k<groupInformation.numGroups; k++)
	{
	  printf("%2u <- ",groupInformation.groupIndex[k][0]);
	  for (l=0; l<groupInformation.numGroups; l++)
	    if (G->connected(l,k))
	      printf(" %3lu",l);
	  printf("\n");
	};

#endif

      // compute the conditional probabilities

      for (k=0; k<groupInformation.numGroups; k++)
	{
	  int  tmpIndexSize;
	  int *tmpIndex=NULL;

	  tmpIndexSize = 0;
	  for (m=0; m<G->getNumIn(k); m++)
	    tmpIndexSize += groupInformation.groupSize[G->getParentList(k)[m]];

	  if (tmpIndexSize>0)
	    tmpIndex = (int*) Calloc(tmpIndexSize,sizeof(int));

	  joinGroupIndexes(&groupInformation,
			   G->getParentList(k),
			   G->getNumIn(k),
			   tmpIndex);
	  
	  groupInformation.indexSize[k] = tmpIndexSize+groupInformation.groupSize[k];
	  groupInformation.index[k]     = (int*) Calloc(groupInformation.indexSize[k],sizeof(int));

	  joinIndexes(tmpIndex,
		      tmpIndexSize,
		      groupInformation.groupIndex[k],
		      groupInformation.groupSize[k],
		      groupInformation.index[k]);

	  groupInformation.frequencyTree[k] = new FrequencyTree();

	  groupInformation.frequencyTree[k]->computeIndexedFrequencies(parents,
								       groupInformation.index[k],
								       groupInformation.indexSize[k]);
		    
	  groupInformation.numInstances[k] = groupInformation.frequencyTree[k]->getNumInstances();

	  groupInformation.instances[k]   = (char**)  Calloc(groupInformation.numInstances[k],sizeof(char*));
	  groupInformation.frequencies[k] = (double*) Calloc(groupInformation.numInstances[k],sizeof(double));

	  for (l=0; l<groupInformation.numInstances[k]; l++)
	    groupInformation.instances[k][l] = (char*) Malloc(groupInformation.indexSize[k]);

	  groupInformation.frequencyTree[k]->getInstancesAndFrequencies(groupInformation.instances[k],
									groupInformation.frequencies[k]);

	  groupInformation.parentIndexSize[k] = groupInformation.indexSize[k]-groupInformation.groupSize[k];

	  if (tmpIndex)
	    Free(tmpIndex);
	};

      // topologically order the graph nodes (groups)

      for (k=0; k<groupInformation.numGroups; k++)
	added[k]=0;
	    
      numAdded=0;
	    
      while (numAdded<groupInformation.numGroups)
	{
	  for (k=0; k<groupInformation.numGroups; k++)
	    if (added[k]==0)
	      {
		canAdd = 1;
		for (l=0; l<G->getNumIn(k); l++)
		  if (added[G->getParentList(k)[l]]==0)
		    canAdd = 0;
				
		if (canAdd)
		  {
		    added[k]          = 1;
		    index[numAdded++] = k;
		  };
	      };
	};

#ifdef HBOA_DISPLAY_TOPOLOGICAL_ORDERING

      printf("Topological ordering: ");
      for (k=0; k<groupInformation.numGroups; k++)
	printf("%3u ",index[k]);
      printf("\n");
#endif

    };

  // generate the jth coordinate (chromosome) for all dhildren

  for (i=0; i<M; i++)
    {
      child = &(children->individual[i]);
      x = child->chromosome;
      // 	  x     = getChromosome(child,j);

      // set the already fixed portion

      for (k=0; k<numFixed; k++)
	x[fixedIndex[k]]=fixedValue[k];

      // generate modeled portion
	  
      for (m=0; m<groupInformation.numGroups; m++)
	{
	  int *tmpIndex;
	  int tmpIndexSize;
	  double dAux;
	  double r;
	  int picked;

	  k=index[m];

	  tmpIndex = (int*) Calloc(groupInformation.numInstances[k],sizeof(int));

	  tmpIndexSize=0;

	  // find out what instances we are dealing with (add them to the tmpIndex index)

	  for (l=0; l<groupInformation.numInstances[k]; l++)
	    if (indexedMatch(x,groupInformation.instances[k][l],groupInformation.index[k],groupInformation.parentIndexSize[k]))
	      tmpIndex[tmpIndexSize++] = l;

	  r = drand();

	  dAux=0;
	  for (l=0; (dAux<=r); l++)
	    {
	      dAux += 
		groupInformation.frequencies[k][tmpIndex[l]]/
		groupInformation.frequencyTree[k]->getFrequency(groupInformation.instances[k][tmpIndex[l]],
								groupInformation.parentIndexSize[k]);

	      //		  printf("Yet another frequency: %f\n",groupInformation.frequencies[k][tmpIndex[l]]/
	      //		    groupInformation.frequencyTree[k]->getFrequency(groupInformation.instances[k][tmpIndex[l]],
	      //								    groupInformation.parentIndexSize[k]));
	    };

	  picked = l-1;

	  //	      printf("This frequency should add up to 1. Real value: %f\n",dAux);

	  //	      printf("And the winner is: %u (%f)\n",picked,r);
	      
	  //	      getchar();

	  // well...copy the instance we picked
	      
	  copyIndexed(x,
		      groupInformation.instances[k][tmpIndex[picked]]+groupInformation.parentIndexSize[k],
		      groupInformation.groupIndex[k],
		      groupInformation.groupSize[k]);

	  Free(tmpIndex);

	};

    };

  // free the memory used by instances and frequencies of each group

  for (k=0; k<groupInformation.numGroups; k++)
    {

      for (l=0; l<groupInformation.numInstances[k]; l++)
	Free(groupInformation.instances[k][l]);
      Free(groupInformation.instances[k]);
      Free(groupInformation.frequencies[k]);
      Free(groupInformation.index[k]);
		    
      delete groupInformation.frequencyTree[k];
    };

  // free memory occupied by dependency graph

  if (numVars>0)
    delete G;

  // the fitness has to be evaluated for all children,

  for (i=0; i<M; i++)
    {
      child = &(children->individual[i]);
      child->fCalculated = 0;
    };

  // offspring has not been evaluated yet

  children->evaluated = 0;

  // free the UMFs
  
  freeUMF(&p0,&p1,n  );

  Free(added);
  Free(index);

  for (k=0; k<n; k++)
    {
      Free(operatorApplicability.addition[k]);
      Free(operatorApplicability.removal[k]);
      Free(operatorApplicability.reversal[k]);
      Free(operatorApplicability.joint[k]);
    };
  
  Free(operatorApplicability.addition);
  Free(operatorApplicability.removal);
  Free(operatorApplicability.reversal);
  Free(operatorApplicability.joint);

  for (k=0; k<n; k++)
    {
      Free(operatorGain.addition[k]);
      Free(operatorGain.removal[k]);
      Free(operatorGain.reversal[k]);
      Free(operatorGain.joint[k]);
    };

  Free(operatorGain.addition);
  Free(operatorGain.removal);
  Free(operatorGain.reversal);
  Free(operatorGain.joint);
  
  Free(fixedIndex);
  Free(fixedValue);

 // free the rest of the memory

  for (k=0; k<n; k++)
      Free(groupInformation.groupIndex[k]);

  Free(groupInformation.groupSize);
  Free(groupInformation.groupIndex);
  Free(groupInformation.frequencyTree);
  Free(groupInformation.numInstances);
  Free(groupInformation.instances);
  Free(groupInformation.frequencies);
  Free(groupInformation.full);
  Free(groupInformation.indexSize);
  Free(groupInformation.index);
  Free(groupInformation.parentIndexSize);
  Free(groupInformation.contribution);


  printf("Leaving recombination.\n");

  // get back

  return 0;
};

//=============================================================================================
// checked completely.

int hboaUpdateAfterAddition(int newFrom,
 			    int newTo,
 			    GroupInformation *groupInformation,
 			    OperatorApplicability *operatorApplicability,
 			    AcyclicOrientedGraph *G,
 			    RecombinationParams *params)
{
  int k,l,n;

  n = G->size();

  // the removal is possible when the edge has been just added

  operatorApplicability->removal[newFrom][newTo] = 1;

  // the reversal might be possible, too
  // WARNING: carefully with this.....would get into trouble if updated addition applicability first....
  // could be also ((G->canReverseEdge[newFrom][newTo])&&(!full[newTo])) - safer and slower

  operatorApplicability->reversal[newFrom][newTo] = operatorApplicability->addition[newTo][newFrom];

  // can't add anymore the added neither the opposite edge

  operatorApplicability->addition[newFrom][newTo] = operatorApplicability->addition[newTo][newFrom] = 0;

  // is the ending node of the new edge full yet?

  groupInformation->full[newTo]=(G->getNumIn(newTo)>=params->maxIncoming);

  if (G->getNumIn(newTo)>=3)
    printf("NumIn = %u\n",G->getNumIn(newTo));

  // if the ending node is full now no more edges can end here

  if (groupInformation->full[newTo])
    for (k=0; k<n; k++)
      {
	operatorApplicability->addition[k][newTo] = 0;
	operatorApplicability->reversal[newTo][k] = 0;
      };

  // check for cycles with the new edge...

  for (k=1; k<n; k++)
    for (l=0; l<k; l++)
      if (((k!=newFrom)||(l!=newTo))&&((k!=newTo)||(l!=newFrom)))
	{
	  // does the new edge forbid creating an edge k,l by means of a path that might create a cycle with this?
 	
	  if ((operatorApplicability->addition[k][l])&&(G->existsPath(l,newFrom)&&(G->existsPath(newTo,k))))
	    {
	      operatorApplicability->addition[k][l] = 0;
	      operatorApplicability->reversal[l][k] = 0;
	    };

	  // does the new edge forbid creating an edge l,k by means of a path that might create a cycle with this?
 	
	  if ((operatorApplicability->addition[l][k])&&(G->existsPath(k,newFrom)&&(G->existsPath(newTo,l))))
	    {
	      operatorApplicability->addition[l][k] = 0;
	      operatorApplicability->reversal[k][l] = 0;
	    };
	};

  // that's all, let's get back

  return 0;
};

//==========================================================
// checked completely.

int hboaUpdateAfterRemoval(int newFrom,
 			   int newTo,
 			   GroupInformation *groupInformation,
 			   OperatorApplicability *operatorApplicability,
 			   AcyclicOrientedGraph *G,
 			   RecombinationParams *params)
{
  int k,l,n;

  // initialize the variables

  n = G->size();

  // can't remove or reverse the thing anymore, can add it back though

  operatorApplicability->removal[newFrom][newTo]  = 0;
  operatorApplicability->reversal[newFrom][newTo] = 0;
  operatorApplicability->addition[newFrom][newTo] = 1;

  // can't be full after removing the edge

  groupInformation->full[newTo]=0;

  // allows the edges that got possible to be added or reversed after removal to be added

  for (k=1; k<n; k++)
    for (l=0; l<k; l++)
      {
	// additions go first

	if ((!operatorApplicability->addition[k][l])&&
	    (!groupInformation->full[l])&&
	    (!G->connected(k,l))&&
	    (!G->existsPath(l,k)))
	  operatorApplicability->addition[k][l] = 1;
 	
	if ((!operatorApplicability->addition[l][k])&&
	    (!groupInformation->full[k])&&
	    (!G->connected(l,k))&&
	    (!G->existsPath(k,l)))
	  operatorApplicability->addition[l][k] = 1;

	// reversals follow

	if (params->allowReversals)
	  {
	    if ((G->connected(k,l))&&
		(!operatorApplicability->reversal[k][l])&&
		(!groupInformation->full[k])&&
		(G->canReverseEdge(k,l)))
	      operatorApplicability->reversal[k][l] = 1;
 	
	    if ((G->connected(l,k))&&
		(!operatorApplicability->reversal[l][k])&&
		(!groupInformation->full[l])&&
		(G->canReverseEdge(l,k)))
	      operatorApplicability->reversal[l][k] = 1;
	  };
      };

  return 0;
};

//==========================================================
// checked completely; can be tricky though!!!

int hboaUpdateAfterReversal(int newFrom,
 			    int newTo,
 			    GroupInformation *groupInformation,
 			    OperatorApplicability *operatorApplicability,
 			    AcyclicOrientedGraph *G,
 			    RecombinationParams *params)
{
  //what has happened after we removed an edge newFrom->newTo

  hboaUpdateAfterRemoval(newFrom,newTo,groupInformation,operatorApplicability,G,params);

   // what has happened after we added edge newTo->newFrom

  hboaUpdateAfterAddition(newTo,newFrom,groupInformation,operatorApplicability,G,params);

  // get back

  return 0;
};

//==========================================================

int hboaUpdateAfterJoint(int newFrom,
			 int newTo,
			 GroupInformation *groupInformation,
			 OperatorApplicability *operatorApplicability,
			 AcyclicOrientedGraph *G,
			 RecombinationParams *params)
{
  register int k,l;
  int          n;

  // the size of the graph (# of nodes)

  n = G->size();

  // full update..could be optimized but does it make too much sense at the moment?

  for (k=0; k<n; k++)
    {
      groupInformation->full[k] = (G->getNumIn(k)>=params->maxIncoming);
      
      for (l=0; l<k; l++)
	{
	  operatorApplicability->addition[k][l] = ((!groupInformation->full[k])&&(!G->connected(k,l)))&&(G->canAddEdge(k,l));
	  
	  if (params->allowRemovals)
	    operatorApplicability->removal[k][l]  = (G->connected(k,l));

	  if (params->allowReversals)
	    operatorApplicability->reversal[k][l] = ((G->connected(k,l))&&(G->canReverseEdge(k,l)));
	};
    };

  // get back
  
  return 0;
};


//==========================================================
// checked completely.

int hboaRecomputeAdditionGains(int from,
 			       int to,
			       int numVars,
 			       GroupInformation *groupInformation,
 			       OperatorGain *operatorGain,
 			       OperatorApplicability *operatorApplicability,
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

  n = G->size();
  parentCount = G->getNumIn(to);

  // initialize the parent list

  parentList = G->getParentList(to);

  // allocate memory for the index of nodes from which the edge is to be updated

  updateIdx = (int*) Calloc(n,sizeof(int));

  // create the index of those nodes (!!! might be optimized...10% chance...lol !!!)

  numUpdated = 0;

  if (!params->allowReversals)
    {
      for (k=0; k<n; k++)
	if ((k!=to)&&
	    (operatorApplicability->addition[k][to]))
	  updateIdx[numUpdated++] = k;
    }
  else
    {
      for (k=0; k<n; k++)
	if ((k!=to)&&
	    (!G->connected(k,to)))
	  updateIdx[numUpdated++] = k;
    };

  // compute the gains for edges from all updateIdx to to

  if (numUpdated>0)
    mdlComputeGroupAdditionGains(to,
				 numVars,
				 groupInformation,
				 operatorGain,
				 G,
				 updateIdx,
				 numUpdated,
				 parentList,
				 parentCount,
				 population,
				 params);

  // if the prior network is used, penalize the additions when unmatched with the prior network

  for (k=0; k<numUpdated; k++)
    operatorGain->addition[updateIdx[k]][to] -= getPenalty(updateIdx[k],to);

  // free memory used by list of beginning nodes of the edges to be updated

  Free(updateIdx);

  // get back

  return 0;
};

//==========================================================
// checked completely.

int hboaRecomputeRemovalGains(int from,
 			      int to,
			      int numVars,
 			      GroupInformation *groupInformation,
 			      OperatorGain *operatorGain,
 			      OperatorApplicability *operatorApplicability,
 			      AcyclicOrientedGraph *G,
 			      Population *population,
 			      RecombinationParams *params)
{
  int parentCount;
  int *parentList;
  int k;

  // initialize the variables

  parentCount = G->getNumIn(to);

  // initialize the parent list

  parentList = (int*) Calloc(parentCount,sizeof(int));
  memcpy(parentList,G->getParentList(to),parentCount*sizeof(int));;

  // compute the gains for edges from all updateIdx to to

  if (parentCount>0)
    mdlComputeGroupRemovalGains(to,
				numVars,
				groupInformation,
				operatorGain,
				G,
				parentList,
				parentCount,
				population,
				params);

  // if the prior network is used, penalize the removals when unmatched with the prior network

  for (k=0; k<parentCount; k++)
    operatorGain->removal[parentList[k]][to] += getPenalty(parentList[k],to);

  // free the parent list

  Free(parentList);

  // get back

  return 0;
};

//==========================================================
// checked completely.

int hboaRecomputeReversalGains(int from,
 			       int to,
			       int numVars,
 			       GroupInformation *groupInformation,
 			       OperatorGain *operatorGain,
 			       OperatorApplicability *operatorApplicability,
 			       AcyclicOrientedGraph *G,
 			       Population *population,
 			       RecombinationParams *params)
{
  int k,n;

  // initialize the variables

  n = G->size();

  // compute reversal gains as a function of removal and addition gains (that are all ready waiting)
  // this is a very delicate thing....big deal....not sure it works!!!!!!
  // probably has to compute all addition gains not just legal...

  for (k=0; k<n; k++)
    if ((G->connected(k,to))&&(operatorApplicability->reversal[k][to]))
      operatorGain->reversal[k][to] = operatorGain->removal[k][to]+operatorGain->addition[to][k];

  for (k=0; k<n; k++)
    if ((G->connected(to,k))&&(operatorApplicability->reversal[to][k]))
      operatorGain->reversal[to][k] = operatorGain->removal[to][k]+operatorGain->addition[k][to];

  for (k=0; k<n; k++)
    if ((G->connected(k,from))&&(operatorApplicability->reversal[k][from]))
      operatorGain->reversal[k][from] = operatorGain->removal[k][from]+operatorGain->addition[from][k];

  for (k=0; k<n; k++)
    if ((G->connected(from,k))&&(operatorApplicability->reversal[from][k]))
      operatorGain->reversal[from][k] = operatorGain->removal[from][k]+operatorGain->addition[k][from];

  // get back

  return 0;
};

// ========================================================

int hboaRecomputeJointGains(int from,
			    int to,
			    int numVars,
			    GroupInformation *groupInformation,
			    OperatorGain *operatorGain,
			    OperatorApplicability *operatorApplicability,
			    AcyclicOrientedGraph *G,
			    Population *population,
			    RecombinationParams *params)
{
  fprintf(stderr,"ERROR: Recompute joint gains not implemented yet!\n");
  exit(-1);

  return 0;
};

// ========================================================

double hBoaGetPenalty(int from, int to)
{
  // add penalty function and stuff !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  return 0;
};

// =========================================================

int joinGroups( int a,
 		int b,
 		GroupInformation *groupInformation,
 		OperatorApplicability *operatorApplicability,
 		AcyclicOrientedGraph *G,
 		AcyclicOrientedGraph *newG,
 		Population *P)
{
  int  k,l;
  FrequencyTree *frequencyTree;
  int shiftedA;
  int shiftedK;
  int shiftedL;

  // add all variables from the second group to the first one

  groupInformation->groupSize[a] += groupInformation->groupSize[b];

  for (k=0; k<groupInformation->groupSize[b]; k++)
    groupInformation->groupIndex[a][groupInformation->groupSize[b]+k+1] = groupInformation->groupIndex[b][k];

  // get all instances of a new group that are in the data

  frequencyTree = new FrequencyTree();

  frequencyTree->computeIndexedFrequencies(P,groupInformation->groupIndex[a],groupInformation->groupSize[a]);

  // add all edges that are the same for the both

  shiftedA = (a>=b)? (a-1):a;

  for (k=0; k<G->size(); k++)
    if ((k!=a)&&(k!=b))
      {
	shiftedK = (k>=b)? (k-1):k;

	if ((G->connected(k,a))&&(G->connected(k,b)))
	  newG->addEdge(shiftedK,shiftedA);

	if ((G->connected(a,k))&&(G->connected(b,k)))
	  newG->addEdge(shiftedA,shiftedK);
      };

  // add all remaining edges not interacting with a or b

  for (k=0; k<G->size(); k++)
    if ((k!=a)&&(k!=b))
      for (l=0; l<k; l++)
	if ((l!=a)&&(l!=b))
	  {
	    shiftedK = (k>b)? (k-1):k;
	    shiftedL = (l>b)? (l-1):l;

	    if (G->connected(k,l))
	      newG->addEdge(shiftedK,shiftedL);
	    else
	      if (G->connected(l,k))
		newG->addEdge(shiftedL,shiftedK);
	  }

  // delete the group b

  for (k=b+1; k<groupInformation->numGroups; k++)
    {
      swapInt(&(groupInformation->groupSize[k-1]),&(groupInformation->groupSize[k]));
      swapInt(&(groupInformation->numInstances[k-1]),&(groupInformation->numInstances[k]));
      swapPointers((void**) &(groupInformation->instances[k-1]),(void**) &(groupInformation->instances[k]));
      swapPointers((void**) &(groupInformation->groupIndex[k-1]),(void**) &(groupInformation->groupIndex[k]));
      groupInformation->full[k-1]=groupInformation->full[k];

      for (l=0; l<groupInformation->numGroups; l++)
	{
	  if (l<b)
	    {
	      operatorApplicability->addition[k-1][l] = operatorApplicability->addition[k][l];
	      operatorApplicability->addition[k-1][l] = operatorApplicability->addition[k][l];
	      operatorApplicability->addition[k-1][l] = operatorApplicability->addition[k][l];
	      operatorApplicability->addition[k-1][l] = operatorApplicability->addition[k][l];
	    }
	  else
	    if (l>b)
	      {
		operatorApplicability->addition[k-1][l-1] = operatorApplicability->addition[k][l];
		operatorApplicability->addition[k-1][l-1] = operatorApplicability->addition[k][l];
		operatorApplicability->addition[k-1][l-1] = operatorApplicability->addition[k][l];
		operatorApplicability->addition[k-1][l-1] = operatorApplicability->addition[k][l];
	      };
	};
    };

  groupInformation->numGroups--;

  // get back with the number of a new group

  return a;
};

// =========================================================

int joinGroupIndexes(GroupInformation *groupInformation, int *index, int n, int *destination)
{
  int i,j,k;

  // copy the groups into one array (all variables together into the destination)

  for (i=0,k=0; i<n; i++)
    for (j=0; j<groupInformation->groupSize[index[i]]; j++)
      destination[k++]=groupInformation->groupIndex[index[i]][j];

  // get back

  return 0;
};

// =========================================================

//int joinIndexes(int *index1, int n1, int *index2, int n2, int *destination)
//{
//  int i,k;

  // copy the two arrays into one

//  for (i=0; i<n1; i++)
//    destination[i]=index1[i];
  
//  for (i=0, k=n1; i<n2; i++, k++)
//    destination[k]=index2[i];

  // get back

//  return 0;
//};

// ===================================================================================

int indexedMatch(char *x, char *y, int *index, int n)
{
  for (int i=0; i<n; i++)
    if (x[index[i]]!=y[i])
      return 0;

  return 1;
};

// ===================================================================================

int copyIndexed(char *dest, char *src, int *index, int n)
{
  for (int i=0; i<n; i++)
    dest[index[i]] = src[i];

  return 0;
};
